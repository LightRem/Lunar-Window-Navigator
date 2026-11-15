//! Уточнение найденного перелёта численным интегрированием.
//!
//! Поиск вариантов ([`crate::trajectory`]) идёт по формулам: до Луны аппарат «чувствует» только
//! Землю, у Луны — только Луну. Здесь лучший вариант пересчитывается шаг за шагом со всеми
//! силами сразу ([`crate::propagate`]):
//!
//! 1. **Опорная орбита и перелёт** — Земля со сжатием (J2), Луна, Солнце. Скорость после разгона
//!    к Луне подбирается «методом стрельбы» (метод Ньютона с численными производными) так, чтобы
//!    периселений пришёлся на высоту окололунной орбиты, а плоскость подлёта прошла через цель.
//! 2. **Окололунная орбита** — поле Луны с масконами (GRAIL, [`crate::gravity`]), притяжение
//!    Земли и Солнца. Отсюда фактический момент схода с орбиты и высота в начале торможения.
//! 3. **Торможение двигателем** — моделирование от PDI до остановки над точкой
//!    ([`crate::lunar_landing::simulate_descent`]) с вращением Луны.
//!
//! Плоскость подлёта поправляется, пока точка остановки не совпадёт с точкой посадки.
//! Момент старта и момент разгона берутся из поиска и не меняются.
//!
//! 4. **Разгон к Луне конечной длительности** — если известен двигатель (разгонная ступень ракеты
//!    или сам аппарат), мгновенный импульс заменяется работой двигателя: постоянная тяга, масса
//!    убывает по расходу топлива, направление тяги поворачивается в плоскости опорной орбиты с
//!    постоянной угловой скоростью (скорость поворота подбирается на наименьший расход). Разгон
//!    центрирован на моменте импульса и подбирается так, чтобы прилёт к Луне был тем же
//!    (тот же момент, высота и плоскость периселения). Разница Δv — гравитационные потери.

use anise::constants::frames::{MOON_J2000, MOON_ME_DE440_ME421_FRAME};
use anise::math::Vector3;
use anise::prelude::Almanac;
use anise::time::{Epoch, Unit};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};
use crate::gravity::GravityField;
use crate::lambert;
use crate::lunar_landing::{
    DescentSim, LandingTarget, MOON_ROTATION_RATE_RAD_S, REASON_LANDING_BOTH, REASON_LANDING_EARTH,
    REASON_LANDING_SUN, REASON_LANDING_WAIT, rotate_about, simulate_descent,
};
use crate::propagate::{
    Dynamics, EarthCentered, Forces, MoonCentered, State, integrate, pos, state, vel,
};
use crate::runlog::{ser_epoch, ser_vec3};
use crate::sky::{MOON_RADIUS_KM, sky_at};
use crate::trajectory::{BudgetCheck, DeltaVBudget, Maneuver, TransferSolution};
use crate::trajlog::TrajectoryPath;
use crate::visibility::BodyVisibility;

/// Параметры численного уточнения (секция `[numerical]` в `data/transfer.toml`).
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct NumericalConfig {
    /// Уточнять найденный перелёт численным интегрированием.
    pub enabled: bool,
    /// Учитывать сжатие Земли (J2).
    pub earth_j2: bool,
    /// Учитывать притяжение Солнца.
    pub sun_gravity: bool,
    /// Степень разложения поля Луны (масконы): 0 — точечная масса, наибольшая — 60.
    pub lunar_gravity_degree: usize,
    /// Относительная точность шага интегрирования.
    pub tolerance: f64,
    /// Допустимый промах точки остановки относительно точки посадки, км.
    pub landing_tolerance_km: f64,
    /// Наибольшее число итераций подбора.
    pub max_iterations: u32,
    /// Уточнять лучший вариант каждого космодрома (`true`) или только лучший из всех (`false`).
    pub refine_all: bool,
}

impl Default for NumericalConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            earth_j2: true,
            sun_gravity: true,
            lunar_gravity_degree: 60,
            tolerance: 1e-10,
            landing_tolerance_km: 0.1,
            max_iterations: 15,
            refine_all: false,
        }
    }
}

impl NumericalConfig {
    pub fn validate(&self) -> Result<()> {
        let bad = |m: &str| Err(Error::Config(format!("numerical: {m}")));
        if self.lunar_gravity_degree > 60 {
            return bad("lunar_gravity_degree must be in [0, 60]");
        }
        if !(self.tolerance >= 1e-13 && self.tolerance <= 1e-6) {
            return bad("tolerance must be in [1e-13, 1e-6]");
        }
        if !(self.landing_tolerance_km >= 1e-3 && self.landing_tolerance_km <= 100.0) {
            return bad("landing_tolerance_km must be in [0.001, 100]");
        }
        if !(1..=100).contains(&self.max_iterations) {
            return bad("max_iterations must be in [1, 100]");
        }
        Ok(())
    }
}

/// Уточнённая посадка.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct RefinedLanding {
    /// Наклонение окололунной орбиты к лунному экватору при выходе на неё, градусы.
    pub orbit_inclination_deg: f64,
    /// Наименьшая и наибольшая высота над сферой 1737.4 км за время ожидания (масконы), км.
    pub orbit_min_altitude_km: f64,
    pub orbit_max_altitude_km: f64,
    #[serde(serialize_with = "ser_epoch")]
    pub doi_epoch: Epoch,
    pub doi_delta_v_m_s: f64,
    /// Время на окололунной орбите от выхода на неё до схода, с.
    pub wait_in_orbit_s: f64,
    /// На сколько витков позже, чем по формулам, выполнен спуск (если на исходном витке не выполнялись условия).
    pub extra_orbits: u32,
    #[serde(serialize_with = "ser_epoch")]
    pub pdi_epoch: Epoch,
    /// Фактическая высота в начале торможения над сферой 1737.4 км, км.
    pub pdi_altitude_km: f64,
    /// Высота в начале торможения над точкой посадки, км.
    pub pdi_height_above_site_km: f64,
    pub pdi_speed_km_s: f64,
    /// Торможение двигателем до остановки над точкой.
    pub braking: DescentSim,
    /// Заход и вертикальный участок после торможения, с.
    pub terminal_duration_s: f64,
    /// Момент касания.
    #[serde(serialize_with = "ser_epoch")]
    pub touchdown_epoch: Epoch,
    /// Сдвиг момента касания относительно расчёта по формулам, с.
    pub touchdown_shift_s: f64,
    /// Весь спуск: сход с орбиты + торможение + добавка на заход, м/с.
    pub descent_delta_v_m_s: f64,
    /// Промах точки остановки относительно точки посадки, км.
    pub miss_km: f64,
    pub sun: BodyVisibility,
    pub earth: BodyVisibility,
    /// Условия посадки (Солнце, Земля) выполнены в уточнённый момент касания.
    pub conditions_met: bool,
    pub unmet_reason: Option<&'static str>,
}

/// Результат численного уточнения.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Refined {
    pub config: NumericalConfig,
    /// Подбор сошёлся в пределах допусков.
    pub converged: bool,
    pub iterations: u32,
    /// Шагов интегрирования всего.
    pub integration_steps: u64,
    /// Сдвиг точки разгона из-за сжатия Земли и притяжения Луны и Солнца на опорной орбите, км.
    pub parking_position_shift_km: f64,
    #[serde(serialize_with = "ser_vec3")]
    pub tli_position_km: Vector3,
    #[serde(serialize_with = "ser_vec3")]
    pub tli_velocity_after_km_s: Vector3,
    pub tli_delta_v_km_s: f64,
    /// Изменение Δv разгона относительно расчёта по формулам, м/с.
    pub tli_delta_v_change_m_s: f64,
    /// Момент периселения (выхода на окололунную орбиту).
    #[serde(serialize_with = "ser_epoch")]
    pub arrival_epoch: Epoch,
    /// Сдвиг момента прилёта относительно расчёта по формулам, с.
    pub arrival_shift_s: f64,
    /// Удалось ли сохранить момент прилёта из поиска (иначе он сдвинулся из-за притяжения Луны).
    pub arrival_time_kept: bool,
    pub time_of_flight_s: f64,
    pub periselene_altitude_km: f64,
    pub periselene_speed_km_s: f64,
    pub loi_delta_v_km_s: f64,
    pub loi_delta_v_change_m_s: f64,
    pub landing: Option<RefinedLanding>,
    /// Бюджет аппарата по уточнённым манёврам.
    pub budget: Option<BudgetCheck>,
    /// Разгон к Луне конечной длительности (если известен двигатель); иначе разгон — мгновенный импульс.
    pub tli_burn: Option<TliBurn>,
    /// Гравитационный параметр Луны в модели поля, км³/с².
    pub moon_gm_km3_s2: f64,
    /// Точки и события траектории — для журнала траектории (в журнал информации не пишутся).
    #[serde(skip)]
    pub path: Option<TrajectoryPath>,
}

/// Ускорение свободного падения для удельного импульса, м/с².
const G0: f64 = 9.806_65;

/// Двигатель, выполняющий разгон к Луне.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Engine {
    /// Чем выполняется разгон (например, «Fregat (Soyuz-2.1b / Fregat)»).
    pub name: String,
    pub thrust_n: f64,
    pub isp_s: f64,
    /// Масса в начале разгона (ступень с остатком топлива и груз), кг.
    pub mass_kg: f64,
    /// Топливо, доступное на разгон, кг.
    pub propellant_kg: f64,
}

/// Разгон к Луне конечной длительности (результат уточнения).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TliBurn {
    pub engine: Engine,
    #[serde(serialize_with = "ser_epoch")]
    pub start_epoch: Epoch,
    #[serde(serialize_with = "ser_epoch")]
    pub end_epoch: Epoch,
    pub duration_s: f64,
    /// Δv, который даёт двигатель (по формуле Циолковского), км/с.
    pub delta_v_km_s: f64,
    /// Δv мгновенного импульса для того же прилёта, км/с.
    pub impulsive_delta_v_km_s: f64,
    /// Гравитационные потери: Δv двигателя минус Δv мгновенного импульса, м/с.
    pub gravity_loss_m_s: f64,
    pub propellant_used_kg: f64,
    /// Наибольший Δv, который двигатель может дать при доступном топливе, км/с.
    pub max_delta_v_km_s: f64,
    /// Топлива на этот разгон не хватает.
    pub exceeds_propellant: bool,
    /// Топлива не хватило бы и на мгновенный импульс (тогда это не потери, а допущение о
    /// начальной орбите: паспортная масса к Луне подразумевает разгон с более высокой орбиты).
    pub impulsive_exceeds_propellant: bool,
    /// Скорость поворота направления тяги в плоскости орбиты, доля от угловой скорости движения по орбите.
    pub steering_rate_fraction: f64,
    /// Угол, пройденный по орбите за время разгона, градусы.
    pub burn_arc_deg: f64,
    #[serde(serialize_with = "ser_vec3")]
    pub end_position_km: Vector3,
    #[serde(serialize_with = "ser_vec3")]
    pub end_velocity_km_s: Vector3,
}

/// Модель работы двигателя на разгоне.
#[derive(Debug, Clone, Copy)]
pub(crate) struct BurnModel {
    thrust_n: f64,
    /// Скорость истечения, м/с.
    c_m_s: f64,
    m0: f64,
    /// Ось поворота направления тяги (нормаль опорной орбиты).
    normal: Vector3,
    /// Скорость поворота направления тяги, рад/с.
    omega: f64,
}

impl BurnModel {
    fn duration(&self, dv_km_s: f64) -> f64 {
        self.propellant(dv_km_s) * self.c_m_s / self.thrust_n
    }

    fn propellant(&self, dv_km_s: f64) -> f64 {
        self.m0 * (1.0 - (-dv_km_s * 1000.0 / self.c_m_s).exp())
    }

    /// Начало и конец разгона Δv, центрированного на `t_mid` (но не раньше `t_min`).
    fn window(&self, dv_km_s: f64, t_mid: f64, t_min: f64) -> (f64, f64) {
        let tb = self.duration(dv_km_s);
        let t_s = (t_mid - 0.5 * tb).max(t_min);
        (t_s, t_s + tb)
    }
}

/// Полёт с работающим двигателем: притяжение как у [`EarthCentered`] плюс тяга.
pub(crate) struct Thrusting<'a> {
    base: &'a EarthCentered<'a>,
    m: BurnModel,
    t_s: f64,
    t_mid: f64,
    dir: Vector3,
}

impl Dynamics for Thrusting<'_> {
    fn acceleration(&self, t: f64, y: &State) -> Result<Vector3> {
        let mass = self.m.m0 - self.m.thrust_n / self.m.c_m_s * (t - self.t_s);
        let d = rotate_about(self.dir, self.m.normal, self.m.omega * (t - self.t_mid));
        // Н/кг = м/с², в км/с² — делить на 1000.
        Ok(self.base.acceleration(t, y)? + d * (self.m.thrust_n / mass / 1000.0))
    }
}

/// Разгон: от точки `anchor` на опорной орбите (раньше начала разгона) до начала разгона,
/// затем работа двигателя. Направление тяги в середине разгона — по `dv`, |dv| — Δv двигателя.
/// Возвращает начало, конец, состояние в конце и число шагов.
fn burn_through(
    base: &EarthCentered,
    anchor: (f64, State),
    m: &BurnModel,
    t_mid: f64,
    dv: &Vector3,
    tol: f64,
) -> Result<(f64, f64, State, u32)> {
    let (t_s, t_e) = m.window(dv.norm(), t_mid, anchor.0);
    let (y_s, n1) = if t_s > anchor.0 {
        let o = integrate(base, anchor.0, &anchor.1, t_s, tol, None, None)?;
        (o.y, o.steps)
    } else {
        (anchor.1, 0)
    };
    let th = Thrusting {
        base,
        m: *m,
        t_s,
        t_mid: 0.5 * (t_s + t_e),
        dir: dv.normalize(),
    };
    let o = integrate(&th, t_s, &y_s, t_e, tol, None, None)?;
    Ok((t_s, t_e, o.y, n1 + o.steps))
}

/// Начальное приближение для конечного разгона: Δv (вектор в середине разгона), при котором
/// орбита сразу после работы двигателя совпадает с орбитой после мгновенного импульса `dv_imp`
/// в `t_mid` — та же энергия и то же направление перигея. Без этого на медленных перелётах
/// разгон «тех же» Δv уводит точку прилёта на десятки тысяч километров, и подбор не сходится.
fn match_post_burn_orbit(
    base: &EarthCentered,
    anchor: (f64, State),
    m: &BurnModel,
    t_mid: f64,
    r1: &Vector3,
    v_park: &Vector3,
    dv_imp: &Vector3,
    tol: f64,
) -> Result<Vector3> {
    let mu = base.mu_earth;
    let elements = |r: &Vector3, v: &Vector3| {
        let energy = 0.5 * v.norm_squared() - mu / r.norm();
        let e_vec = v.cross(&r.cross(v)) / mu - r.normalize();
        (energy, e_vec)
    };
    let (e_goal, ecc_goal) = elements(r1, &(v_park + dv_imp));
    let normal = m.normal;
    let mut dv = *dv_imp;
    for _ in 0..8 {
        let (_, _, y, _) = burn_through(base, anchor, m, t_mid, &dv, tol)?;
        let (r, v) = (pos(&y), vel(&y));
        let (e_now, ecc_now) = elements(&r, &v);
        // dE ≈ v·dv: поправка величины; поворот — на угол между направлениями перигея.
        let d_mag = (e_goal - e_now) / v.norm();
        let (a, b) = (
            ecc_now - normal * normal.dot(&ecc_now),
            ecc_goal - normal * normal.dot(&ecc_goal),
        );
        let angle = a.cross(&b).dot(&normal).atan2(a.dot(&b));
        dv = rotate_about(dv, normal, angle);
        dv += dv.normalize() * d_mag;
        if d_mag.abs() < 1e-7 && angle.abs() < 1e-9 {
            break;
        }
    }
    Ok(dv)
}

#[derive(Clone, Copy)]
struct Peri {
    t: f64,
    r: Vector3,
    v: Vector3,
    steps: u32,
}

struct Shooter<'a> {
    dynamics: &'a EarthCentered<'a>,
    tol: f64,
    t_tli: f64,
    t_max: f64,
    r1: Vector3,
    r_target: f64,
    /// Заданный момент периселения, с от опорной эпохи.
    t_target: f64,
    /// Закрепить момент периселения (иначе он свободен).
    fix_time: bool,
    /// Скорость на опорной орбите в момент разгона.
    v_park: Vector3,
    /// Δv разгона по формулам, км/с.
    dv_ref: f64,
    /// Режим грубого наведения: нормаль желаемой плоскости подлёта.
    coarse_plane: Option<Vector3>,
    /// Выбор «закреплять ли момент прилёта» уже сделан.
    decided: bool,
    /// Разгон конечной длительности: модель двигателя и точка опоры на опорной орбите.
    /// `None` — мгновенный импульс в `t_tli`.
    burn: Option<(BurnModel, (f64, State))>,
}

/// Решает A·x = b (A — m×m, m ≤ 3) методом Гаусса с выбором главного элемента.
fn solve_small(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let m = b.len();
    for i in 0..m {
        let piv = (i..m).max_by(|&x, &y| a[x][i].abs().total_cmp(&a[y][i].abs()))?;
        if a[piv][i].abs() < 1e-300 {
            return None;
        }
        a.swap(i, piv);
        b.swap(i, piv);
        for k in (i + 1)..m {
            let f = a[k][i] / a[i][i];
            for j in i..m {
                a[k][j] -= f * a[i][j];
            }
            b[k] -= f * b[i];
        }
    }
    let mut x = vec![0.0; m];
    for i in (0..m).rev() {
        let sum: f64 = ((i + 1)..m).map(|j| a[i][j] * x[j]).sum();
        x[i] = (b[i] - sum) / a[i][i];
    }
    Some(x)
}

impl Shooter<'_> {
    /// Интегрирует перелёт со скоростью `u` после разгона до периселения.
    fn fly(&self, u: &Vector3) -> Result<Option<Peri>> {
        let d = self.dynamics;
        let mut near_moon = |t: f64, y: &State| -> Result<f64> {
            let (rm, vm) = d.moon(t)?;
            let (r, v) = (pos(y) - rm, vel(y) - vm);
            // Ближайшая к Луне точка; у Земли (до 200 тыс. км) не ищем.
            Ok(if pos(y).norm() < 200_000.0 {
                1.0
            } else {
                -r.dot(&v) / (r.norm() * v.norm())
            })
        };
        let (t_start, y_start, burn_steps) = match &self.burn {
            Some((m, anchor)) => {
                let (_, t_e, y, n) =
                    burn_through(d, *anchor, m, self.t_tli, &(u - self.v_park), self.tol)?;
                (t_e, y, n)
            }
            None => (self.t_tli, state(&self.r1, u), 0),
        };
        let out = integrate(
            d,
            t_start,
            &y_start,
            self.t_max,
            self.tol,
            Some(&mut near_moon),
            None,
        )?;
        if !out.event {
            return Ok(None);
        }
        let (rm, vm) = d.moon(out.t)?;
        Ok(Some(Peri {
            t: out.t,
            r: pos(&out.y) - rm,
            v: vel(&out.y) - vm,
            steps: out.steps + burn_steps,
        }))
    }

    /// Невязки, км:
    /// 1. прицельная дальность |B| минус та, при которой периселений на нужной высоте
    ///    (|B| почти линейно зависит от скорости после разгона, в отличие от высоты периселения);
    /// 2. если есть цель `aim` — выход плоскости подлёта из цели;
    /// 3. если момент прилёта закреплён — отклонение момента периселения (1 с считается за 1 км).
    fn residuals(&self, p: &Peri, aim: Option<&Vector3>) -> Vec<f64> {
        let mu = self.dynamics.mu_moon;
        let v_inf2 = (p.v.norm_squared() - 2.0 * mu / p.r.norm()).max(1e-6);
        let b_now = p.r.cross(&p.v).norm() / v_inf2.sqrt();
        let b_goal = self.r_target * (1.0 + 2.0 * mu / (self.r_target * v_inf2)).sqrt();
        if let Some(h_des) = self.coarse_plane {
            // Грубое наведение: точка прицеливания B в плоскости поперёк направления подлёта.
            // B зависит от скорости после разгона почти линейно, поэтому подбор сходится издалека,
            // даже если начальная траектория проходит у самого центра Луны.
            let (r_hat, v_hat) = (p.r.normalize(), p.v.normalize());
            let ecc = 1.0 + p.r.norm() * v_inf2 / mu;
            let s_hat = r_hat / ecc + v_hat * ((ecc * ecc - 1.0).sqrt() / ecc);
            let b_vec = s_hat.cross(&r_hat.cross(&v_hat)) * b_now;
            let t_axis = s_hat.cross(&h_des).normalize();
            let r_axis = s_hat.cross(&t_axis);
            return vec![b_vec.dot(&t_axis) - b_goal, b_vec.dot(&r_axis)];
        }
        let mut f = vec![b_now - b_goal];
        if let Some(s) = aim {
            f.push(p.r.cross(&p.v).normalize().dot(s) * self.r_target);
        }
        if self.fix_time {
            f.push(p.t - self.t_target);
        }
        f
    }

    /// Подбирает скорость после разгона. Возвращает скорость, периселений, число итераций,
    /// признак сходимости и число шагов интегрирования.
    fn solve(
        &self,
        mut u: Vector3,
        aim: Option<&Vector3>,
        max_iter: u32,
    ) -> Result<std::result::Result<Shot, String>> {
        const MISS: &str = "the integrated trajectory does not reach the Moon";
        let Some(mut p) = self.fly(&u)? else {
            return Ok(Err(MISS.into()));
        };
        let mut steps = p.steps as u64;
        for it in 0..max_iter {
            let f = self.residuals(&p, aim);
            if f.iter().all(|x| x.abs() < 1e-2) {
                return Ok(Ok((u, p, it, true, steps)));
            }
            // Численные производные по трём составляющим скорости (1 мм/с), параллельно.
            let du = 1e-6;
            let cols: Vec<Option<(Vec<f64>, u32)>> = (0..3)
                .into_par_iter()
                .map(|i| {
                    let mut e = Vector3::zeros();
                    e[i] = du;
                    Ok(self
                        .fly(&(u + e))?
                        .map(|q| (self.residuals(&q, aim), q.steps)))
                })
                .collect::<Result<_>>()?;
            let mut jac = vec![Vector3::zeros(); f.len()];
            for (i, col) in cols.iter().enumerate() {
                let Some((fi, n)) = col else {
                    return Ok(Err(MISS.into()));
                };
                steps += *n as u64;
                for (row, (a, b)) in jac.iter_mut().zip(fi.iter().zip(&f)) {
                    row[i] = (a - b) / du;
                }
            }
            // Шаг наименьшей нормы: Δu = −Jᵀ (J Jᵀ)⁻¹ f.
            let gram: Vec<Vec<f64>> = jac
                .iter()
                .map(|a| jac.iter().map(|b| a.dot(b)).collect())
                .collect();
            let Some(lambda) = solve_small(gram, f.clone()) else {
                return Ok(Err(
                    "targeting is degenerate (the landing site lies on the approach line)".into(),
                ));
            };
            let mut step = Vector3::zeros();
            for (row, l) in jac.iter().zip(&lambda) {
                step -= row * *l;
            }
            // Полный шаг Ньютона принимается, если невязка уменьшилась; иначе шаг укорачивается.
            let merit = |f: &[f64]| f.iter().map(|x| x * x).sum::<f64>();
            let m0 = merit(&f);
            let mut accepted = None;
            for _ in 0..12 {
                if let Some(q) = self.fly(&(u + step))? {
                    steps += q.steps as u64;
                    if merit(&self.residuals(&q, aim)) < m0 {
                        accepted = Some(q);
                        break;
                    }
                }
                step *= 0.5;
            }
            let Some(q) = accepted else { break };
            u += step;
            p = q;
        }
        let ok = self.residuals(&p, aim).iter().all(|x| x.abs() < 1e-2);
        Ok(Ok((u, p, max_iter, ok, steps)))
    }
}

type Shot = (Vector3, Peri, u32, bool, u64);

/// Первый подбор — с моментом прилёта из поиска. У части перелётов притяжение Луны заметно
/// ускоряет подлёт, и прежний момент либо недостижим, либо требует лишнего топлива — тогда
/// момент отпускается (аппарат прилетает в другое время). Выбор сохраняется.
fn solve_auto(
    sh: &mut Shooter,
    u: Vector3,
    aim: Option<&Vector3>,
    max_iter: u32,
) -> Result<std::result::Result<Shot, String>> {
    if sh.decided {
        return sh.solve(u, aim, max_iter);
    }
    sh.decided = true;
    sh.fix_time = true;
    // Момент прилёта из поиска сохраняется, если это стоит не больше 30 м/с разгона сверх расчёта.
    if let Ok(fixed) = sh.solve(u, aim, max_iter)? {
        if fixed.3 && (fixed.0 - sh.v_park).norm() - sh.dv_ref <= 0.030 {
            return Ok(Ok(fixed));
        }
    }
    sh.fix_time = false;
    sh.solve(u, aim, max_iter)
}

/// Угол по ходу движения от аппарата до направления `s`, отсчитанный в плоскости орбиты, рад.
fn angle_ahead(y: &State, s: &Vector3) -> f64 {
    let (r, v) = (pos(y), vel(y));
    let h = r.cross(&v).normalize();
    let r_hat = r.normalize();
    r_hat.cross(s).dot(&h).atan2(r_hat.dot(s))
}

struct LunarOut {
    /// Состояние до и после импульса схода с орбиты (относительно Луны).
    doi_pre: State,
    doi_post: State,
    inclination_deg: f64,
    min_alt: f64,
    max_alt: f64,
    t_doi: f64,
    doi_dv_m_s: f64,
    t_pdi: f64,
    pdi: State,
    sim: DescentSim,
    /// Момент остановки над точкой, с.
    t_stop: f64,
    /// Ошибка: направление на точку посадки минус направление на точку остановки.
    error: Vector3,
    /// Боковая составляющая ошибки (вдоль нормали к орбите), рад.
    cross: f64,
    steps: u64,
}

struct LunarPhase<'a> {
    almanac: &'a Almanac,
    dynamics: &'a MoonCentered<'a>,
    target: &'a LandingTarget,
    tol: f64,
    mu: f64,
    /// Δv, израсходованный аппаратом до схода с орбиты, м/с.
    spent_m_s: f64,
}

impl LunarPhase<'_> {
    fn site_dir(&self, t: f64) -> Result<Vector3> {
        let e = self.dynamics.epoch(t);
        let dcm = self
            .almanac
            .rotate(MOON_ME_DE440_ME421_FRAME, MOON_J2000, e)
            .map_err(|err| Error::Almanac(format!("Moon orientation at {e}: {err}")))?;
        Ok(dcm * self.target.site.site.up())
    }

    /// Полёт от выхода на орбиту (`y0` при t = 0) до остановки над поверхностью.
    /// `t_doi_min` — раньше этого момента сход с орбиты не выполняется.
    fn run(
        &self,
        y0: &State,
        t_doi_min: f64,
    ) -> Result<std::result::Result<LunarOut, &'static str>> {
        let d = self.dynamics;
        let site_radius = self.target.site.site.position_km().norm();
        let model = &self.target.descent;
        let r_orbit = pos(y0).norm();
        let plan = match model.plan(self.mu, r_orbit, self.spent_m_s, site_radius)? {
            Ok(p) => p,
            Err(reason) => return Ok(Err(reason)),
        };
        let period = std::f64::consts::TAU * (r_orbit.powi(3) / self.mu).sqrt();
        let pole = {
            let e = d.epoch(0.0);
            self.almanac
                .rotate(MOON_ME_DE440_ME421_FRAME, MOON_J2000, e)
                .map_err(|err| Error::Almanac(format!("Moon orientation at {e}: {err}")))?
                * Vector3::new(0.0, 0.0, 1.0)
        };
        let inclination_deg = pos(y0)
            .cross(&vel(y0))
            .normalize()
            .dot(&pole)
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees();
        let (mut min_alt, mut max_alt) = (f64::INFINITY, f64::NEG_INFINITY);
        let mut track = |_t: f64, y: &State| {
            let alt = pos(y).norm() - MOON_RADIUS_KM;
            min_alt = min_alt.min(alt);
            max_alt = max_alt.max(alt);
        };
        let mut steps = 0u64;

        // Ожидание на орбите, затем сход за (180° + дуга торможения) до точки посадки.
        let wait = integrate(d, 0.0, y0, t_doi_min, self.tol, None, Some(&mut track))?;
        steps += wait.steps as u64;
        let lead = plan.total_s();
        let mut rho = plan.range_angle_rad;
        let mut doi_event = |t: f64, y: &State| -> Result<f64> {
            Ok((angle_ahead(y, &self.site_dir(t + lead)?) - std::f64::consts::PI - rho).sin())
        };
        let doi = integrate(
            d,
            wait.t,
            &wait.y,
            wait.t + 1.6 * period,
            self.tol,
            Some(&mut doi_event),
            Some(&mut track),
        )?;
        steps += doi.steps as u64;
        if !doi.event {
            return Err(Error::Config(
                "numerical refinement: no descent opportunity found within one orbit".into(),
            ));
        }
        // Импульс схода: точка схода становится апоцентром эллипса с перицентром на высоте PDI.
        let (r, v) = (pos(&doi.y), vel(&doi.y));
        let r_peri = MOON_RADIUS_KM + model.config.pdi_altitude_km;
        let sma = 0.5 * (r.norm() + r_peri);
        let speed = (self.mu * (2.0 / r.norm() - 1.0 / sma)).sqrt();
        let v_new = r.cross(&v).normalize().cross(&r.normalize()) * speed;
        let doi_dv_m_s = (v - v_new).norm() * 1000.0;
        let mass = model.mass_after(self.spent_m_s + doi_dv_m_s);
        let (vehicle, _) = model.vehicle();

        // Первые 10 минут эллипса — общие для всех итераций подбора момента PDI.
        let first = integrate(
            d,
            doi.t,
            &state(&r, &v_new),
            doi.t + 600.0,
            self.tol,
            None,
            None,
        )?;
        steps += first.steps as u64;
        let omega = pole * MOON_ROTATION_RATE_RAD_S;
        let mut braking_s = plan.braking_duration_s;
        let mut result = None;
        for _ in 0..8 {
            let mut pdi_event = |t: f64, y: &State| -> Result<f64> {
                Ok((angle_ahead(y, &self.site_dir(t + braking_s)?) - rho).sin())
            };
            let pdi = integrate(
                d,
                first.t,
                &first.y,
                first.t + period,
                self.tol,
                Some(&mut pdi_event),
                None,
            )?;
            steps += pdi.steps as u64;
            if !pdi.event {
                return Err(Error::Config(
                    "numerical refinement: descent ellipse never reaches the braking point".into(),
                ));
            }
            let sim = match simulate_descent(
                &pos(&pdi.y),
                &vel(&pdi.y),
                mass,
                &vehicle,
                model.config.braking_throttle,
                &omega,
                self.mu,
                site_radius,
            ) {
                Ok(s) => s,
                Err(reason) => return Ok(Err(reason)),
            };
            let t_stop = pdi.t + sim.duration_s;
            let s = self.site_dir(t_stop)?;
            let p_hat = sim.touchdown_km.normalize();
            let h = pos(&pdi.y).cross(&vel(&pdi.y)).normalize();
            // Недолёт (+) или перелёт (−) вдоль трассы: сдвигаем начало торможения.
            let along = p_hat.cross(&s).dot(&h).atan2(p_hat.dot(&s));
            braking_s = sim.duration_s;
            let done = along.abs() * site_radius < 0.02;
            result = Some(LunarOut {
                doi_pre: doi.y,
                doi_post: state(&r, &v_new),
                inclination_deg,
                min_alt: 0.0,
                max_alt: 0.0,
                t_doi: doi.t,
                doi_dv_m_s,
                t_pdi: pdi.t,
                pdi: pdi.y,
                sim,
                t_stop,
                error: s - p_hat,
                cross: (s - p_hat).dot(&h),
                steps: 0,
            });
            if done {
                break;
            }
            rho -= along;
        }
        let mut out = result.expect("at least one iteration");
        (out.min_alt, out.max_alt, out.steps) = (min_alt, max_alt, steps);
        Ok(Ok(out))
    }
}

/// Уточняет найденный перелёт. `Ok(Err(причина))` — уточнить не удалось.
pub fn refine(
    almanac: &Almanac,
    cfg: &NumericalConfig,
    sol: &TransferSolution,
    landing: Option<&LandingTarget>,
    budget: Option<&DeltaVBudget>,
    tli_engine: Option<&Engine>,
) -> Result<std::result::Result<Refined, String>> {
    cfg.validate()?;
    let forces = Forces {
        earth_j2: cfg.earth_j2,
        sun: cfg.sun_gravity,
    };
    let (mu_e, mu_m_conic) = (sol.mu_earth_km3_s2, sol.mu_moon_km3_s2);
    let field = GravityField::moon_grail(cfg.lunar_gravity_degree)?;
    let mu_m = field.gm;
    let t0 = sol.launch.epoch;
    let earth = EarthCentered::new(almanac, t0, mu_e, mu_m_conic, forces)?;
    let dep = &sol.departure;
    let arr = &sol.arrival;
    let r_llo = arr.lunar_orbit_radius_km;

    // Опорная орбита: от выведения до разгона.
    // Опорная орбита: от выхода на неё (после подъёма) до разгона.
    let h_park = sol.parking_orbit.plane_normal;
    let y_ins = state(
        &sol.ascent.insertion_position_km,
        &sol.ascent.insertion_velocity_km_s,
    );
    let t_ins = sol.ascent.duration_s;
    let t_tli = (dep.epoch - t0).to_seconds();
    let park = integrate(&earth, t_ins, &y_ins, t_tli, cfg.tolerance, None, None)?;
    let (r1, v_park) = (pos(&park.y), vel(&park.y));
    let mut steps_total = park.steps as u64;

    // Начальное приближение: прицел не в центр Луны, а мимо него — на нужную высоту и плоскость.
    let tof = (arr.epoch - dep.epoch).to_seconds();
    let v_inf = arr.v_infinity_km_s;
    let v_hat = v_inf.normalize();
    let pole = almanac
        .rotate(MOON_ME_DE440_ME421_FRAME, MOON_J2000, arr.epoch)
        .map_err(|e| Error::Almanac(format!("Moon orientation at {}: {e}", arr.epoch)))?
        * Vector3::new(0.0, 0.0, 1.0);
    let conic_landing = sol.landing.as_ref();
    let h_lunar = match conic_landing {
        Some(l) => l.lunar_orbit_normal_j2000,
        None => (pole - v_hat * v_hat.dot(&pole)).normalize(),
    };
    let b = r_llo * (1.0 + 2.0 * mu_m / (r_llo * v_inf.norm_squared())).sqrt();
    let aim_point = arr.moon_position_km + v_hat.cross(&h_lunar) * b;
    // Если пересчёт начального приближения неустойчив (перелёт почти на 180°), берётся скорость из поиска.
    let u_conic = dep.velocity_transfer_km_s + (v_park - dep.velocity_parking_km_s);
    let u0 = lambert::solve(&r1, &aim_point, tof, mu_e, &h_park)
        .map(|s| s.v1)
        .ok()
        .filter(|v| (v - u_conic).norm() < 0.05)
        .unwrap_or(u_conic);

    let mut shooter = Shooter {
        dynamics: &earth,
        tol: cfg.tolerance,
        t_tli,
        t_max: t_tli + 2.0 * tof,
        r1,
        r_target: r_llo,
        t_target: t_tli + tof,
        fix_time: true,
        v_park,
        dv_ref: dep.delta_v_magnitude_km_s,
        coarse_plane: None,
        decided: false,
        burn: None,
    };
    let site_dir_at = |e: Epoch| -> Result<Vector3> {
        let dcm = almanac
            .rotate(MOON_ME_DE440_ME421_FRAME, MOON_J2000, e)
            .map_err(|err| Error::Almanac(format!("Moon orientation at {e}: {err}")))?;
        Ok(dcm * landing.expect("landing target").site.site.up())
    };

    let spent_before_descent = |dv_tli: f64, dv_loi: f64| -> f64 {
        budget.map_or(0.0, |b| {
            b.maneuvers
                .iter()
                .map(|m| match m {
                    Maneuver::Tli => dv_tli * 1000.0,
                    Maneuver::Loi => dv_loi * 1000.0,
                    Maneuver::Descent => 0.0,
                })
                .sum()
        })
    };

    // Грубое наведение на нужную высоту и плоскость, затем точное.
    shooter.coarse_plane = Some(h_lunar);
    let (mut u, mut iterations) = match shooter.solve(u0, None, cfg.max_iterations)? {
        Ok((u, _, it, _, n)) => {
            steps_total += n;
            (u, it)
        }
        Err(reason) => return Ok(Err(reason)),
    };
    shooter.coarse_plane = None;
    let converged;
    let peri;
    let mut refined_landing: Option<RefinedLanding> = None;
    let mut lunar_out: Option<LunarOut> = None;
    let mut final_aim: Option<Vector3> = None;

    match (landing, conic_landing) {
        (Some(target), Some(conic)) => {
            let window = &target.window;
            let period = std::f64::consts::TAU * (r_llo.powi(3) / mu_m).sqrt();
            let max_wait = window.max_wait_in_orbit_h * 3600.0;
            let min_wait = target.descent.config.min_time_in_orbit_h * 3600.0;
            // Цель плоскости подлёта: направление на точку посадки в расчётный момент, в плоскости из поиска.
            let s_conic = site_dir_at(conic.epoch)?;
            let aim0 = (s_conic - h_lunar * h_lunar.dot(&s_conic)).normalize();
            let mut extra = 0u32;
            let mut kept: Option<(RefinedLanding, Vector3, Peri, bool, LunarOut, Vector3)> = None;
            loop {
                let mut last: Option<(LunarOut, Epoch, f64, Vector3, Peri, Vector3)> = None;
                let mut reached = false;
                // Поворот плоскости подлёта вокруг направления подлёта подбирается методом секущих.
                let (mut phi, mut prev): (f64, Option<(f64, f64)>) = (0.0, None);
                for _ in 0..cfg.max_iterations {
                    let aim = rotate_about(aim0, v_hat, phi);
                    let (u_new, p, it, ok, n) =
                        match solve_auto(&mut shooter, u, Some(&aim), cfg.max_iterations)? {
                            Ok(x) => x,
                            Err(reason) => return Ok(Err(reason)),
                        };
                    u = u_new;
                    iterations += it + 1;
                    steps_total += n;
                    if !ok {
                        return Ok(Err("trajectory targeting did not converge".into()));
                    }
                    // Выход на круговую орбиту в периселении.
                    let t_peri = t0 + p.t * Unit::Second;
                    let v_circ = (mu_m / p.r.norm()).sqrt();
                    let dv_loi = p.v.norm() - v_circ;
                    let y_llo = state(&p.r, &(p.v.normalize() * v_circ));
                    let moon = MoonCentered::new(almanac, t_peri, &field, mu_e, forces);
                    let phase = LunarPhase {
                        almanac,
                        dynamics: &moon,
                        target,
                        tol: cfg.tolerance,
                        mu: mu_m,
                        spent_m_s: spent_before_descent((u - v_park).norm(), dv_loi),
                    };
                    // Тот же виток, что и по формулам (плюс дополнительные, если условия не выполнились).
                    let t_doi_min = min_wait
                        .max((conic.doi_epoch - t_peri).to_seconds() - 0.5 * period)
                        + extra as f64 * period;
                    let out = match phase.run(&y_llo, t_doi_min)? {
                        Ok(o) => o,
                        Err(reason) => return Ok(Err(reason.into())),
                    };
                    steps_total += out.steps;
                    let cross = out.cross;
                    let miss = out.error.norm() * MOON_RADIUS_KM;

                    last = Some((out, t_peri, miss, u, p, aim));
                    if miss < cfg.landing_tolerance_km {
                        reached = true;
                        break;
                    }
                    let next = match prev {
                        Some((phi0, cross0)) if (cross - cross0).abs() > 1e-12 => {
                            phi - cross * (phi - phi0) / (cross - cross0)
                        }
                        // Первый пробный поворот — на величину бокового промаха.
                        _ => phi + if cross.abs() > 1e-9 { cross } else { 1e-4 },
                    };
                    prev = Some((phi, cross));
                    phi = next.clamp(phi - 0.2, phi + 0.2);
                }
                let (out, t_peri, miss, u_at, p_at, aim_at) = last.expect("at least one iteration");
                let touchdown = t_peri
                    + (out.t_stop + target.descent.config.terminal_duration_s) * Unit::Second;
                let vis = target
                    .site
                    .evaluate(&sky_at(almanac, &target.site.site, touchdown)?);
                let sun_ok = !window.require_sun_visible
                    || (vis.sun.visible
                        && vis.sun.elevation_deg >= window.sun_elevation_min_deg
                        && vis.sun.elevation_deg <= window.sun_elevation_max_deg);
                let earth_ok = !window.require_earth_visible || vis.earth.visible;
                // Лимит ожидания на орбите (минута допуска на численный шум).
                let wait_ok = out.t_doi <= max_wait + 60.0;
                let pdi_r = pos(&out.pdi).norm();
                let info = RefinedLanding {
                    orbit_inclination_deg: out.inclination_deg,
                    orbit_min_altitude_km: out.min_alt,
                    orbit_max_altitude_km: out.max_alt,
                    doi_epoch: t_peri + out.t_doi * Unit::Second,
                    doi_delta_v_m_s: out.doi_dv_m_s,
                    wait_in_orbit_s: out.t_doi,
                    extra_orbits: extra,
                    pdi_epoch: t_peri + out.t_pdi * Unit::Second,
                    pdi_altitude_km: pdi_r - MOON_RADIUS_KM,
                    pdi_height_above_site_km: pdi_r - target.site.site.position_km().norm(),
                    pdi_speed_km_s: vel(&out.pdi).norm(),
                    braking: DescentSim {
                        samples: Vec::new(),
                        ..out.sim.clone()
                    },
                    terminal_duration_s: target.descent.config.terminal_duration_s,
                    touchdown_epoch: touchdown,
                    touchdown_shift_s: (touchdown - conic.epoch).to_seconds(),
                    descent_delta_v_m_s: out.doi_dv_m_s
                        + out.sim.delta_v_m_s
                        + target.descent.config.terminal_allowance_m_s,
                    miss_km: miss,
                    sun: vis.sun,
                    earth: vis.earth,
                    conditions_met: sun_ok && earth_ok && wait_ok,
                    unmet_reason: match (sun_ok, earth_ok) {
                        (true, true) if !wait_ok => Some(REASON_LANDING_WAIT),
                        (true, true) => None,
                        (false, true) => Some(REASON_LANDING_SUN),
                        (true, false) => Some(REASON_LANDING_EARTH),
                        (false, false) => Some(REASON_LANDING_BOTH),
                    },
                };
                let met = info.conditions_met;
                let next_wait = info.wait_in_orbit_s + period;
                if met || kept.is_none() {
                    kept = Some((info, u_at, p_at, reached, out, aim_at));
                }
                // Условия не выполнились на этом витке: пробуем следующие, пока позволяет ожидание.
                if met || !reached || next_wait > max_wait || extra >= 12 {
                    break;
                }
                extra += 1;
            }
            let (info, u_kept, p_kept, reached, out_kept, aim_kept) =
                kept.expect("at least one attempt");
            final_aim = Some(aim_kept);
            u = u_kept;
            converged = reached;
            peri = p_kept;
            refined_landing = Some(info);
            lunar_out = Some(out_kept);
        }
        _ => {
            let (u_new, p, it, ok, n) = match solve_auto(&mut shooter, u, None, cfg.max_iterations)?
            {
                Ok(x) => x,
                Err(reason) => return Ok(Err(reason)),
            };
            u = u_new;
            iterations += it;
            steps_total += n;
            converged = ok;
            peri = p;
        }
    }

    let p = peri;
    let arrival_epoch = t0 + p.t * Unit::Second;
    if ((u - v_park).norm() - dep.delta_v_magnitude_km_s).abs() > 0.3 {
        return Ok(Err("trajectory targeting diverged".into()));
    }

    // Разгон конечной длительности: тот же прилёт (момент, высота и плоскость периселения),
    // что у мгновенного импульса, — поэтому всё, что дальше (окололунная орбита, посадка),
    // остаётся в силе.
    let u_impulsive = u;
    let mut burn_record: Option<(BurnModel, (f64, State), Vector3)> = None;
    let mut tli_burn: Option<TliBurn> = None;
    if let Some(eng) = tli_engine {
        let dv_imp = (u - v_park).norm();
        let base_model = BurnModel {
            thrust_n: eng.thrust_n,
            c_m_s: eng.isp_s * G0,
            m0: eng.mass_kg,
            normal: h_park.normalize(),
            omega: 0.0,
        };
        // Точка опоры — раньше самого раннего возможного начала разгона.
        let t_a = (t_tli - 0.6 * base_model.duration(dv_imp * 1.3)).max(t_ins);
        let anchor = (
            t_a,
            integrate(&earth, t_ins, &y_ins, t_a, cfg.tolerance, None, None)?.y,
        );
        let n_orbit = v_park.norm() / r1.norm();
        shooter.fix_time = true;
        shooter.decided = true;
        shooter.coarse_plane = None;
        shooter.t_target = p.t;
        let aim = final_aim;
        // Скорость поворота тяги — доля угловой скорости движения по орбите; подбирается
        // на наименьший Δv (сетка, затем уточнение золотым сечением).
        // Неудачная скорость поворота может увести пробную траекторию в Землю — это просто
        // «не годится», а не ошибка всего расчёта. Подбор продолжается от лучшей найденной скорости.
        let mut u_guess = u_impulsive;
        let h_approach = p.r.cross(&p.v).normalize();
        let mut coarse_steps = 0u64;
        let mut solve_at =
            |k: f64, guess: Option<Vector3>| -> Result<Option<(f64, Vector3, Peri, u64)>> {
                let model = BurnModel {
                    omega: k * n_orbit,
                    ..base_model
                };
                // Без готового приближения — орбита после разгона подгоняется под импульсную.
                let guess = match guess {
                    Some(g) => g,
                    None => {
                        v_park
                            + match_post_burn_orbit(
                                &earth,
                                anchor,
                                &model,
                                t_tli,
                                &r1,
                                &v_park,
                                &(u_impulsive - v_park),
                                cfg.tolerance,
                            )?
                    }
                };
                shooter.burn = Some((model, anchor));
                // Как и для мгновенного импульса: грубое наведение по картинной плоскости (на ту же
                // плоскость подлёта), затем точное — с прежним моментом прилёта.
                shooter.coarse_plane = Some(h_approach);
                let guess = match shooter.solve(guess, None, cfg.max_iterations) {
                    Ok(Ok((uc, _, _, _, n))) => {
                        coarse_steps += n;
                        uc
                    }
                    _ => guess,
                };
                shooter.coarse_plane = None;
                Ok(
                    match shooter.solve(guess, aim.as_ref(), cfg.max_iterations) {
                        Ok(Ok((uf, pf, _, true, n))) => Some(((uf - v_park).norm(), uf, pf, n)),
                        _ => None,
                    },
                )
            };
        // Δv почти параболически зависит от скорости поворота: три точки, вершина параболы.
        let mut best: Option<(f64, f64, Vector3, Peri)> = None;
        let mut pts: Vec<(f64, f64)> = Vec::new();
        for k in [0.4, 0.7, 1.0] {
            if let Some((dv, uf, pf, n)) = solve_at(k, None)? {
                steps_total += n;
                pts.push((k, dv));
                if best.as_ref().is_none_or(|b| dv < b.1) {
                    best = Some((k, dv, uf, pf));
                    u_guess = uf;
                }
            }
        }
        if pts.len() == 3 {
            let ((x0, y0), (x1, y1), (x2, y2)) = (pts[0], pts[1], pts[2]);
            let den = (x0 - x1) * (x0 - x2) * (x1 - x2);
            let a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / den;
            let b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / den;
            if a > 0.0 {
                let kv = (-b / (2.0 * a)).clamp(0.0, 1.5);
                if let Some((dv, uf, pf, n)) = solve_at(kv, Some(u_guess))? {
                    steps_total += n;
                    if best.as_ref().is_none_or(|bb| dv < bb.1) {
                        best = Some((kv, dv, uf, pf));
                    }
                }
            }
        }
        if best.is_none() {
            // Запасной вариант: направление тяги не поворачивается.
            if let Some((dv, uf, pf, n)) = solve_at(0.0, None)? {
                steps_total += n;
                best = Some((0.0, dv, uf, pf));
            }
        }
        steps_total += coarse_steps;
        let Some(best) = best else {
            return Ok(Err(
                "targeting with a finite trans-lunar injection burn did not converge".into(),
            ));
        };
        let (k, dv_fin, u_fin, _p_fin) = best;
        let model = BurnModel {
            omega: k * n_orbit,
            ..base_model
        };
        let (t_s, t_e, y_e, _) = burn_through(
            &earth,
            anchor,
            &model,
            t_tli,
            &(u_fin - v_park),
            cfg.tolerance,
        )?;
        let max_dv =
            eng.isp_s * G0 * (eng.mass_kg / (eng.mass_kg - eng.propellant_kg).max(1e-9)).ln()
                / 1000.0;
        let used = model.propellant(dv_fin);
        tli_burn = Some(TliBurn {
            engine: eng.clone(),
            start_epoch: t0 + t_s * Unit::Second,
            end_epoch: t0 + t_e * Unit::Second,
            duration_s: t_e - t_s,
            delta_v_km_s: dv_fin,
            impulsive_delta_v_km_s: dv_imp,
            gravity_loss_m_s: (dv_fin - dv_imp) * 1000.0,
            propellant_used_kg: used,
            max_delta_v_km_s: max_dv,
            exceeds_propellant: used > eng.propellant_kg,
            impulsive_exceeds_propellant: model.propellant(dv_imp) > eng.propellant_kg,
            steering_rate_fraction: k,
            burn_arc_deg: (n_orbit * (t_e - t_s)).to_degrees(),
            end_position_km: pos(&y_e),
            end_velocity_km_s: vel(&y_e),
        });
        burn_record = Some((model, anchor, u_fin - v_park));
        u = u_fin;
    }
    let dv_tli = (u - v_park).norm();
    let dv_loi = p.v.norm() - (mu_m / p.r.norm()).sqrt();
    let tof_s = p.t - t_tli;
    let path = record_path(
        almanac,
        &earth,
        &field,
        mu_e,
        forces,
        cfg.tolerance,
        &Launch {
            y_ins,
            t_ins,
            ascent: &sol.ascent,
            plane_normal: h_park,
            parking_radius_km: sol.parking_orbit.radius_km,
            launch_r: sol.launch.site_direction_j2000 * sol.ascent.pad_radius_km,
            t_tli,
            v_park,
            u: u_impulsive,
            burn: burn_record,
        },
        &p,
        lunar_out.as_ref().zip(landing),
    )?;
    let refined_budget = budget.map(|b| {
        let tof_limit = sol
            .budget
            .as_ref()
            .is_some_and(|c| c.time_of_flight_at_lower_limit);
        b.check(
            dv_tli,
            dv_loi,
            refined_landing.as_ref().map(|l| l.descent_delta_v_m_s),
            tof_limit,
        )
    });
    Ok(Ok(Refined {
        config: cfg.clone(),
        converged,
        iterations,
        integration_steps: steps_total,
        parking_position_shift_km: (r1 - dep.position_km).norm(),
        tli_position_km: r1,
        tli_velocity_after_km_s: u,
        tli_delta_v_km_s: dv_tli,
        tli_delta_v_change_m_s: (dv_tli - dep.delta_v_magnitude_km_s) * 1000.0,
        arrival_epoch,
        arrival_shift_s: (arrival_epoch - arr.epoch).to_seconds(),
        arrival_time_kept: shooter.fix_time,
        time_of_flight_s: tof_s,
        periselene_altitude_km: p.r.norm() - MOON_RADIUS_KM,
        periselene_speed_km_s: p.v.norm(),
        loi_delta_v_km_s: dv_loi,
        loi_delta_v_change_m_s: (dv_loi - arr.delta_v_loi_km_s) * 1000.0,
        landing: refined_landing,
        budget: refined_budget,
        tli_burn,
        moon_gm_km3_s2: mu_m,
        path: Some(path),
    }))
}

/// Старт и разгон — для записи траектории.
struct Launch<'a> {
    /// Состояние при выходе на опорную орбиту.
    y_ins: State,
    /// Момент выхода на опорную орбиту, с от старта.
    t_ins: f64,
    ascent: &'a crate::trajectory::AscentInfo,
    plane_normal: Vector3,
    parking_radius_km: f64,
    /// Точка старта на поверхности, J2000.
    launch_r: Vector3,
    t_tli: f64,
    v_park: Vector3,
    /// Скорость после разгона (при мгновенном импульсе).
    u: Vector3,
    /// Разгон конечной длительности: модель, точка опоры, Δv двигателя (вектор в середине разгона).
    burn: Option<(BurnModel, (f64, State), Vector3)>,
}

/// Интегрирует участок, записывая состояние через промежутки `dt(t, y)`. Возвращает конечное состояние.
fn sampled(
    d: &dyn crate::propagate::Dynamics,
    t_start: f64,
    y0: &State,
    t_end: f64,
    tol: f64,
    dt: impl Fn(f64, &State) -> Result<f64>,
    mut record: impl FnMut(f64, &State) -> Result<()>,
) -> Result<State> {
    let (mut t, mut y) = (t_start, *y0);
    record(t, &y)?;
    while t < t_end - 1e-6 {
        let step = dt(t, &y)?.min(t_end - t);
        let out = integrate(d, t, &y, t + step, tol, None, None)?;
        (t, y) = (out.t, out.y);
        record(t, &y)?;
    }
    Ok(y)
}

/// Повторяет итоговый полёт и записывает точки и события для журнала траектории.
#[allow(clippy::too_many_arguments)]
fn record_path(
    almanac: &Almanac,
    earth: &EarthCentered,
    field: &GravityField,
    mu_e: f64,
    forces: Forces,
    tol: f64,
    launch: &Launch,
    p: &Peri,
    lunar: Option<(&LunarOut, &LandingTarget)>,
) -> Result<TrajectoryPath> {
    use crate::trajlog::{Bodies, Phase, km};
    const EARTH_ROTATION_RAD_S: f64 = 7.292_115e-5;
    let bodies = Bodies {
        almanac,
        t0: earth.t0,
    };
    let mut path = TrajectoryPath::default();
    let (ev, sm) = (&mut path.events, &mut path.samples);

    // Старт и подъём до опорной орбиты (кинематическая модель, см. AscentInfo).
    let v_ground = (earth.pole * EARTH_ROTATION_RAD_S).cross(&launch.launch_r);
    ev.push(bodies.event(
        "launch",
        "launch from the pad",
        0.0,
        &launch.launch_r,
        &v_ground,
        None,
        false,
    )?);
    match (&launch.ascent.simulation, launch.ascent.to_j2000) {
        (Some(sim), Some(m)) => {
            // Моделирование ракеты: точки каждую секунду и события подъёма.
            for p in &sim.points {
                sm.push(bodies.sample(Phase::Ascent, p.t_s, &(m * p.r), &(m * p.v), false)?);
            }
            for e in sim
                .events
                .iter()
                .filter(|e| !matches!(e.name, "liftoff" | "insertion"))
            {
                ev.push(bodies.event(
                    e.name,
                    e.detail.clone(),
                    e.t_s,
                    &(m * e.r),
                    &(m * e.v),
                    None,
                    false,
                )?);
            }
        }
        _ => {
            let n_asc = (launch.t_ins / 10.0).ceil().max(1.0) as usize;
            for i in 0..=n_asc {
                let t = launch.t_ins * i as f64 / n_asc as f64;
                let (r, v) = if i == 0 {
                    (launch.launch_r, v_ground)
                } else {
                    launch
                        .ascent
                        .state_at(t, &launch.plane_normal, launch.parking_radius_km)
                };
                sm.push(bodies.sample(Phase::Ascent, t, &r, &v, false)?);
            }
        }
    }
    let y = launch.y_ins;
    ev.push(bodies.event(
        "parking_insertion",
        "parking orbit insertion",
        launch.t_ins,
        &pos(&y),
        &vel(&y),
        None,
        false,
    )?);
    // Разгон к Луне: мгновенный импульс или работа двигателя.
    let (t_transfer, y_transfer) = match &launch.burn {
        None => {
            sampled(
                earth,
                launch.t_ins,
                &y,
                launch.t_tli,
                tol,
                |_, _| Ok(30.0),
                |t, y| {
                    sm.push(bodies.sample(Phase::ParkingOrbit, t, &pos(y), &vel(y), false)?);
                    Ok(())
                },
            )?;
            let r1 = integrate(
                earth,
                launch.t_ins,
                &launch.y_ins,
                launch.t_tli,
                tol,
                None,
                None,
            )
            .map(|o| pos(&o.y))?;
            ev.push(bodies.event(
                "tli",
                "trans-lunar injection burn",
                launch.t_tli,
                &r1,
                &launch.u,
                Some(launch.u - launch.v_park),
                false,
            )?);
            (launch.t_tli, state(&r1, &launch.u))
        }
        Some((m, anchor, dv)) => {
            // Опорная орбита до точки опоры, от неё до начала разгона — так же, как при подборе.
            let y_a = sampled(
                earth,
                launch.t_ins,
                &y,
                anchor.0,
                tol,
                |_, _| Ok(30.0),
                |t, y| {
                    sm.push(bodies.sample(Phase::ParkingOrbit, t, &pos(y), &vel(y), false)?);
                    Ok(())
                },
            )?;
            let _ = y_a;
            let (t_s, t_e) = m.window(dv.norm(), launch.t_tli, anchor.0);
            let y_s = if t_s > anchor.0 {
                let o = integrate(earth, anchor.0, &anchor.1, t_s, tol, None, None)?;
                sm.push(bodies.sample(Phase::ParkingOrbit, t_s, &pos(&o.y), &vel(&o.y), false)?);
                o.y
            } else {
                anchor.1
            };
            ev.push(bodies.event(
                "tli",
                format!(
                    "trans-lunar injection burn starts (engine burn {:.0} s, delta-v {:.0} m/s)",
                    t_e - t_s,
                    dv.norm() * 1000.0
                ),
                t_s,
                &pos(&y_s),
                &vel(&y_s),
                None,
                false,
            )?);
            let th = Thrusting {
                base: earth,
                m: *m,
                t_s,
                t_mid: 0.5 * (t_s + t_e),
                dir: dv.normalize(),
            };
            let y_e = sampled(
                &th,
                t_s,
                &y_s,
                t_e,
                tol,
                |_, _| Ok(10.0),
                |t, y| {
                    sm.push(bodies.sample(Phase::TliBurn, t, &pos(y), &vel(y), false)?);
                    Ok(())
                },
            )?;
            ev.push(bodies.event(
                "tli_end",
                "trans-lunar injection burn ends",
                t_e,
                &pos(&y_e),
                &vel(&y_e),
                None,
                false,
            )?);
            (t_e, y_e)
        }
    };

    // Перелёт до периселения.
    let near =
        |t: f64, y: &State| -> Result<bool> { Ok((pos(y) - earth.moon(t)?.0).norm() < 70_000.0) };
    sampled(
        earth,
        t_transfer,
        &y_transfer,
        p.t,
        tol,
        |t, y| Ok(if near(t, y)? { 30.0 } else { 300.0 }),
        |t, y| {
            let mut s = bodies.sample(Phase::Transfer, t, &pos(y), &vel(y), false)?;
            if near(t, y)? {
                s.r_moon = Some(km(&(pos(y) - earth.moon(t)?.0)));
            }
            sm.push(s);
            Ok(())
        },
    )?;

    // Выход на окололунную орбиту.
    let mu_m = field.gm;
    let v_circ = p.v.normalize() * (mu_m / p.r.norm()).sqrt();
    ev.push(bodies.event(
        "loi",
        "lunar orbit insertion burn (periselene)",
        p.t,
        &p.r,
        &v_circ,
        Some(v_circ - p.v),
        true,
    )?);
    let t_peri = earth.t0 + p.t * Unit::Second;
    let moon = MoonCentered::new(almanac, t_peri, field, mu_e, forces);
    let y_llo = state(&p.r, &v_circ);
    let lunar_sample = |sm: &mut Vec<_>, phase: Phase, t: f64, y: &State| -> Result<()> {
        sm.push(bodies.sample(phase, p.t + t, &pos(y), &vel(y), true)?);
        Ok(())
    };

    let Some((out, target)) = lunar else {
        // Без точки посадки: один виток окололунной орбиты.
        let period = std::f64::consts::TAU * (p.r.norm().powi(3) / mu_m).sqrt();
        sampled(
            &moon,
            0.0,
            &y_llo,
            period,
            tol,
            |_, _| Ok(60.0),
            |t, y| lunar_sample(sm, Phase::LunarOrbit, t, y),
        )?;
        return Ok(path);
    };

    // Ожидание на орбите до схода.
    sampled(
        &moon,
        0.0,
        &y_llo,
        out.t_doi,
        tol,
        |_, _| Ok(60.0),
        |t, y| lunar_sample(sm, Phase::LunarOrbit, t, y),
    )?;
    let (pre, post) = (out.doi_pre, out.doi_post);
    ev.push(bodies.event(
        "doi",
        "descent orbit insertion burn",
        p.t + out.t_doi,
        &pos(&post),
        &vel(&post),
        Some(vel(&post) - vel(&pre)),
        true,
    )?);
    // Эллипс спуска до начала торможения.
    sampled(
        &moon,
        out.t_doi,
        &post,
        out.t_pdi,
        tol,
        |_, _| Ok(10.0),
        |t, y| lunar_sample(sm, Phase::DescentCoast, t, y),
    )?;
    let pdi = out.pdi;
    ev.push(bodies.event(
        "pdi",
        "powered descent initiation",
        p.t + out.t_pdi,
        &pos(&pdi),
        &vel(&pdi),
        None,
        true,
    )?);
    // Торможение.
    for (ts, r, v) in &out.sim.samples {
        lunar_sample(sm, Phase::Braking, out.t_pdi + ts, &state(r, v))?;
    }
    let (t_end, r_end, v_end) = *out.sim.samples.last().expect("braking samples");
    ev.push(bodies.event(
        "braking_end",
        "braking stops above the landing site",
        p.t + out.t_pdi + t_end,
        &r_end,
        &v_end,
        None,
        true,
    )?);
    // Заход и вертикальный участок: точка посадки вращается вместе с Луной.
    let site_radius = target.site.site.position_km().norm();
    let site_up = target.site.site.up();
    let t_stop = out.t_pdi + t_end;
    let terminal = target.descent.config.terminal_duration_s;
    let omega_m = |t: f64| -> Result<Vector3> {
        Ok(bodies.moon_dcm(p.t + t)?
            * Vector3::new(0.0, 0.0, 1.0)
            * crate::lunar_landing::MOON_ROTATION_RATE_RAD_S)
    };
    // Заход не моделируется двигателем: точки плавно ведут от остановки после торможения
    // (несколько метров над точкой и в стороне, малая скорость) к месту посадки. Кривая —
    // кубический сплайн Эрмита в осях, вращающихся с Луной: в начале положение и скорость
    // как в конце торможения, в конце — точка посадки и нулевая скорость относительно грунта.
    let site_body = site_up * site_radius;
    let dcm_stop = bodies.moon_dcm(p.t + t_stop)?;
    let d0 = dcm_stop.transpose() * r_end - site_body;
    let u0 = dcm_stop.transpose() * (v_end - omega_m(t_stop)?.cross(&r_end));
    // Длительность сглаживания: не дольше, чем нужно, чтобы кривая не проскочила точку посадки
    // (у сплайна Эрмита это |скорость| · τ ≤ 3 · |смещение|, по вертикали — отдельно, чтобы не
    // уйти под грунт). Остаток запаса времени аппарат стоит на месте посадки.
    let up_body = site_up.normalize();
    let (d_up, u_up) = (d0.dot(&up_body), u0.dot(&up_body));
    let mut tau = terminal;
    if u0.norm() > 1e-9 {
        tau = tau.min(3.0 * d0.norm() / u0.norm());
    }
    if u_up < 0.0 && d_up > 0.0 {
        tau = tau.min(3.0 * d_up / -u_up);
    }
    let tau = tau.max(1.0).min(terminal.max(1.0));
    let site_state = |t: f64| -> Result<State> {
        let s = ((t - t_stop) / tau).clamp(0.0, 1.0);
        let (h00, h10) = (
            2.0 * s.powi(3) - 3.0 * s * s + 1.0,
            s.powi(3) - 2.0 * s * s + s,
        );
        let (dh00, dh10) = (
            (6.0 * s * s - 6.0 * s) / tau,
            (3.0 * s * s - 4.0 * s + 1.0) / tau,
        );
        let p_body = site_body + d0 * h00 + u0 * (h10 * tau);
        let v_body = d0 * dh00 + u0 * (dh10 * tau);
        let dcm = bodies.moon_dcm(p.t + t)?;
        let r = dcm * p_body;
        Ok(state(&r, &(dcm * v_body + omega_m(t)?.cross(&r))))
    };
    // Точки: каждые 0,5 с на сглаживании, дальше каждые 5 с.
    let mut t = t_stop;
    while t < t_stop + terminal - 1e-9 {
        lunar_sample(sm, Phase::Terminal, t, &site_state(t)?)?;
        t += if t < t_stop + tau { 0.5 } else { 5.0 };
        t = t.min(t_stop + terminal);
    }
    lunar_sample(
        sm,
        Phase::Terminal,
        t_stop + terminal,
        &site_state(t_stop + terminal)?,
    )?;
    let y_td = site_state(t_stop + terminal)?;
    ev.push(bodies.event(
        "touchdown",
        "touchdown at the landing site",
        p.t + t_stop + terminal,
        &pos(&y_td),
        &vel(&y_td),
        None,
        true,
    )?);
    Ok(path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn config_validation() {
        assert!(NumericalConfig::default().validate().is_ok());
        assert!(
            NumericalConfig {
                lunar_gravity_degree: 61,
                ..Default::default()
            }
            .validate()
            .is_err()
        );
        assert!(
            NumericalConfig {
                tolerance: 1.0,
                ..Default::default()
            }
            .validate()
            .is_err()
        );
        assert!(
            NumericalConfig {
                max_iterations: 0,
                ..Default::default()
            }
            .validate()
            .is_err()
        );
    }

    #[test]
    fn small_linear_solver() {
        let x = solve_small(
            vec![
                vec![2.0, 1.0, 0.0],
                vec![1.0, 3.0, 1.0],
                vec![0.0, 1.0, 4.0],
            ],
            vec![3.0, 5.0, 5.0],
        )
        .unwrap();
        for (a, b) in x.iter().zip([1.0, 1.0, 1.0]) {
            assert!((a - b).abs() < 1e-12);
        }
        assert!(
            solve_small(vec![vec![1.0, 2.0], vec![2.0, 4.0]], vec![1.0, 2.0])
                .is_none_or(|x| x.iter().any(|v| !v.is_finite()))
        );
    }

    use crate::ephemeris::{AlmanacSource, load_almanac_blocking};
    use crate::lunar_landing::{DescentModel, LandingWindowConfig};
    use crate::spaceport::SpaceportCatalog;
    use crate::trajectory::{DateConstraint, DirectTransfer, Trajectory};
    use crate::visibility::{SiteVisibility, VisibilityConfig};

    fn real() -> (Almanac, SpaceportCatalog) {
        let almanac =
            load_almanac_blocking(AlmanacSource::LocalDir("data/kernels".into())).unwrap();
        let text = std::fs::read_to_string("data/spaceports.toml").unwrap();
        (
            almanac,
            SpaceportCatalog::parse(&text, std::path::Path::new("data/spaceports.toml")).unwrap(),
        )
    }

    fn target(lat: f64, lon: f64) -> LandingTarget {
        LandingTarget {
            site: SiteVisibility::new(None, lat, lon, 0.0, VisibilityConfig::default()).unwrap(),
            window: LandingWindowConfig::default(),
            descent: DescentModel::default(),
        }
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_refinement_reaches_lunar_orbit() {
        let (almanac, catalog) = real();
        let ksc = catalog.get("ksc").unwrap();
        let day = Epoch::from_gregorian_utc_at_midnight(2026, 12, 10);
        let rep = DirectTransfer::default()
            .solve(&almanac, ksc, DateConstraint::LaunchDay(day), None, None)
            .unwrap();
        let sol = rep.solution.unwrap();
        let cfg = NumericalConfig::default();
        let r = refine(&almanac, &cfg, &sol, None, None, None)
            .unwrap()
            .unwrap();
        assert!(r.converged, "{r:?}");
        assert!(
            (r.periselene_altitude_km - sol.arrival.lunar_orbit_altitude_km).abs() < 0.05,
            "{}",
            r.periselene_altitude_km
        );
        // Поправки к расчёту по формулам небольшие: десятки м/с.
        assert!(
            r.tli_delta_v_change_m_s.abs() < 40.0,
            "{}",
            r.tli_delta_v_change_m_s
        );
        assert!(
            r.loi_delta_v_change_m_s.abs() < 120.0,
            "{}",
            r.loi_delta_v_change_m_s
        );
        assert!(r.landing.is_none() && r.budget.is_none());
        // Без точки посадки в журнале траектории — до одного витка окололунной орбиты.
        let path = r.path.as_ref().unwrap();
        assert_eq!(
            path.events.iter().map(|e| e.name).collect::<Vec<_>>(),
            ["launch", "parking_insertion", "tli", "loi"]
        );
        assert!(
            path.samples
                .iter()
                .any(|s| s.phase == crate::trajlog::Phase::LunarOrbit)
        );

        // Независимая проверка: полёт с найденной скоростью после разгона приходит на ту же высоту.
        let forces = Forces {
            earth_j2: true,
            sun: true,
        };
        let earth = EarthCentered::new(
            &almanac,
            sol.departure.epoch,
            sol.mu_earth_km3_s2,
            sol.mu_moon_km3_s2,
            forces,
        )
        .unwrap();
        let t_end = (r.arrival_epoch - sol.departure.epoch).to_seconds();
        let out = integrate(
            &earth,
            0.0,
            &state(&r.tli_position_km, &r.tli_velocity_after_km_s),
            t_end,
            1e-11,
            None,
            None,
        )
        .unwrap();
        let moon = earth.moon(t_end).unwrap().0;
        let alt = (pos(&out.y) - moon).norm() - MOON_RADIUS_KM;
        assert!((alt - r.periselene_altitude_km).abs() < 0.5, "{alt}");

        // Без сжатия Земли и Солнца ответ другой: силы действительно учитываются.
        let plain = NumericalConfig {
            earth_j2: false,
            sun_gravity: false,
            ..Default::default()
        };
        let r2 = refine(&almanac, &plain, &sol, None, None, None)
            .unwrap()
            .unwrap();
        assert!(r2.converged);
        assert!((r2.tli_velocity_after_km_s - r.tli_velocity_after_km_s).norm() > 1e-5);
        assert!(r2.parking_position_shift_km < r.parking_position_shift_km);
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_finite_tli_burn_keeps_the_arrival_and_costs_gravity_loss() {
        let (almanac, catalog) = real();
        let ksc = catalog.get("ksc").unwrap();
        let day = Epoch::from_gregorian_utc_at_midnight(2026, 12, 10);
        let sol = DirectTransfer::default()
            .solve(&almanac, ksc, DateConstraint::LaunchDay(day), None, None)
            .unwrap()
            .solution
            .unwrap();
        let cfg = NumericalConfig::default();
        let imp = refine(&almanac, &cfg, &sol, None, None, None)
            .unwrap()
            .unwrap();
        assert!(imp.tli_burn.is_none());
        // Разгонная ступень уровня ICPS с грузом ~11,6 т.
        let eng = Engine {
            name: "test stage".into(),
            thrust_n: 110_100.0,
            isp_s: 462.0,
            mass_kg: 44_360.0,
            propellant_kg: 28_948.0,
        };
        let r = refine(&almanac, &cfg, &sol, None, None, Some(&eng))
            .unwrap()
            .unwrap();
        let b = r.tli_burn.as_ref().expect("finite burn");
        // Прилёт тот же, что при мгновенном импульсе.
        assert!((r.arrival_epoch - imp.arrival_epoch).to_seconds().abs() < 1.0);
        assert!((r.periselene_altitude_km - imp.periselene_altitude_km).abs() < 0.05);
        assert!((b.impulsive_delta_v_km_s - imp.tli_delta_v_km_s).abs() < 1e-9);
        // Потери положительные и реалистичные для разгона ~15 мин; Δv двигателя = импульс + потери.
        assert!(
            b.gravity_loss_m_s > 20.0 && b.gravity_loss_m_s < 200.0,
            "{}",
            b.gravity_loss_m_s
        );
        assert!((r.tli_delta_v_km_s - b.delta_v_km_s).abs() < 1e-12);
        // Длительность и топливо — по формуле Циолковского; разгон охватывает момент импульса.
        let c = eng.isp_s * G0;
        let used = eng.mass_kg * (1.0 - (-b.delta_v_km_s * 1000.0 / c).exp());
        assert!((b.propellant_used_kg - used).abs() < 1e-6);
        assert!((b.duration_s - used * c / eng.thrust_n).abs() < 1e-6);
        assert!(b.start_epoch < sol.departure.epoch && b.end_epoch > sol.departure.epoch);
        assert!(!b.exceeds_propellant);
        // Траектория: этап tli_burn между событиями tli и tli_end, без скачков скорости на стыках.
        let path = r.path.as_ref().unwrap();
        let names: Vec<&str> = path.events.iter().map(|e| e.name).collect();
        assert_eq!(
            names,
            ["launch", "parking_insertion", "tli", "tli_end", "loi"]
        );
        use crate::trajlog::Phase;
        let sm = &path.samples;
        assert!(sm.iter().filter(|s| s.phase == Phase::TliBurn).count() > 50);
        for w in sm
            .windows(2)
            .filter(|w| w[0].t_s == w[1].t_s && w[0].phase != w[1].phase)
        {
            if matches!(
                (w[0].phase, w[1].phase),
                (Phase::ParkingOrbit, Phase::TliBurn) | (Phase::TliBurn, Phase::Transfer)
            ) {
                let dv: f64 = (0..3)
                    .map(|i| (w[0].v[i] - w[1].v[i]).powi(2))
                    .sum::<f64>()
                    .sqrt();
                assert!(dv < 1e-5, "jump {dv} at {}", w[0].t_s);
            }
        }
        // Независимая проверка: от конца разгона полёт приходит на ту же высоту.
        let forces = Forces {
            earth_j2: true,
            sun: true,
        };
        let earth = EarthCentered::new(
            &almanac,
            b.end_epoch,
            sol.mu_earth_km3_s2,
            sol.mu_moon_km3_s2,
            forces,
        )
        .unwrap();
        let t_end = (r.arrival_epoch - b.end_epoch).to_seconds();
        let out = integrate(
            &earth,
            0.0,
            &state(&b.end_position_km, &b.end_velocity_km_s),
            t_end,
            1e-11,
            None,
            None,
        )
        .unwrap();
        let alt = (pos(&out.y) - earth.moon(t_end).unwrap().0).norm() - MOON_RADIUS_KM;
        assert!((alt - r.periselene_altitude_km).abs() < 0.5, "{alt}");
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_simulated_ascent_flows_into_the_trajectory() {
        use crate::launch_vehicle::{LaunchVehicleLibrary, LaunchVehicleParams};
        use crate::trajectory::AscentSetup;
        let (almanac, catalog) = real();
        let ksc = catalog.get("ksc").unwrap();
        let lib = LaunchVehicleLibrary::presets_only().unwrap();
        let rk = lib
            .resolve("sls_block1", &LaunchVehicleParams::default())
            .unwrap();
        let setup = AscentSetup {
            rocket_id: "sls_block1".into(),
            rocket_name: rk.config.name.clone(),
            vehicle: rk.config.ascent.clone().unwrap(),
            upper: rk.config.upper_stage(),
            payload_kg: 16103.0,
        };
        let dt = DirectTransfer::default()
            .with_ascent([("ksc".to_string(), setup)].into_iter().collect());
        let day = Epoch::from_gregorian_utc_at_midnight(2027, 1, 1);
        let t = target(5.0, -10.0);
        let rep = dt
            .solve(
                &almanac,
                ksc,
                DateConstraint::ArrivalDay(day),
                None,
                Some(&t),
            )
            .unwrap();
        let sol = rep.solution.expect("solution");
        let a = &sol.ascent;
        assert_eq!(a.model, "simulation");
        let sim = a.simulation.as_ref().unwrap();
        assert!((sim.perigee_km - 185.0).abs() < 2.0 && (sim.apogee_km - 185.0).abs() < 2.0);
        // Точка выхода на орбиту в J2000 — на высоте орбиты, в плоскости опорной орбиты.
        let h = sol.parking_orbit.plane_normal;
        assert!((a.insertion_position_km.norm() - sol.parking_orbit.radius_km).abs() < 2.0);
        assert!(
            a.insertion_position_km.normalize().dot(&h).abs() < 2e-3,
            "{}",
            a.insertion_position_km.normalize().dot(&h)
        );
        assert!(a.insertion_velocity_km_s.normalize().dot(&h).abs() < 2e-3);
        assert!(sol.departure.coast_s.to_seconds() > a.duration_s);
        let r = refine(
            &almanac,
            &NumericalConfig::default(),
            &sol,
            Some(&t),
            None,
            None,
        )
        .unwrap()
        .unwrap();
        assert!(r.converged && r.landing.as_ref().unwrap().miss_km < 0.1);
        let path = r.path.unwrap();
        // В журнале траектории — точки подъёма из моделирования и события ступеней.
        let names: Vec<_> = path.events.iter().map(|e| e.name).collect();
        assert!(
            names.contains(&"booster_separation") && names.contains(&"max_q"),
            "{names:?}"
        );
        assert!(
            path.events.windows(2).all(|w| w[0].t_s <= w[1].t_s),
            "{names:?}"
        );
        let asc: Vec<_> = path
            .samples
            .iter()
            .filter(|s| s.phase == crate::trajlog::Phase::Ascent)
            .collect();
        assert!(asc.len() > 400);
        let first_park = path
            .samples
            .iter()
            .find(|s| s.phase == crate::trajlog::Phase::ParkingOrbit)
            .unwrap();
        let last = asc.last().unwrap();
        let d: f64 = (0..3)
            .map(|i| (first_park.r[i] - last.r[i]).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(d < 0.01, "ascent/parking joint {d} km");
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_refinement_lands_on_the_site() {
        let (almanac, catalog) = real();
        let ksc = catalog.get("ksc").unwrap();
        let day = Epoch::from_gregorian_utc_at_midnight(2027, 1, 1);
        for (lat, lon) in [(5.0, -10.0), (0.0, -30.0)] {
            let t = target(lat, lon);
            let rep = DirectTransfer::default()
                .solve(
                    &almanac,
                    ksc,
                    DateConstraint::ArrivalDay(day),
                    None,
                    Some(&t),
                )
                .unwrap();
            let sol = rep.solution.unwrap();
            let conic = sol.landing.unwrap();
            let r = refine(
                &almanac,
                &NumericalConfig::default(),
                &sol,
                Some(&t),
                None,
                None,
            )
            .unwrap()
            .unwrap();
            let l = r.landing.clone().expect("refined landing");
            assert!(r.converged, "{lat} {lon}: {r:?}");
            assert!(l.miss_km < 0.1, "miss {}", l.miss_km);
            assert!(l.conditions_met, "{:?}", l.unmet_reason);
            assert!(
                r.tli_delta_v_change_m_s.abs() < 40.0,
                "{}",
                r.tli_delta_v_change_m_s
            );
            // Порядок событий и разумные величины.
            assert!(
                r.arrival_epoch < l.doi_epoch
                    && l.doi_epoch < l.pdi_epoch
                    && l.pdi_epoch < l.touchdown_epoch
            );
            assert!(
                (l.pdi_altitude_km - 15.0).abs() < 4.0,
                "pdi {}",
                l.pdi_altitude_km
            );
            assert!(l.orbit_min_altitude_km > 90.0 && l.orbit_max_altitude_km < 110.0);
            assert!(
                l.orbit_max_altitude_km - l.orbit_min_altitude_km > 0.3,
                "mascons change the altitude"
            );
            assert!((l.doi_delta_v_m_s - conic.descent.doi_delta_v_m_s).abs() < 3.0);
            // Журнал траектории: события по порядку, этапы непрерывны, касание — в точке посадки.
            let path = r.path.as_ref().expect("trajectory path");
            let names: Vec<_> = path.events.iter().map(|e| e.name).collect();
            assert_eq!(
                names,
                [
                    "launch",
                    "parking_insertion",
                    "tli",
                    "loi",
                    "doi",
                    "pdi",
                    "braking_end",
                    "touchdown"
                ]
            );
            assert!(
                path.samples
                    .first()
                    .is_some_and(|s| s.phase == crate::trajlog::Phase::Ascent)
            );
            assert!((path.events[1].t_s - sol.ascent.duration_s).abs() < 1e-6);
            assert!(path.events.windows(2).all(|w| w[0].t_s <= w[1].t_s));
            assert!(path.samples.windows(2).all(|w| w[0].t_s <= w[1].t_s));
            let td = path.events.last().unwrap();
            let site = t.site.site.position_km().norm();
            let rm = td.r_moon.unwrap();
            assert!(((rm[0] * rm[0] + rm[1] * rm[1] + rm[2] * rm[2]).sqrt() - site).abs() < 0.01);
            let loi = &path.events[3];
            let dv = loi.delta_v.unwrap();
            assert!(
                ((dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]).sqrt() - r.loi_delta_v_km_s).abs()
                    < 1e-5
            );
            assert!(
                (l.braking.delta_v_m_s - conic.descent.powered_descent_delta_v_m_s).abs() < 40.0
            );
            // Касание — в пределах нескольких часов от расчёта по формулам.
            assert!(
                l.touchdown_shift_s.abs() < 4.0 * 3600.0,
                "{}",
                l.touchdown_shift_s
            );

            // Без масконов высота орбиты почти не меняется, а момент касания другой.
            let smooth = NumericalConfig {
                lunar_gravity_degree: 0,
                ..Default::default()
            };
            let r0 = refine(&almanac, &smooth, &sol, Some(&t), None, None)
                .unwrap()
                .unwrap();
            let l0 = r0.landing.unwrap();
            assert!(
                l0.orbit_max_altitude_km - l0.orbit_min_altitude_km
                    < l.orbit_max_altitude_km - l.orbit_min_altitude_km
            );
            assert!((l0.touchdown_epoch - l.touchdown_epoch).to_seconds().abs() > 0.1);
        }
    }
}
