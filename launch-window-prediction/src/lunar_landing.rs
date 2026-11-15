//! Посадка в заданную точку после прилёта к Луне.
//!
//! Связывает координаты точки посадки с датами миссии:
//!
//! 1. **Плоскость орбиты.** Аппарат подходит к Луне по гиперболе со скоростью v∞ и тормозит в
//!    периселении (LOI) на круговую орбиту высотой `lunar_orbit_altitude_km`. Плоскость гиперболы
//!    можно выбрать любой из содержащих направление v∞ — это прицеливание на подлёте (коррекция
//!    в пути в несколько м/с), поэтому выход на орбиту всегда компланарный. Плоскость выбирается
//!    один раз — при LOI — и дальше не меняется манёврами, поэтому поворот плоскости не нужен:
//!    она с самого начала нацелена так, чтобы в момент посадки пройти через точку.
//! 2. **Сжатие Луны (J2).** Неправильность поля тяготения поворачивает плоскость низкой орбиты
//!    вокруг оси Луны (до ~1.2°/сутки у экваториальной орбиты) и немного меняет период.
//!    Плоскость при LOI выбирается с поправкой на этот поворот к моменту посадки.
//! 3. **Спуск не мгновенный.** Импульс схода с орбиты (DOI) переводит аппарат на эллипс
//!    до перицентра `pdi_altitude_km`; полвитка пассивного полёта; затем торможение двигателем
//!    (PDI → касание), которое моделируется численно ([`simulate_descent`]): отсюда его длительность,
//!    дальность, расход топлива и нужная тяга. DOI выполняется за (180° + дуга торможения) до точки;
//!    моментом посадки считается касание.
//! 4. **Условия посадки.** В момент касания Солнце видно над рельефом и стоит в разрешённом
//!    диапазоне высот, Земля видна для связи. Иначе аппарат ждёт следующего витка
//!    (~2 ч), но DOI не позже `max_wait_in_orbit_h` после прилёта.
//!
//! Это быстрый расчёт по формулам для перебора вариантов. Найденный лучший вариант затем
//! уточняется численным интегрированием ([`crate::refine`]) с масконами и притяжением Земли и Солнца.

use anise::constants::frames::{MOON_J2000, MOON_ME_DE440_ME421_FRAME};
use anise::math::Vector3;
use anise::prelude::Almanac;
use anise::time::{Epoch, Unit};
use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};
use crate::runlog::{ser_epoch, ser_vec3};
use crate::sky::{MOON_RADIUS_KM, sky_at};
use crate::spacecraft::G0;
use crate::visibility::{BodyVisibility, SiteVisibility};

/// Безразмерный коэффициент сжатия Луны J2 (GRAIL GL0660B: C20 = −9.088e-5, J2 = −C20·√5).
pub const MOON_J2: f64 = 2.032e-4;
/// Опорный радиус модели гравитационного поля GRAIL, км.
pub const MOON_J2_REFERENCE_RADIUS_KM: f64 = 1738.0;
/// Аппарат по умолчанию для расчёта торможения, если `--spacecraft` не задан:
/// посадочная ступень Apollo LM с взлётной ступенью (масса, тяга, удельный импульс из пресета).
pub const APOLLO_LM_VEHICLE: DescentVehicle =
    DescentVehicle { initial_mass_kg: 16103.0, max_thrust_n: 46706.0, isp_s: 311.0 };
/// Угловая скорость вращения Луны, рад/с (сидерический период 27.321661 сут).
pub const MOON_ROTATION_RATE_RAD_S: f64 = 2.661_699_5e-6;

/// Условия посадки (секция `[landing_window]` в `data/transfer.toml`).
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct LandingWindowConfig {
    /// Солнце должно быть видно над рельефом (по правилам из `visibility.toml`).
    pub require_sun_visible: bool,
    /// Допустимая высота Солнца над ровным горизонтом в момент посадки, градусы.
    pub sun_elevation_min_deg: f64,
    pub sun_elevation_max_deg: f64,
    /// Земля должна быть видна над рельефом (связь при посадке).
    pub require_earth_visible: bool,
    /// Наибольшее время на окололунной орбите от прилёта до схода с орбиты (DOI), часы.
    pub max_wait_in_orbit_h: f64,
}

impl Default for LandingWindowConfig {
    fn default() -> Self {
        Self {
            require_sun_visible: true,
            sun_elevation_min_deg: 0.0,
            sun_elevation_max_deg: 20.0,
            require_earth_visible: true,
            max_wait_in_orbit_h: 48.0,
        }
    }
}

impl LandingWindowConfig {
    pub fn validate(&self) -> Result<()> {
        let bad = |m: &str| Err(Error::Config(format!("landing_window: {m}")));
        if !(self.sun_elevation_min_deg >= -10.0
            && self.sun_elevation_max_deg <= 90.0
            && self.sun_elevation_min_deg < self.sun_elevation_max_deg)
        {
            return bad("need -10 <= sun_elevation_min_deg < sun_elevation_max_deg <= 90");
        }
        if !(self.max_wait_in_orbit_h >= 0.0 && self.max_wait_in_orbit_h <= 24.0 * 30.0) {
            return bad("max_wait_in_orbit_h must be in [0, 720]");
        }
        Ok(())
    }
}

/// Параметры спуска с окололунной орбиты (секция `[descent]` в `data/transfer.toml`).
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct DescentConfig {
    /// Высота перицентра орбиты спуска, где начинается торможение (PDI), км. Apollo: ~15 км.
    pub pdi_altitude_km: f64,
    /// Наименьшее время на окололунной орбите от прилёта до DOI (проверка систем), часы.
    pub min_time_in_orbit_h: f64,
    /// Тяга при торможении в долях от максимальной. Apollo: основной участок на ~94% (9870 из 10500 фунтов).
    pub braking_throttle: f64,
    /// Добавка Δv к смоделированному торможению на заход и вертикальный участок у поверхности, м/с.
    /// Подобрана по Apollo: автоматический спуск LM требовал 6827 фут/с = 2081 м/с (NTRS 19700024568),
    /// моделирование торможения LM до остановки даёт ~1795 м/с.
    pub terminal_allowance_m_s: f64,
    /// Длительность захода и вертикального участка после торможения, с. Apollo: весь автоматический
    /// спуск 11 мин 20 с, торможение по моделированию ~8 мин 15 с — остаётся ~3 мин.
    pub terminal_duration_s: f64,
    /// Учитывать сжатие Луны при поиске по формулам (поворот плоскости орбиты и поправку периода).
    pub lunar_j2: bool,
}

impl Default for DescentConfig {
    fn default() -> Self {
        Self {
            pdi_altitude_km: 15.0,
            min_time_in_orbit_h: 2.0,
            braking_throttle: 0.94,
            terminal_allowance_m_s: 285.0,
            terminal_duration_s: 185.0,
            lunar_j2: true,
        }
    }
}

impl DescentConfig {
    pub fn validate(&self) -> Result<()> {
        let bad = |m: &str| Err(Error::Config(format!("descent: {m}")));
        if !(self.pdi_altitude_km >= 1.0 && self.pdi_altitude_km <= 100.0) {
            return bad("pdi_altitude_km must be in [1, 100]");
        }
        if !(self.min_time_in_orbit_h >= 0.0 && self.min_time_in_orbit_h <= 24.0 * 30.0) {
            return bad("min_time_in_orbit_h must be in [0, 720]");
        }
        if !(self.braking_throttle > 0.0 && self.braking_throttle <= 1.0) {
            return bad("braking_throttle must be in (0, 1]");
        }
        if !(self.terminal_allowance_m_s >= 0.0 && self.terminal_allowance_m_s <= 2000.0) {
            return bad("terminal_allowance_m_s must be in [0, 2000]");
        }
        if !(self.terminal_duration_s >= 0.0 && self.terminal_duration_s <= 3600.0) {
            return bad("terminal_duration_s must be in [0, 3600]");
        }
        Ok(())
    }
}

/// Двигатель и масса лунного аппарата — для моделирования торможения.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct DescentVehicle {
    /// Масса аппарата до его первого манёвра, кг.
    pub initial_mass_kg: f64,
    pub max_thrust_n: f64,
    pub isp_s: f64,
}

/// Всё, что нужно для расчёта спуска.
#[derive(Debug, Clone, PartialEq, Serialize, Default)]
pub struct DescentModel {
    pub config: DescentConfig,
    /// `None` — берётся Apollo LM ([`APOLLO_LM_VEHICLE`]).
    pub vehicle: Option<DescentVehicle>,
}

/// Смоделированное торможение двигателем от PDI до касания.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct DescentSim {
    pub duration_s: f64,
    /// Δv по формуле Циолковского (с гравитационными потерями), м/с.
    pub delta_v_m_s: f64,
    /// Тяга в долях от максимальной (меньше заданной, если двигатель избыточно мощный).
    pub throttle: f64,
    pub initial_mass_kg: f64,
    pub final_mass_kg: f64,
    /// Точка касания (в системе, где задано начальное состояние), км.
    #[serde(serialize_with = "ser_vec3")]
    pub touchdown_km: Vector3,
    /// Угол по поверхности от PDI до касания, рад.
    pub range_angle_rad: f64,
    pub range_km: f64,
    /// Точки торможения каждые 5 с: (секунды от начала, положение, скорость) — для журнала траектории.
    #[serde(skip)]
    pub samples: Vec<(f64, Vector3, Vector3)>,
}

pub const REASON_DESCENT_THRUST: &str =
    "powered descent is impossible: the engine thrust is too low to stop above the surface (raise thrust or pdi_altitude_km)";

/// Моделирует торможение двигателем от PDI до остановки над точкой.
///
/// Тяга постоянна (`throttle` от максимальной), как на основном участке торможения Apollo.
/// Направление тяги задаёт простое наведение: вертикальная составляющая выводит аппарат на
/// радиус `target_radius_km` с нулевой вертикальной скоростью к моменту остановки (квадратичный
/// закон, как в наведении Apollo), остальная тяга гасит горизонтальную скорость относительно
/// поверхности. Движение — в инерциальной системе с началом в центре Луны (точечная масса `mu`),
/// `omega` — вектор вращения Луны, рад/с. Интегрирование — метод Рунге — Кутты 4-го порядка
/// с шагом 0.5 с, масса убывает по расходу топлива.
pub fn simulate_descent(
    r0: &Vector3,
    v0: &Vector3,
    mass_kg: f64,
    vehicle: &DescentVehicle,
    throttle: f64,
    omega: &Vector3,
    mu: f64,
    target_radius_km: f64,
) -> std::result::Result<DescentSim, &'static str> {
    // Мощный двигатель гасит скорость быстрее, чем аппарат успевает снизиться: тогда тяга
    // уменьшается (дросселирование), пока торможение не приведёт к поверхности.
    let mut th = throttle;
    for _ in 0..60 {
        match descent_run(r0, v0, mass_kg, vehicle, th, omega, mu, target_radius_km) {
            Ok(sim) => return Ok(sim),
            Err(true) => th *= 0.93,
            Err(false) => return Err(REASON_DESCENT_THRUST),
        }
    }
    Err(REASON_DESCENT_THRUST)
}

/// Один прогон торможения с заданной тягой. `Err(true)` — остановились выше точки
/// (тяга избыточна), `Err(false)` — тяги не хватило.
#[allow(clippy::too_many_arguments)]
fn descent_run(
    r0: &Vector3,
    v0: &Vector3,
    mass_kg: f64,
    vehicle: &DescentVehicle,
    throttle: f64,
    omega: &Vector3,
    mu: f64,
    target_radius_km: f64,
) -> std::result::Result<DescentSim, bool> {
    let ve = G0 * vehicle.isp_s; // м/с
    let thrust = vehicle.max_thrust_n * throttle;
    let mdot = thrust / ve;
    // Ускорение (км/с²) и признак «тяги не хватает держать высоту».
    let f = |r: &Vector3, v: &Vector3, m: f64| -> (Vector3, bool) {
        let rn = r.norm();
        let up = r / rn;
        let v_rel = v - omega.cross(r);
        let v_up = v_rel.dot(&up);
        let v_h = v_rel - up * v_up;
        let vh = v_h.norm();
        let v_h_inertial = (v - up * v.dot(&up)).norm();
        let a_t = thrust / m / 1000.0;
        // Сколько ещё тормозить: время, за которое двигатель гасит горизонтальную скорость.
        let t_burn = (m / mdot) * (1.0 - (-vh * 1000.0 / ve).exp());
        let balance = mu / (rn * rn) - v_h_inertial * v_h_inertial / rn;
        let vertical = |t_go: f64| 6.0 * (target_radius_km - rn) / (t_go * t_go) - 4.0 * v_up / t_go + balance;
        // Часть тяги уходит на вертикаль, поэтому торможение дольше: уточняем время один раз.
        let first = vertical(t_burn + 5.0).clamp(-0.5 * a_t, 0.95 * a_t);
        let t_go = t_burn * a_t / (a_t * a_t - first * first).sqrt() + 5.0;
        let need = vertical(t_go);
        // Вниз двигатель толкает не более чем вполсилы: снижение в основном обеспечивает тяготение.
        let a_up = need.clamp(-0.5 * a_t, 0.95 * a_t);
        let a_h = (a_t * a_t - a_up * a_up).sqrt();
        let brake = if vh > 1e-9 { -v_h / vh * a_h } else { Vector3::zeros() };
        (-r * (mu / rn.powi(3)) + up * a_up + brake, need > 0.95 * a_t)
    };
    let h_vel = |r: &Vector3, v: &Vector3| {
        let up = r.normalize();
        let v_rel = v - omega.cross(r);
        v_rel - up * v_rel.dot(&up)
    };
    let dt = 0.5;
    let (mut r, mut v, mut m, mut t) = (*r0, *v0, mass_kg, 0.0);
    let mut saturated_s = 0.0;
    let mut vh = h_vel(&r, &v);
    let mut samples = vec![(0.0, r, v)];
    // До остановки: горизонтальная скорость относительно поверхности погашена (или сменила знак).
    while vh.norm() > 5e-4 {
        let (k1v, sat) = f(&r, &v, m);
        let (r2, v2) = (r + v * (dt / 2.0), v + k1v * (dt / 2.0));
        let k2v = f(&r2, &v2, m - mdot * dt / 2.0).0;
        let (r3, v3) = (r + v2 * (dt / 2.0), v + k2v * (dt / 2.0));
        let k3v = f(&r3, &v3, m - mdot * dt / 2.0).0;
        let (r4, v4) = (r + v3 * dt, v + k3v * dt);
        let k4v = f(&r4, &v4, m - mdot * dt).0;
        r += (v + v2 * 2.0 + v3 * 2.0 + v4) * (dt / 6.0);
        v += (k1v + k2v * 2.0 + k3v * 2.0 + k4v) * (dt / 6.0);
        m -= mdot * dt;
        t += dt;
        saturated_s += if sat { dt } else { 0.0 };
        if (t / dt).round() as u64 % 10 == 0 {
            samples.push((t, r, v));
        }
        let vh_new = h_vel(&r, &v);
        // Остаток меньше, чем двигатель гасит за один шаг, шагом уже не различить.
        if vh_new.dot(&vh) <= 0.0 || vh_new.norm() < thrust / m / 1000.0 * dt {
            break;
        }
        vh = vh_new;
        // Упали на поверхность, не погасив скорость; либо кончилась масса или время.
        if r.norm() < target_radius_km - 0.2 || m < 0.02 * mass_kg || t > 4.0 * 3600.0 {
            return Err(false);
        }
    }
    // Остановились заметно выше точки: не успели снизиться. Ниже или с нехваткой тяги на
    // удержание высоты — тяги мало.
    if r.norm() - target_radius_km > 0.3 {
        return Err(true);
    }
    if r.norm() - target_radius_km < -0.3 || saturated_s > 0.5 * t {
        return Err(false);
    }
    let range_angle = r0.cross(&r).norm().atan2(r0.dot(&r));
    Ok(DescentSim {
        duration_s: t,
        delta_v_m_s: ve * (mass_kg / m).ln(),
        throttle,
        initial_mass_kg: mass_kg,
        final_mass_kg: m,
        touchdown_km: r,
        range_angle_rad: range_angle,
        range_km: range_angle * target_radius_km,
        samples: {
            if samples.last().is_none_or(|l| l.0 < t) {
                samples.push((t, r, v));
            }
            samples
        },
    })
}

/// Рассчитанный спуск (не зависит от даты).
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct DescentPlan {
    /// Импульс схода с круговой орбиты на эллипс спуска, м/с.
    pub doi_delta_v_m_s: f64,
    /// Δv торможения от PDI до касания по моделированию, м/с.
    pub powered_descent_delta_v_m_s: f64,
    /// Добавка на заход и вертикальный участок, м/с.
    pub terminal_allowance_m_s: f64,
    /// Весь спуск: DOI + торможение + добавка, м/с.
    pub total_delta_v_m_s: f64,
    /// Полёт по эллипсу от DOI до PDI (половина витка), с.
    pub coast_s: f64,
    /// Торможение двигателем до остановки над точкой (по моделированию), с.
    pub braking_duration_s: f64,
    /// Всё от PDI до касания: торможение + заход и вертикальный участок, с.
    pub powered_descent_duration_s: f64,
    /// Тяга при торможении в долях от максимальной.
    pub throttle: f64,
    /// Масса в начале торможения, кг.
    pub pdi_mass_kg: f64,
    /// Чьи параметры двигателя: "spacecraft" или "apollo_lm".
    pub vehicle_source: &'static str,
    /// Скорость в PDI, км/с.
    pub pdi_speed_km_s: f64,
    /// Дальность торможения по поверхности, км.
    pub powered_descent_range_km: f64,
    /// Угловая дальность торможения, рад.
    #[serde(skip)]
    pub range_angle_rad: f64,
}

impl DescentPlan {
    pub fn total_s(&self) -> f64 {
        self.coast_s + self.powered_descent_duration_s
    }
}

impl DescentModel {
    pub fn vehicle(&self) -> (DescentVehicle, &'static str) {
        self.vehicle.map_or((APOLLO_LM_VEHICLE, "apollo_lm"), |v| (v, "spacecraft"))
    }

    /// Масса в начале торможения после манёвров на `spent_delta_v_m_s`.
    pub fn mass_after(&self, spent_delta_v_m_s: f64) -> f64 {
        let (v, _) = self.vehicle();
        v.initial_mass_kg * (-spent_delta_v_m_s / (G0 * v.isp_s)).exp()
    }

    /// Номинальный спуск с круговой орбиты радиуса `r_orbit_km` на радиус `target_radius_km`.
    /// `prior_delta_v_m_s` — Δv, уже израсходованный аппаратом до спуска (для массы в PDI).
    pub fn plan(
        &self,
        mu_moon: f64,
        r_orbit_km: f64,
        prior_delta_v_m_s: f64,
        target_radius_km: f64,
    ) -> Result<std::result::Result<DescentPlan, &'static str>> {
        let c = &self.config;
        let r_p = MOON_RADIUS_KM + c.pdi_altitude_km;
        if !(r_p < r_orbit_km) {
            return Err(Error::Config(format!(
                "descent: pdi_altitude_km ({}) must be below the lunar orbit altitude ({:.1} km)",
                c.pdi_altitude_km,
                r_orbit_km - MOON_RADIUS_KM
            )));
        }
        if !(target_radius_km < r_p) {
            return Err(Error::Config(format!(
                "descent: the landing site ({:.1} km above the reference sphere) is above pdi_altitude_km ({})",
                target_radius_km - MOON_RADIUS_KM,
                c.pdi_altitude_km
            )));
        }
        let a = 0.5 * (r_orbit_km + r_p);
        let v_c = (mu_moon / r_orbit_km).sqrt();
        let v_apo = (mu_moon * (2.0 / r_orbit_km - 1.0 / a)).sqrt();
        let v_p = (mu_moon * (2.0 / r_p - 1.0 / a)).sqrt();
        let doi = (v_c - v_apo) * 1000.0;
        let (vehicle, source) = self.vehicle();
        let m_pdi = self.mass_after(prior_delta_v_m_s + doi);
        let sim = match simulate_descent(
            &Vector3::new(r_p, 0.0, 0.0),
            &Vector3::new(0.0, v_p, 0.0),
            m_pdi,
            &vehicle,
            c.braking_throttle,
            &Vector3::zeros(),
            mu_moon,
            target_radius_km,
        ) {
            Ok(s) => s,
            Err(reason) => return Ok(Err(reason)),
        };
        Ok(Ok(DescentPlan {
            doi_delta_v_m_s: doi,
            powered_descent_delta_v_m_s: sim.delta_v_m_s,
            terminal_allowance_m_s: c.terminal_allowance_m_s,
            total_delta_v_m_s: doi + sim.delta_v_m_s + c.terminal_allowance_m_s,
            coast_s: std::f64::consts::PI * (a.powi(3) / mu_moon).sqrt(),
            braking_duration_s: sim.duration_s,
            powered_descent_duration_s: sim.duration_s + c.terminal_duration_s,
            throttle: sim.throttle,
            pdi_mass_kg: m_pdi,
            vehicle_source: source,
            pdi_speed_km_s: v_p,
            powered_descent_range_km: sim.range_km,
            range_angle_rad: sim.range_angle_rad,
        }))
    }
}

/// Точка посадки вместе с условиями и моделью спуска.
#[derive(Debug, Clone)]
pub struct LandingTarget {
    /// Точка, её рельеф и правила видимости.
    pub site: SiteVisibility,
    pub window: LandingWindowConfig,
    pub descent: DescentModel,
}

/// Найденная посадка.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct LandingInfo {
    /// Момент касания поверхности.
    #[serde(serialize_with = "ser_epoch")]
    pub epoch: Epoch,
    /// Сход с орбиты (DOI).
    #[serde(serialize_with = "ser_epoch")]
    pub doi_epoch: Epoch,
    /// Начало торможения (PDI).
    #[serde(serialize_with = "ser_epoch")]
    pub pdi_epoch: Epoch,
    /// Время на окололунной орбите от прилёта до DOI, с.
    pub wait_after_arrival_s: f64,
    /// Номер витка, на котором выполняется спуск (0 — первая возможность после прилёта).
    pub orbit_pass: u32,
    pub lunar_orbit_period_s: f64,
    /// Наклонение окололунной орбиты к лунному экватору, градусы.
    pub lunar_orbit_inclination_deg: f64,
    /// Нормаль плоскости орбиты при выходе на неё (LOI), J2000.
    #[serde(serialize_with = "ser_vec3")]
    pub lunar_orbit_normal_j2000: Vector3,
    /// Нормаль плоскости в момент посадки (после поворота из-за J2), J2000.
    #[serde(serialize_with = "ser_vec3")]
    pub lunar_orbit_normal_at_landing_j2000: Vector3,
    /// Поворот плоскости орбиты вокруг оси Луны от LOI до посадки (J2), градусы.
    pub nodal_precession_deg: f64,
    /// Угол по орбите от периселения (LOI) до точки посадки, градусы.
    pub angle_periselene_to_site_deg: f64,
    /// Эксцентриситет гиперболы подлёта.
    pub approach_eccentricity: f64,
    pub descent: DescentPlan,
    pub sun: BodyVisibility,
    pub earth: BodyVisibility,
}

/// Почему посадки нет.
pub const REASON_LANDING_SUN: &str =
    "no landing opportunity within max_wait_in_orbit_h: the Sun is hidden or outside the allowed elevation at the site";
pub const REASON_LANDING_EARTH: &str =
    "no landing opportunity within max_wait_in_orbit_h: the Earth is not visible from the site";
pub const REASON_LANDING_BOTH: &str =
    "no landing opportunity within max_wait_in_orbit_h: neither lighting nor Earth visibility conditions are met at the site";
pub const REASON_LANDING_DAY: &str = "no landing opportunity at the site within the requested day";
/// Уточнённый расчёт: аппарат прилетел раньше, и до схода с орбиты пришлось бы ждать дольше лимита.
pub const REASON_LANDING_WAIT: &str =
    "the landing pass needs a longer wait in lunar orbit than max_wait_in_orbit_h allows";

/// Поворот вектора `v` вокруг единичной оси `k` на угол `a` (формула Родрига).
pub(crate) fn rotate_about(v: Vector3, k: Vector3, a: f64) -> Vector3 {
    let (s, c) = a.sin_cos();
    v * c + k.cross(&v) * s + k * (k.dot(&v) * (1.0 - c))
}

/// Ищет первую посадку после прилёта в момент `t_arrival` со скоростью подлёта `v_inf`
/// (относительно Луны, J2000). `window` — допустимый интервал момента касания (например,
/// заданные сутки). `mu_moon` — км³/с², `r_orbit_km` — радиус окололунной орбиты,
/// `prior_delta_v_m_s` — Δv, израсходованный аппаратом до спуска.
#[allow(clippy::too_many_arguments)]
pub fn find_landing(
    almanac: &Almanac,
    target: &LandingTarget,
    t_arrival: Epoch,
    v_inf: Vector3,
    mu_moon: f64,
    r_orbit_km: f64,
    prior_delta_v_m_s: f64,
    window: Option<(Epoch, Epoch)>,
) -> Result<std::result::Result<LandingInfo, &'static str>> {
    let mut passes = LandingPasses::default();
    find_landing_in(&mut passes, almanac, target, t_arrival, v_inf, mu_moon, r_orbit_km, prior_delta_v_m_s, window)
}

/// Витки, уже посчитанные для одного прилёта (момент, скорость подлёта, орбита, израсходованный
/// Δv): план спуска, геометрия каждого витка и видимость в момент касания от суток не зависят.
/// Экономный режим просматривает один и тот же прилёт для нескольких соседних суток — с
/// [`LandingCache`] витки считаются один раз, а для каждых суток применяется только их окно.
#[derive(Debug, Default)]
pub struct LandingPasses {
    plan: Option<std::result::Result<DescentPlan, &'static str>>,
    pole: Option<Vector3>,
    /// Геометрия витков по порядку; `None` — на этом витке ожидание уже больше лимита (дальше не считается).
    geo: Vec<Option<PassGeo>>,
    /// Видимость в момент касания на витке: (Солнце, Земля), посчитанная по мере надобности.
    vis: Vec<Option<(BodyVisibility, BodyVisibility)>>,
}

#[derive(Debug, Clone, Copy)]
struct PassGeo {
    t: Epoch,
    h: Vector3,
    theta_site: f64,
    d_node: f64,
    t_doi: Epoch,
    wait_s: f64,
}

/// Общий кэш [`LandingPasses`] по прилётам (ключ — точные значения входных данных).
#[derive(Debug, Default)]
pub struct LandingCache {
    entries: std::sync::Mutex<std::collections::HashMap<[u64; 7], std::sync::Arc<std::sync::Mutex<LandingPasses>>>>,
}

impl LandingCache {
    /// [`find_landing`] с запоминанием витков этого прилёта.
    #[allow(clippy::too_many_arguments)]
    pub fn find(
        &self,
        almanac: &Almanac,
        target: &LandingTarget,
        t_arrival: Epoch,
        v_inf: Vector3,
        mu_moon: f64,
        r_orbit_km: f64,
        prior_delta_v_m_s: f64,
        window: Option<(Epoch, Epoch)>,
    ) -> Result<std::result::Result<LandingInfo, &'static str>> {
        let key = [
            t_arrival.to_tai_duration().total_nanoseconds() as u64,
            v_inf.x.to_bits(),
            v_inf.y.to_bits(),
            v_inf.z.to_bits(),
            mu_moon.to_bits(),
            r_orbit_km.to_bits(),
            prior_delta_v_m_s.to_bits(),
        ];
        let entry = self.entries.lock().unwrap_or_else(|e| e.into_inner()).entry(key).or_default().clone();
        // Без ожидания: если витки этого прилёта сейчас считает другой поток (или этот же, выше по
        // стеку, — rayon может подсунуть задачу ждущему потоку), считаем без кэша. Ждать здесь
        // нельзя — это может стать взаимной блокировкой.
        let mut passes = match entry.try_lock() {
            Ok(g) => g,
            Err(std::sync::TryLockError::Poisoned(e)) => e.into_inner(),
            Err(std::sync::TryLockError::WouldBlock) => {
                let mut fresh = LandingPasses::default();
                return find_landing_in(&mut fresh, almanac, target, t_arrival, v_inf, mu_moon, r_orbit_km, prior_delta_v_m_s, window);
            }
        };
        find_landing_in(&mut passes, almanac, target, t_arrival, v_inf, mu_moon, r_orbit_km, prior_delta_v_m_s, window)
    }
}

#[allow(clippy::too_many_arguments)]
fn find_landing_in(
    memo: &mut LandingPasses,
    almanac: &Almanac,
    target: &LandingTarget,
    t_arrival: Epoch,
    v_inf: Vector3,
    mu_moon: f64,
    r_orbit_km: f64,
    prior_delta_v_m_s: f64,
    window: Option<(Epoch, Epoch)>,
) -> Result<std::result::Result<LandingInfo, &'static str>> {
    let tau = std::f64::consts::TAU;
    let v_hat = v_inf.normalize();
    let vinf = v_inf.norm();
    let ecc = 1.0 + r_orbit_km * vinf * vinf / mu_moon;
    let nu_inf = (-1.0 / ecc).acos();
    let n = (mu_moon / r_orbit_km.powi(3)).sqrt();
    let period_s = tau / n;
    let w = &target.window;
    let dc = &target.descent.config;
    let max_wait_s = w.max_wait_in_orbit_h * 3600.0;
    let min_wait_s = dc.min_time_in_orbit_h * 3600.0;
    let plan = match memo.plan {
        Some(p) => p,
        None => {
            let p = target.descent.plan(mu_moon, r_orbit_km, prior_delta_v_m_s, target.site.site.position_km().norm())?;
            memo.plan = Some(p);
            p
        }
    };
    let plan = match plan {
        Ok(p) => p,
        Err(reason) => return Ok(Err(reason)),
    };
    let descent_s = plan.total_s();
    let site_fixed = target.site.site.up();
    let j2k = if dc.lunar_j2 { MOON_J2 * (MOON_J2_REFERENCE_RADIUS_KM / r_orbit_km).powi(2) } else { 0.0 };

    let moon_dcm = |t: Epoch| {
        almanac
            .rotate(MOON_ME_DE440_ME421_FRAME, MOON_J2000, t)
            .map_err(|e| Error::Almanac(format!("Moon orientation at {t}: {e}")))
    };
    let pole = match memo.pole {
        Some(p) => p,
        None => {
            let p = moon_dcm(t_arrival)? * Vector3::new(0.0, 0.0, 1.0);
            memo.pole = Some(p);
            p
        }
    };

    // Геометрия витка `k` (момент касания уточняется итерациями: точка поворачивается вместе с
    // Луной, а плоскость орбиты — из-за J2; скорость поворота зависит от наклонения).
    let pass_geo = |k: u32| -> Result<PassGeo> {
        let mut t = t_arrival + (min_wait_s + k as f64 * period_s + descent_s) * Unit::Second;
        let mut cos_i = 0.0;
        let mut geo: Option<PassGeo> = None;
        for _ in 0..12 {
            let s = moon_dcm(t)? * site_fixed;
            let dt = (t - t_arrival).to_seconds();
            let node_rate = -1.5 * n * j2k * cos_i;
            let n_eff = n * (1.0 + 0.75 * j2k * (8.0 * cos_i * cos_i - 2.0));
            let d_node = node_rate * dt;
            // Точка в «плоскости при LOI»: убираем будущий поворот плоскости.
            let s_eff = rotate_about(s, pole, -d_node);
            let mut h0 = v_hat.cross(&s_eff);
            if h0.norm() < 1e-9 {
                // Точка на линии подлёта: годится любая плоскость, берём ближайшую к экватору.
                h0 = v_hat.cross(&pole);
                if h0.norm() < 1e-9 {
                    h0 = v_hat.cross(&Vector3::new(1.0, 0.0, 0.0));
                }
            }
            let mut best: Option<(Vector3, f64, f64)> = None;
            for sign in [1.0, -1.0] {
                let h = h0.normalize() * sign;
                let incoming = -v_hat;
                let p = incoming * nu_inf.cos() + h.cross(&incoming) * nu_inf.sin();
                let theta_site = p.cross(&s_eff).dot(&h).atan2(p.dot(&s_eff)).rem_euclid(tau);
                // DOI за полвитка эллипса и дугу торможения до точки.
                let mut theta_doi = theta_site - std::f64::consts::PI - plan.range_angle_rad;
                while theta_doi < n_eff * min_wait_s {
                    theta_doi += tau;
                }
                theta_doi += tau * k as f64;
                if best.is_none_or(|(_, _, b)| theta_doi < b) {
                    best = Some((h, theta_site, theta_doi));
                }
            }
            let (h, theta_site, theta_doi) = best.expect("two candidates");
            cos_i = h.dot(&pole);
            let t_doi = t_arrival + (theta_doi / n_eff) * Unit::Second;
            let t_new = t_doi + descent_s * Unit::Second;
            let converged = (t_new - t).to_seconds().abs() < 0.5;
            t = t_new;
            geo = Some(PassGeo { t, h, theta_site, d_node, t_doi, wait_s: (t_doi - t_arrival).to_seconds() });
            if converged {
                break;
            }
        }
        Ok(geo.expect("at least one iteration"))
    };

    let (mut sun_fail, mut earth_fail, mut any_pass) = (0u32, 0u32, false);
    let max_passes = ((max_wait_s - min_wait_s).max(0.0) / period_s).ceil() as u32 + 1;
    for k in 0..=max_passes {
        let ki = k as usize;
        if memo.geo.len() <= ki {
            // Витки считаются по порядку; после витка сверх лимита ожидания — не считаются.
            if memo.geo.last().is_some_and(|g| g.is_none()) {
                break;
            }
            let g = pass_geo(k)?;
            memo.geo.push((g.wait_s <= max_wait_s).then_some(g));
            memo.vis.push(None);
        }
        // Ожидание до схода с орбиты больше лимита: дальше витки ещё позже.
        let Some(g) = memo.geo[ki] else { break };
        let t = g.t;
        if let Some((lo, hi)) = window {
            if t < lo {
                continue;
            }
            if t >= hi {
                return Ok(Err(if any_pass { worst(sun_fail, earth_fail) } else { REASON_LANDING_DAY }));
            }
        }
        any_pass = true;

        let (sun, earth) = match memo.vis[ki] {
            Some(v) => v,
            None => {
                let sky = sky_at(almanac, &target.site.site, t)?;
                let vis = target.site.evaluate(&sky);
                memo.vis[ki] = Some((vis.sun, vis.earth));
                (vis.sun, vis.earth)
            }
        };
        // Пределы высоты Солнца действуют только когда Солнце требуется (без него - посадка в любое время).
        let sun_ok = !w.require_sun_visible
            || (sun.visible && sun.elevation_deg >= w.sun_elevation_min_deg && sun.elevation_deg <= w.sun_elevation_max_deg);
        let earth_ok = !w.require_earth_visible || earth.visible;
        if sun_ok && earth_ok {
            return Ok(Ok(LandingInfo {
                epoch: t,
                doi_epoch: g.t_doi,
                pdi_epoch: g.t_doi + plan.coast_s * Unit::Second,
                wait_after_arrival_s: g.wait_s,
                orbit_pass: k,
                lunar_orbit_period_s: period_s,
                lunar_orbit_inclination_deg: g.h.dot(&pole).clamp(-1.0, 1.0).acos().to_degrees(),
                lunar_orbit_normal_j2000: g.h,
                lunar_orbit_normal_at_landing_j2000: rotate_about(g.h, pole, g.d_node),
                nodal_precession_deg: g.d_node.to_degrees(),
                angle_periselene_to_site_deg: g.theta_site.to_degrees(),
                approach_eccentricity: ecc,
                descent: plan,
                sun,
                earth,
            }));
        }
        sun_fail += u32::from(!sun_ok);
        earth_fail += u32::from(!earth_ok);
    }
    Ok(Err(if !any_pass && window.is_some() { REASON_LANDING_DAY } else { worst(sun_fail, earth_fail) }))
}

fn worst(sun_fail: u32, earth_fail: u32) -> &'static str {
    match (sun_fail > 0, earth_fail > 0) {
        (true, true) => REASON_LANDING_BOTH,
        (false, true) => REASON_LANDING_EARTH,
        _ => REASON_LANDING_SUN,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ephemeris::{AlmanacSource, load_almanac_blocking};
    use crate::terrain::TerrainSet;
    use crate::visibility::VisibilityConfig;

    #[test]
    fn config_validation() {
        assert!(LandingWindowConfig::default().validate().is_ok());
        let bad = LandingWindowConfig { sun_elevation_min_deg: 30.0, ..Default::default() };
        assert!(bad.validate().is_err());
        let bad = LandingWindowConfig { max_wait_in_orbit_h: -1.0, ..Default::default() };
        assert!(bad.validate().is_err());
        assert!(DescentConfig::default().validate().is_ok());
        assert!(DescentConfig { terminal_allowance_m_s: -1.0, ..Default::default() }.validate().is_err());
        assert!(DescentConfig { pdi_altitude_km: 0.0, ..Default::default() }.validate().is_err());
    }

    const MU: f64 = 4902.8;
    const R_LLO: f64 = MOON_RADIUS_KM + 100.0;

    #[test]
    fn descent_simulation_matches_apollo() {
        // Без аппарата берётся Apollo LM: сход 100 → 15 км ≈ 20 м/с, полвитка ≈ 57 мин.
        let p = DescentModel::default().plan(MU, R_LLO, 0.0, MOON_RADIUS_KM).unwrap().unwrap();
        assert_eq!(p.vehicle_source, "apollo_lm");
        assert!((18.0..23.0).contains(&p.doi_delta_v_m_s), "{}", p.doi_delta_v_m_s);
        assert!((3300.0..3500.0).contains(&p.coast_s), "{}", p.coast_s);
        // Торможение: Apollo — 12-12.7 мин, ~480 км, автоматический спуск 2081 м/с.
        assert!((450.0..560.0).contains(&p.braking_duration_s), "{}", p.braking_duration_s);
        assert!((p.powered_descent_duration_s - 680.0).abs() < 30.0, "{}", p.powered_descent_duration_s);
        assert!((p.powered_descent_range_km - 481.0).abs() < 50.0, "{}", p.powered_descent_range_km);
        assert!((p.total_delta_v_m_s - p.doi_delta_v_m_s - 2081.0).abs() < 30.0, "{}", p.total_delta_v_m_s);
        assert_eq!(p.throttle, 0.94);
        assert!((p.total_delta_v_m_s - p.doi_delta_v_m_s - p.powered_descent_delta_v_m_s - 285.0).abs() < 1e-9);

        // Скорость гасится ровно на заданном радиусе, топливо расходуется по Циолковскому.
        let v = APOLLO_LM_VEHICLE;
        let r0 = Vector3::new(MOON_RADIUS_KM + 15.0, 0.0, 0.0);
        let sim = simulate_descent(&r0, &Vector3::new(0.0, 1.69, 0.0), 15000.0, &v, 0.94, &Vector3::zeros(), MU, MOON_RADIUS_KM - 2.0)
            .unwrap();
        assert!((sim.touchdown_km.norm() - (MOON_RADIUS_KM - 2.0)).abs() < 0.05, "{}", sim.touchdown_km.norm());
        let burned = v.max_thrust_n * sim.throttle / (G0 * v.isp_s) * sim.duration_s;
        assert!((sim.initial_mass_kg - sim.final_mass_kg - burned).abs() < 1.0);
        // Гравитационные потери: Δv больше погашенной скорости.
        assert!(sim.delta_v_m_s > 1690.0);

        // После торможения у Луны аппарат легче: торможение короче и потери меньше.
        let after_loi = DescentModel::default().plan(MU, R_LLO, 900.0, MOON_RADIUS_KM).unwrap().unwrap();
        assert!(after_loi.pdi_mass_kg < p.pdi_mass_kg);
        assert!(after_loi.braking_duration_s < p.braking_duration_s);
        assert!(after_loi.powered_descent_delta_v_m_s < p.powered_descent_delta_v_m_s);

        // Слабый двигатель: чем меньше тяговооружённость, тем больше потери; совсем слабый не сядет.
        let weak = DescentModel {
            vehicle: Some(DescentVehicle { max_thrust_n: 30000.0, ..APOLLO_LM_VEHICLE }),
            ..Default::default()
        };
        let pw = weak.plan(MU, R_LLO, 0.0, MOON_RADIUS_KM).unwrap().unwrap();
        assert!(pw.powered_descent_delta_v_m_s > p.powered_descent_delta_v_m_s + 20.0);
        let too_weak = DescentModel {
            vehicle: Some(DescentVehicle { max_thrust_n: 6000.0, ..APOLLO_LM_VEHICLE }),
            ..Default::default()
        };
        assert_eq!(too_weak.plan(MU, R_LLO, 0.0, MOON_RADIUS_KM).unwrap(), Err(REASON_DESCENT_THRUST));
        // Перицентр выше орбиты и точка выше перицентра — ошибки настроек.
        let high = DescentModel { config: DescentConfig { pdi_altitude_km: 100.0, ..Default::default() }, ..Default::default() };
        assert!(high.plan(MU, MOON_RADIUS_KM + 50.0, 0.0, MOON_RADIUS_KM).is_err());
        assert!(DescentModel::default().plan(MU, R_LLO, 0.0, MOON_RADIUS_KM + 16.0).is_err());
    }

    #[test]
    fn rodrigues_rotation() {
        let z = Vector3::new(0.0, 0.0, 1.0);
        let r = rotate_about(Vector3::new(1.0, 0.0, 0.0), z, std::f64::consts::FRAC_PI_2);
        assert!((r - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-12);
    }

    fn real_almanac() -> Almanac {
        load_almanac_blocking(AlmanacSource::LocalDir("data/kernels".into())).unwrap()
    }

    fn flat_target(lat: f64, lon: f64, window: LandingWindowConfig) -> LandingTarget {
        let site = SiteVisibility::new(None, lat, lon, 0.0, VisibilityConfig::default()).unwrap();
        LandingTarget { site, window, descent: DescentModel::default() }
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_landing_cache_gives_the_same_answers_for_every_day() {
        let almanac = real_almanac();
        let t_a = Epoch::from_gregorian_utc_hms(2026, 12, 3, 19, 42, 0);
        let v_inf = Vector3::new(0.3, -1.1, 0.5);
        let target = flat_target(10.0, 30.0, LandingWindowConfig { max_wait_in_orbit_h: 72.0, ..Default::default() });
        let cache = LandingCache::default();
        let day0 = Epoch::from_gregorian_utc_at_midnight(2026, 12, 3);
        // Сутки в разном порядке (кэш заполняется не по порядку) и без окна.
        let mut windows: Vec<Option<(Epoch, Epoch)>> =
            [3, 0, 1, 4, 2, 5].iter().map(|&d| Some((day0 + d * Unit::Day, day0 + (d + 1) * Unit::Day))).collect();
        windows.push(None);
        for w in windows {
            let plain = find_landing(&almanac, &target, t_a, v_inf, 4902.8, 1837.4, 900.0, w).unwrap();
            let cached = cache.find(&almanac, &target, t_a, v_inf, 4902.8, 1837.4, 900.0, w).unwrap();
            assert_eq!(plain, cached, "{w:?}");
        }
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_landing_is_over_the_site_and_meets_conditions() {
        let almanac = real_almanac();
        let t_a = Epoch::from_gregorian_utc_hms(2026, 12, 3, 19, 42, 0);
        let v_inf = Vector3::new(0.3, -1.1, 0.5);
        // Без условий: спуск на первой возможности — не позже витка после минимального времени на орбите.
        let free = LandingWindowConfig {
            require_sun_visible: false,
            require_earth_visible: false,
            sun_elevation_min_deg: -10.0,
            sun_elevation_max_deg: 90.0,
            ..Default::default()
        };
        let target = flat_target(10.0, 30.0, free);
        let r = 1837.4;
        let info = find_landing(&almanac, &target, t_a, v_inf, 4902.8, r, 0.0, None).unwrap().unwrap();
        assert_eq!(info.orbit_pass, 0);
        let min_wait = DescentConfig::default().min_time_in_orbit_h * 3600.0;
        assert!(info.wait_after_arrival_s >= min_wait - 1.0);
        assert!(info.wait_after_arrival_s <= min_wait + info.lunar_orbit_period_s * 1.01);
        // Спуск занимает время: DOI → PDI (полвитка эллипса) → касание.
        let d = info.descent;
        assert!(((info.pdi_epoch - info.doi_epoch).to_seconds() - d.coast_s).abs() < 1e-3);
        assert!(((info.epoch - info.pdi_epoch).to_seconds() - d.powered_descent_duration_s).abs() < 1.0);
        assert!((info.epoch - t_a).to_seconds() > 3600.0 * 2.0 + d.coast_s);
        // В момент касания точка лежит в плоскости орбиты, повёрнутой за счёт J2.
        let dcm = almanac.rotate(MOON_ME_DE440_ME421_FRAME, MOON_J2000, info.epoch).unwrap();
        let s = dcm * target.site.site.up();
        assert!(s.dot(&info.lunar_orbit_normal_at_landing_j2000).abs() < 2e-3, "{}", s.dot(&info.lunar_orbit_normal_at_landing_j2000));
        assert!(info.nodal_precession_deg.abs() > 0.0 && info.nodal_precession_deg.abs() < 1.0);

        // С условиями найденная посадка им удовлетворяет.
        let target = flat_target(10.0, 30.0, LandingWindowConfig { max_wait_in_orbit_h: 400.0, ..Default::default() });
        let info = find_landing(&almanac, &target, t_a, v_inf, 4902.8, r, 0.0, None).unwrap().unwrap();
        assert!(info.sun.visible && info.earth.visible);
        assert!(info.sun.elevation_deg <= 20.0);
        // Выход на орбиту компланарный (плоскость содержит v∞), и на любом витке ожидания
        // точка в момент касания лежит в плоскости, выбранной при LOI и повёрнутой J2.
        assert!(info.lunar_orbit_normal_j2000.dot(&v_inf.normalize()).abs() < 1e-9);
        let dcm = almanac.rotate(MOON_ME_DE440_ME421_FRAME, MOON_J2000, info.epoch).unwrap();
        let s = dcm * target.site.site.up();
        assert!(s.dot(&info.lunar_orbit_normal_at_landing_j2000).abs() < 2e-3, "pass {}", info.orbit_pass);

        // Без J2 плоскость не поворачивается.
        let mut no_j2 = flat_target(10.0, 30.0, free_window());
        no_j2.descent.config.lunar_j2 = false;
        let info = find_landing(&almanac, &no_j2, t_a, v_inf, 4902.8, r, 0.0, None).unwrap().unwrap();
        assert_eq!(info.nodal_precession_deg, 0.0);
        assert!((info.lunar_orbit_normal_at_landing_j2000 - info.lunar_orbit_normal_j2000).norm() < 1e-12);
    }

    fn free_window() -> LandingWindowConfig {
        LandingWindowConfig {
            require_sun_visible: false,
            require_earth_visible: false,
            sun_elevation_min_deg: -10.0,
            sun_elevation_max_deg: 90.0,
            ..Default::default()
        }
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_far_side_never_sees_earth() {
        let almanac = real_almanac();
        let t_a = Epoch::from_gregorian_utc_hms(2026, 12, 3, 19, 42, 0);
        let target = flat_target(0.0, 180.0, LandingWindowConfig::default());
        let r = find_landing(&almanac, &target, t_a, Vector3::new(0.3, -1.1, 0.5), 4902.8, 1837.4, 0.0, None).unwrap();
        assert!(matches!(r, Err(REASON_LANDING_EARTH) | Err(REASON_LANDING_BOTH)), "{r:?}");
    }

    #[test]
    #[ignore = "needs kernels and LOLA maps"]
    fn real_terrain_site_uses_horizon_mask() {
        let almanac = real_almanac();
        let terrain = TerrainSet::open_dir("data/terrain").unwrap();
        let site = SiteVisibility::new(Some(&terrain), -89.484, 29.38, 0.0, VisibilityConfig::default()).unwrap();
        let target = LandingTarget {
            site,
            window: LandingWindowConfig { max_wait_in_orbit_h: 24.0 * 20.0, ..Default::default() },
            descent: DescentModel::default(),
        };
        let t_a = Epoch::from_gregorian_utc_hms(2027, 1, 1, 21, 7, 0);
        let r = find_landing(&almanac, &target, t_a, Vector3::new(0.3, -1.1, 0.5), 4902.8, 1837.4, 0.0, None).unwrap();
        if let Ok(info) = r {
            assert!(info.sun.visible && info.earth.visible);
            // На этом месте Солнце и Земля скрыты первые трое суток после прилёта.
            assert!(info.wait_after_arrival_s > 72.0 * 3600.0, "wait {}", info.wait_after_arrival_s);
        }
    }
}
