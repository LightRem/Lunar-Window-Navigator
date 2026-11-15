//! Подъём ракеты от стартового стола до опорной орбиты — численное моделирование.
//!
//! Модель:
//! - **Земля:** эллипсоид WGS84, притяжение с учётом сжатия (J2), вращение Земли; атмосфера
//!   вращается вместе с Землёй.
//! - **Атмосфера:** стандартная атмосфера США 1976 г. (до 86 км — по слоям, выше — по таблице
//!   плотности до 1000 км). Сопротивление: коэффициент Cd в зависимости от числа Маха,
//!   площадь — сечение центрального блока и ускорителей.
//! - **Ракета:** ускорители (твердотопливные — с профилем тяги по времени) и ступени по порядку;
//!   тяга зависит от давления воздуха (F = F_vac − A_e·p), расход постоянен; сброс отработавших
//!   ступеней и головного обтекателя. Разгонная ступень (та же, что разгоняет к Луне) может
//!   довыводить на орбиту — тогда остаток её топлива идёт на разгон к Луне.
//! - **Наведение:** вертикальный подъём, отклонение по тангажу в сторону азимута пуска,
//!   гравитационный разворот (тяга вдоль скорости относительно воздуха, нулевой угол атаки),
//!   после ухода из плотной атмосферы — замкнутое наведение: на каждом шаге вертикальная
//!   составляющая тяги и составляющая поперёк плоскости выбираются так, чтобы к выключению высота
//!   стала равна высоте орбиты, а вертикальная скорость и выход из плоскости — нулю (как в
//!   терминальном наведении Apollo); остальная тяга — по горизонтали. Выключение двигателя — когда
//!   большая полуось достигла радиуса опорной орбиты. Угол отклонения подбирается так, чтобы масса
//!   на орбите была наибольшей.
//!
//! Не учитываются: ветер, дросселирование на максимальном скоростном напоре, ограничение
//! перегрузки, упругость и управляемость корпуса, притяжение Луны и Солнца (за 10 минут
//! подъёма — менее миллионной доли ускорения).

use anise::math::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};
use crate::spacecraft::G0;

/// WGS84.
pub const EARTH_EQUATORIAL_RADIUS_KM: f64 = 6378.137;
pub const EARTH_FLATTENING: f64 = 1.0 / 298.257_223_563;
/// Угловая скорость вращения Земли, рад/с.
pub const EARTH_ROTATION_RAD_S: f64 = 7.292_115e-5;
/// Сжатие Земли J2 (EGM2008).
pub const EARTH_J2: f64 = 1.082_626_68e-3;
/// Давление на уровне моря, Па.
const P0: f64 = 101_325.0;
/// Шаг интегрирования, с.
const DT: f64 = 0.2;

/// Ускорители или ступень.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct StageSpec {
    pub name: String,
    /// Число одинаковых блоков (у ускорителей), по умолчанию 1. Массы и тяга — на один блок.
    #[serde(default = "one")]
    pub count: u32,
    pub dry_mass_kg: f64,
    pub propellant_kg: f64,
    /// Тяга в пустоте и у земли, Н (на блок). Для профиля тяги — средние значения не нужны:
    /// задаются `burn_time_s` и `thrust_profile`, а наибольшая тяга вычисляется из массы топлива.
    #[serde(default)]
    pub thrust_vac_n: Option<f64>,
    #[serde(default)]
    pub thrust_sl_n: Option<f64>,
    pub isp_vac_s: f64,
    pub isp_sl_s: f64,
    /// Время работы, с (обязательно вместе с профилем тяги).
    #[serde(default)]
    pub burn_time_s: Option<f64>,
    /// Профиль тяги (твердотопливные): пары [доля времени работы, доля наибольшей тяги].
    #[serde(default)]
    pub thrust_profile: Option<Vec<[f64; 2]>>,
    /// Диаметр блока, м (для площади сопротивления).
    pub diameter_m: f64,
}

fn one() -> u32 {
    1
}

/// Данные ракеты для моделирования подъёма (секция `[preset.ascent]` в каталоге ракет).
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct AscentVehicle {
    /// Ускорители, работающие со старта параллельно с первой ступенью (если есть).
    #[serde(default)]
    pub boosters: Option<StageSpec>,
    /// Ступени по порядку. Первая включается на старте.
    pub stages: Vec<StageSpec>,
    /// Разгонная ступень ракеты (та, что разгоняет к Луне) может довыводить на орбиту.
    #[serde(default)]
    pub upper_stage_in_ascent: bool,
    /// Диаметр разгонной ступени, м.
    #[serde(default = "default_upper_diameter")]
    pub upper_stage_diameter_m: f64,
    /// Масса головного обтекателя (или другого сбрасываемого груза), кг, и высота сброса, км.
    #[serde(default)]
    pub fairing_mass_kg: f64,
    #[serde(default = "default_fairing_alt")]
    pub fairing_jettison_altitude_km: f64,
    /// Вертикальный подъём до начала разворота, с.
    #[serde(default = "default_vertical")]
    pub vertical_rise_s: f64,
    /// Пауза между выключением ступени и включением следующей, с.
    #[serde(default = "default_coast")]
    pub stage_coast_s: f64,
    /// Начальное приближение угла отклонения по тангажу после вертикального подъёма, градусы.
    /// Сам угол подбирается так, чтобы масса на орбите была наибольшей.
    #[serde(default = "default_kick")]
    pub pitch_kick_deg: f64,
    /// Предел скоростного напора, кПа: траектории с большим напором отбрасываются при подборе.
    #[serde(default = "default_max_q")]
    pub max_q_limit_kpa: f64,
    /// Источники данных.
    #[serde(default)]
    pub source: Option<String>,
}

fn default_upper_diameter() -> f64 {
    5.0
}
fn default_fairing_alt() -> f64 {
    110.0
}
fn default_vertical() -> f64 {
    10.0
}
fn default_coast() -> f64 {
    5.0
}
fn default_kick() -> f64 {
    3.0
}
fn default_max_q() -> f64 {
    // Типичный расчётный предел для ракет-носителей (у большинства фактический max-Q 30-40 кПа).
    45.0
}

impl StageSpec {
    fn validate(&self, what: &str) -> Vec<String> {
        let mut e = Vec::new();
        let pos = |v: f64| v.is_finite() && v > 0.0;
        if self.count == 0 || self.count > 16 {
            e.push(format!("{what} '{}': count must be in [1, 16]", self.name));
        }
        if !(self.dry_mass_kg >= 0.0) || !pos(self.propellant_kg) {
            e.push(format!("{what} '{}': masses must be positive", self.name));
        }
        if !pos(self.isp_vac_s) || !pos(self.isp_sl_s) || self.isp_sl_s > self.isp_vac_s {
            e.push(format!("{what} '{}': need 0 < isp_sl_s <= isp_vac_s", self.name));
        }
        if !pos(self.diameter_m) {
            e.push(format!("{what} '{}': diameter_m must be positive", self.name));
        }
        match (&self.thrust_profile, self.burn_time_s, self.thrust_vac_n) {
            (Some(p), Some(t), _) => {
                if p.len() < 2 || !pos(t) || p.iter().any(|q| !(q[0] >= 0.0 && q[0] <= 1.0 && q[1] >= 0.0)) {
                    e.push(format!("{what} '{}': bad thrust_profile or burn_time_s", self.name));
                }
            }
            (Some(_), None, _) => e.push(format!("{what} '{}': thrust_profile needs burn_time_s", self.name)),
            (None, _, Some(f)) if pos(f) => {}
            _ => e.push(format!("{what} '{}': thrust_vac_n (or thrust_profile + burn_time_s) is required", self.name)),
        }
        e
    }

    /// Площадь сопла (умноженная на давление у земли даёт разницу тяги в пустоте и у земли), Н/Па.
    fn exit_area(&self, f_vac: f64) -> f64 {
        let f_sl = self.thrust_sl_n.unwrap_or(f_vac * self.isp_sl_s / self.isp_vac_s);
        ((f_vac - f_sl) / P0).max(0.0)
    }

    /// Интеграл профиля тяги по доле времени.
    fn profile_integral(p: &[[f64; 2]]) -> f64 {
        p.windows(2).map(|w| 0.5 * (w[0][1] + w[1][1]) * (w[1][0] - w[0][0])).sum()
    }

    /// Наибольшая тяга в пустоте на блок, Н.
    pub fn peak_thrust_vac(&self) -> f64 {
        match (&self.thrust_profile, self.burn_time_s) {
            (Some(p), Some(t)) => self.propellant_kg * G0 * self.isp_vac_s / (t * Self::profile_integral(p)),
            _ => self.thrust_vac_n.unwrap_or(0.0),
        }
    }

    /// Тяга в пустоте на блок в момент `t` после включения, Н.
    fn thrust_vac_at(&self, t: f64) -> f64 {
        let peak = self.peak_thrust_vac();
        match (&self.thrust_profile, self.burn_time_s) {
            (Some(p), Some(tb)) => {
                let x = t / tb;
                let k = p
                    .windows(2)
                    .find(|w| x >= w[0][0] && x <= w[1][0])
                    .map(|w| w[0][1] + (w[1][1] - w[0][1]) * (x - w[0][0]) / (w[1][0] - w[0][0]).max(1e-12))
                    .unwrap_or(0.0);
                peak * k
            }
            _ => peak,
        }
    }

    fn area_m2(&self) -> f64 {
        self.count as f64 * std::f64::consts::PI * self.diameter_m * self.diameter_m / 4.0
    }
}

impl AscentVehicle {
    pub fn validate(&self) -> Vec<String> {
        let mut e = Vec::new();
        if let Some(b) = &self.boosters {
            e.extend(b.validate("booster"));
        }
        if self.stages.is_empty() {
            e.push("ascent: at least one stage is required".into());
        }
        for s in &self.stages {
            e.extend(s.validate("stage"));
        }
        if !(self.fairing_mass_kg >= 0.0) || !(self.fairing_jettison_altitude_km > 0.0) {
            e.push("ascent: bad fairing parameters".into());
        }
        if !(self.vertical_rise_s >= 0.0 && self.vertical_rise_s <= 60.0) {
            e.push("ascent: vertical_rise_s must be in [0, 60]".into());
        }
        if !(self.stage_coast_s >= 0.0 && self.stage_coast_s <= 60.0) {
            e.push("ascent: stage_coast_s must be in [0, 60]".into());
        }
        if !(self.pitch_kick_deg > 0.0 && self.pitch_kick_deg < 30.0) {
            e.push("ascent: pitch_kick_deg must be in (0, 30)".into());
        }
        e
    }
}

/// Разгонная ступень (из основных данных ракеты).
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct UpperStage {
    pub dry_mass_kg: f64,
    pub propellant_kg: f64,
    pub thrust_n: f64,
    pub isp_s: f64,
}

/// Атмосфера: плотность (кг/м³), давление (Па), скорость звука (м/с) на геометрической высоте `h_km`.
pub fn atmosphere(h_km: f64) -> (f64, f64, f64) {
    // US Standard Atmosphere 1976: слои по геопотенциальной высоте.
    const R: f64 = 287.053;
    const GAMMA: f64 = 1.4;
    const LAYERS: [(f64, f64, f64, f64); 8] = [
        // (база, км; температура, K; давление, Па; градиент, K/км)
        (0.0, 288.15, 101_325.0, -6.5),
        (11.0, 216.65, 22_632.06, 0.0),
        (20.0, 216.65, 5_474.889, 1.0),
        (32.0, 228.65, 868.0187, 2.8),
        (47.0, 270.65, 110.9063, 0.0),
        (51.0, 270.65, 66.93887, -2.8),
        (71.0, 214.65, 3.956_420, -2.0),
        (84.852, 186.946, 0.373_384, 0.0),
    ];
    let h = h_km.max(-1.0);
    if h < 86.0 {
        let hg = 6356.766 * h / (6356.766 + h);
        let l = LAYERS.iter().rev().find(|l| hg >= l.0).unwrap_or(&LAYERS[0]);
        let (hb, tb, pb, lr) = *l;
        let t = tb + lr * (hg - hb);
        let g0r = 9.80665 / R * 1000.0; // на км
        let p = if lr == 0.0 { pb * (-g0r * (hg - hb) / tb).exp() } else { pb * (tb / t).powf(g0r / lr) };
        return (p / (R * t), p, (GAMMA * R * t).sqrt());
    }
    // Выше 86 км: плотность по таблице стандартной атмосферы, логарифмическая интерполяция.
    const HI: [(f64, f64); 12] = [
        (86.0, 6.958e-6),
        (100.0, 5.604e-7),
        (110.0, 9.708e-8),
        (120.0, 2.222e-8),
        (130.0, 8.152e-9),
        (150.0, 2.076e-9),
        (180.0, 5.194e-10),
        (200.0, 2.541e-10),
        (250.0, 6.073e-11),
        (300.0, 1.916e-11),
        (500.0, 5.215e-13),
        (1000.0, 3.561e-15),
    ];
    if h >= 1000.0 {
        return (0.0, 0.0, 300.0);
    }
    let i = HI.iter().rposition(|x| h >= x.0).unwrap_or(0).min(HI.len() - 2);
    let (h0, r0) = HI[i];
    let (h1, r1) = HI[i + 1];
    let rho = r0 * ((r1 / r0).ln() * (h - h0) / (h1 - h0)).exp();
    (rho, rho * R * 186.87, (GAMMA * R * 186.87).sqrt())
}

/// Коэффициент лобового сопротивления по числу Маха (типичный для ракет-носителей).
pub fn drag_coefficient(mach: f64) -> f64 {
    const T: [(f64, f64); 9] =
        [(0.0, 0.30), (0.6, 0.30), (0.8, 0.35), (1.0, 0.55), (1.2, 0.60), (1.5, 0.50), (2.0, 0.42), (3.0, 0.33), (5.0, 0.28)];
    if mach >= 5.0 {
        return 0.28;
    }
    let i = T.iter().rposition(|x| mach >= x.0).unwrap_or(0).min(T.len() - 2);
    let (m0, c0) = T[i];
    let (m1, c1) = T[i + 1];
    c0 + (c1 - c0) * (mach - m0) / (m1 - m0)
}

/// Событие подъёма.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct AscentEvent {
    /// liftoff, pitch_kick, max_q, booster_separation, stage_cutoff, stage_ignition,
    /// fairing_jettison, upper_stage_ignition, insertion.
    pub name: &'static str,
    pub detail: String,
    pub t_s: f64,
    pub altitude_km: f64,
    /// Скорость относительно Земли (воздуха), м/с.
    pub ground_speed_m_s: f64,
    /// Скорость в инерциальной системе, м/с.
    pub inertial_speed_m_s: f64,
    pub mass_kg: f64,
    /// Положение и скорость в местной системе моделирования.
    #[serde(skip)]
    pub r: Vector3,
    #[serde(skip)]
    pub v: Vector3,
}

/// Точка подъёма (в местной инерциальной системе, см. [`AscentResult`]).
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct AscentPoint {
    pub t_s: f64,
    #[serde(skip)]
    pub r: Vector3,
    #[serde(skip)]
    pub v: Vector3,
    pub altitude_km: f64,
    pub mass_kg: f64,
    pub thrust_n: f64,
    pub dynamic_pressure_pa: f64,
    pub mach: f64,
    /// Тангаж — угол оси над местным горизонтом, градусы.
    pub pitch_deg: f64,
}

/// Результат моделирования подъёма.
///
/// Местная инерциальная система: начало в центре Земли, ось z — ось вращения Земли, космодром
/// в момент старта лежит в плоскости xz (x > 0). Перевод в J2000: см. [`AscentResult::to_j2000`].
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct AscentResult {
    pub converged: bool,
    pub iterations: u32,
    /// Подобранный угол отклонения по тангажу (наибольшая масса на орбите), градусы.
    pub pitch_kick_deg: f64,
    /// От старта до выключения двигателя на орбите, с.
    pub duration_s: f64,
    /// Угол в плоскости орбиты от космодрома (в момент старта) до точки выхода, рад.
    pub downrange_angle_rad: f64,
    /// Дальность по поверхности, км.
    pub downrange_km: f64,
    pub insertion_altitude_km: f64,
    pub insertion_speed_km_s: f64,
    /// Перигей и апогей орбиты после выключения над экватором, км.
    pub perigee_km: f64,
    pub apogee_km: f64,
    /// Отклонение плоскости от заданной, градусы.
    pub plane_error_deg: f64,
    pub max_q_pa: f64,
    pub max_q_t_s: f64,
    /// Наибольшая перегрузка (тяга + сопротивление), g.
    pub max_acceleration_g: f64,
    /// Наибольший скоростной напор превышает предел ракеты (ни одна траектория в него не уложилась).
    pub max_q_above_limit: bool,
    /// Потери скорости: гравитационные, на сопротивление, на управление, м/с.
    pub gravity_loss_m_s: f64,
    pub drag_loss_m_s: f64,
    pub steering_loss_m_s: f64,
    /// Потеря тяги от давления воздуха на срез сопла, м/с.
    pub pressure_loss_m_s: f64,
    /// Идеальный разгон по тяге в пустоте (∫F_vac/m dt), м/с.
    pub ideal_delta_v_m_s: f64,
    /// Скорость на старте (вращение Земли), м/с.
    pub initial_speed_m_s: f64,
    /// Разгонная ступень включалась на подъёме.
    pub upper_stage_used: bool,
    /// Остаток топлива разгонной ступени после выхода на орбиту, кг.
    pub upper_stage_propellant_left_kg: f64,
    /// Δv, который может дать разгонная ступень после выхода на орбиту, м/с.
    pub upper_stage_delta_v_left_m_s: f64,
    /// Масса на орбите (разгонная ступень + груз), кг.
    pub mass_in_orbit_kg: f64,
    pub liftoff_mass_kg: f64,
    pub liftoff_thrust_to_weight: f64,
    pub events: Vec<AscentEvent>,
    /// Точки каждую секунду (для журнала траектории; в журнал информации не пишутся).
    #[serde(skip)]
    pub points: Vec<AscentPoint>,
    /// Состояние на выходе в местной системе.
    #[serde(skip)]
    pub insertion_r: Vector3,
    #[serde(skip)]
    pub insertion_v: Vector3,
    #[serde(skip)]
    pub site_r: Vector3,
    #[serde(skip)]
    pub plane_normal: Vector3,
}

/// Входные данные моделирования.
#[derive(Debug, Clone)]
pub struct AscentInput<'a> {
    pub vehicle: &'a AscentVehicle,
    pub upper: UpperStage,
    /// Груз над разгонной ступенью (лунный аппарат), кг.
    pub payload_kg: f64,
    /// Геодезическая широта и высота космодрома.
    pub lat_deg: f64,
    pub alt_km: f64,
    /// Инерциальный азимут плоскости опорной орбиты в точке старта, градусы.
    pub azimuth_deg: f64,
    /// Радиус круговой опорной орбиты, км.
    pub target_radius_km: f64,
    pub mu_km3_s2: f64,
}

/// Геоцентрическая широта и радиус эллипсоида для геодезической широты.
fn ellipsoid_point(lat_deg: f64, alt_km: f64) -> (f64, f64) {
    let lat = lat_deg.to_radians();
    let e2 = EARTH_FLATTENING * (2.0 - EARTH_FLATTENING);
    let n = EARTH_EQUATORIAL_RADIUS_KM / (1.0 - e2 * lat.sin().powi(2)).sqrt();
    let x = (n + alt_km) * lat.cos();
    let z = (n * (1.0 - e2) + alt_km) * lat.sin();
    (z.atan2(x), x.hypot(z))
}

/// Высота над эллипсоидом (приближённо, по геоцентрической широте), км.
fn altitude_km(r: &Vector3) -> f64 {
    let rn = r.norm();
    let s = r.z / rn;
    let b = EARTH_EQUATORIAL_RADIUS_KM * (1.0 - EARTH_FLATTENING);
    let a = EARTH_EQUATORIAL_RADIUS_KM;
    let c = (1.0 - s * s).max(0.0).sqrt();
    rn - a * b / ((b * c).powi(2) + (a * s).powi(2)).sqrt()
}

#[derive(Clone, Copy, PartialEq)]
enum Burning {
    /// Работает ступень с номером (`stages.len()` — разгонная ступень).
    Stage(usize),
    /// Пауза перед включением ступени.
    Coast(usize, f64),
}

struct Sim<'a> {
    inp: &'a AscentInput<'a>,
    kick: f64,
    record: bool,
}

struct Outcome {
    /// Вышли на орбиту (большая полуось достигла заданной).
    inserted: bool,
    t: f64,
    r: Vector3,
    v: Vector3,
    upper_used: bool,
    upper_left: f64,
    mass: f64,
    events: Vec<AscentEvent>,
    points: Vec<AscentPoint>,
    max_q: (f64, f64),
    max_acc: f64,
    losses: (f64, f64, f64),
    /// Потери от давления воздуха на тягу и идеальный разгон (по тяге в пустоте), м/с; начальная скорость, км/с.
    p_loss: f64,
    ideal: f64,
    v_start: f64,
    /// tan угла траектории в начале наведения.
    #[allow(dead_code)]
    tan_g_switch: Option<f64>,
    /// Вся масса в момент выключения (с остатком топлива работающей ступени), кг.
    total_mass: f64,
}

impl Sim<'_> {
    /// Полёт с заданными параметрами наведения.
    fn run(&self) -> Outcome {
        let inp = self.inp;
        let veh = inp.vehicle;
        let mu = inp.mu_km3_s2;
        let omega = Vector3::new(0.0, 0.0, EARTH_ROTATION_RAD_S);
        let (lat_c, r_pad) = ellipsoid_point(inp.lat_deg, inp.alt_km);
        let up0 = Vector3::new(lat_c.cos(), 0.0, lat_c.sin());
        let east0 = Vector3::new(0.0, 1.0, 0.0);
        let north0 = up0.cross(&east0);
        let az = inp.azimuth_deg.to_radians();
        let d_az = north0 * az.cos() + east0 * az.sin();
        let h = up0.cross(&d_az).normalize();
        let e_target = -mu / (2.0 * inp.target_radius_km);

        let mut r = up0 * r_pad;
        let mut v = omega.cross(&r);
        let ns = veh.stages.len();
        let b_count = veh.boosters.as_ref().map_or(0.0, |b| b.count as f64);
        let mut booster_prop = veh.boosters.as_ref().map_or(0.0, |b| b.propellant_kg * b.count as f64);
        let mut boosters_on = veh.boosters.is_some();
        let mut stage_prop: Vec<f64> = veh.stages.iter().map(|s| s.propellant_kg * s.count as f64).collect();
        stage_prop.push(inp.upper.propellant_kg);
        let mut stage_t0 = vec![0.0f64; ns + 1];
        let mut fairing_on = veh.fairing_mass_kg > 0.0;
        let mut burning = Burning::Stage(0);
        let mut t = 0.0;
        let mut events = Vec::new();
        let mut points = Vec::new();
        let (mut max_q, mut max_q_t, mut max_acc) = (0.0f64, 0.0f64, 0.0f64);
        let (mut g_loss, mut d_loss, mut s_loss) = (0.0, 0.0, 0.0);
        let (mut p_loss, mut ideal) = (0.0, 0.0);
        let v_start = v.norm();
        let mut upper_used = false;
        let mut guidance: Option<(f64, f64, f64)> = None; // (t_start, tan γ_start, T_est)
        let mut kick_done = false;
        let mut max_q_logged = false;

        let mass = |booster_prop: f64, boosters_on: bool, stage_prop: &[f64], first_stage: usize, fairing_on: bool| {
            let mut m = inp.payload_kg + inp.upper.dry_mass_kg + stage_prop[ns];
            for (i, s) in veh.stages.iter().enumerate().skip(first_stage) {
                m += s.dry_mass_kg * s.count as f64 + stage_prop[i];
            }
            if boosters_on {
                let b = veh.boosters.as_ref().expect("boosters");
                m += b.dry_mass_kg * b.count as f64 + booster_prop;
            }
            if fairing_on {
                m += veh.fairing_mass_kg;
            }
            m
        };
        let mut first_stage = 0usize;
        let ev = |name: &'static str, detail: String, t: f64, r: &Vector3, v: &Vector3, m: f64| AscentEvent {
            name,
            detail,
            t_s: (t * 10.0).round() / 10.0,
            altitude_km: altitude_km(r),
            ground_speed_m_s: (v - omega.cross(r)).norm() * 1000.0,
            inertial_speed_m_s: v.norm() * 1000.0,
            mass_kg: m,
            r: *r,
            v: *v,
        };
        let m0 = mass(booster_prop, boosters_on, &stage_prop, 0, fairing_on);
        if self.record {
            events.push(ev("liftoff", "liftoff".into(), 0.0, &r, &v, m0));
        }

        // Тяга (Н), расход (кг/с), площадь сопротивления (м²) в момент t на высоте h.
        let propulsion = |t: f64, p_amb: f64, burning: Burning, boosters_on: bool, stage_t0: &[f64], first_stage: usize| {
            let (mut f, mut mdot, mut mdot_b) = (0.0, 0.0, 0.0);
            let mut area = 0.0f64;
            if boosters_on {
                let b = veh.boosters.as_ref().expect("boosters");
                let fv = b.thrust_vac_at(t);
                f += b.count as f64 * (fv - b.exit_area(b.peak_thrust_vac()) * p_amb).max(0.0);
                mdot_b = b.count as f64 * fv / (G0 * b.isp_vac_s);
                area += b.area_m2();
            }
            if let Burning::Stage(i) = burning {
                if i < ns {
                    let s = &veh.stages[i];
                    let fv = s.thrust_vac_at(t - stage_t0[i]);
                    f += s.count as f64 * (fv - s.exit_area(s.peak_thrust_vac()) * p_amb).max(0.0);
                    mdot += s.count as f64 * fv / (G0 * s.isp_vac_s);
                } else {
                    f += inp.upper.thrust_n;
                    mdot += inp.upper.thrust_n / (G0 * inp.upper.isp_s);
                }
            }
            let core = veh.stages.get(first_stage).map_or(veh.upper_stage_diameter_m, |s| s.diameter_m);
            area = area.max(0.0) + std::f64::consts::PI * core * core / 4.0;
            (f, mdot, mdot_b, area)
        };

        loop {
            if t > 3000.0 {
                break;
            }
            let p_amb = atmosphere(altitude_km(&r)).1;
            let (thrust, mdot, mdot_b, _) = propulsion(t, p_amb, burning, boosters_on, &stage_t0, first_stage);
            let guidance_now = guidance;
            let stage_prop_now = stage_prop.clone();
            let m_now = mass(booster_prop, boosters_on, &stage_prop, first_stage, fairing_on);
            // Время работы двигателей, нужное на прирост скорости dv (м/с), по оставшимся ступеням.
            let t_go_for = |dv: f64, m: f64, first: usize, prop: &[f64]| -> f64 {
                let (mut t_go, mut m_est, mut left) = (0.0, m, dv);
                let current = match burning {
                    Burning::Stage(i) | Burning::Coast(i, _) => i,
                };
                for i in current.max(first)..=ns {
                    if left <= 0.0 {
                        break;
                    }
                    let (fv, isp, dry) = if i < ns {
                        let st = &veh.stages[i];
                        (st.peak_thrust_vac() * st.count as f64, st.isp_vac_s, st.dry_mass_kg * st.count as f64)
                    } else if veh.upper_stage_in_ascent {
                        (inp.upper.thrust_n, inp.upper.isp_s, 0.0)
                    } else {
                        break;
                    };
                    let ve = G0 * isp;
                    let dv_stage = ve * (m_est / (m_est - prop[i]).max(1.0)).ln();
                    let use_dv = left.min(dv_stage);
                    t_go += m_est * (1.0 - (-use_dv / ve).exp()) / (fv / ve);
                    m_est -= prop[i] + dry;
                    left -= use_dv;
                }
                t_go
            };
            let m = mass(booster_prop, boosters_on, &stage_prop, first_stage, fairing_on);

            // Направление тяги.
            let dir = |r: &Vector3, v: &Vector3, t: f64| -> Vector3 {
                let up = r.normalize();
                let downrange = h.cross(&up).normalize();
                let v_air = v - omega.cross(r);
                if t < veh.vertical_rise_s {
                    return up;
                }
                if guidance_now.is_some() {
                    // Замкнутое наведение: вертикальное ускорение и ускорение поперёк плоскости
                    // выбираются так, чтобы к выключению (через t_go) высота стала равна высоте орбиты,
                    // а вертикальная скорость и выход из плоскости — нулю; остальная тяга — по горизонтали.
                    let rn = r.norm();
                    let v_r = v.dot(&up);
                    let v_h = v.dot(&downrange);
                    let y = r.dot(&h);
                    let yd = v.dot(&h);
                    let v_t = (mu / inp.target_radius_km).sqrt();
                    let dv_go = ((v_t - v_h).powi(2) + v_r * v_r + yd * yd).sqrt() * 1000.0;
                    let t_go = t_go_for(dv_go, m_now, first_stage, &stage_prop_now).max(3.0);
                    let a_t = (thrust / m_now / 1000.0).max(1e-9);
                    let g_eff = mu / (rn * rn) - v_h * v_h / rn;
                    let a_r = (g_eff + 6.0 * (inp.target_radius_km - rn) / (t_go * t_go) - 4.0 * v_r / t_go).clamp(-0.95 * a_t, 0.95 * a_t);
                    let a_y = (-(6.0 * y / (t_go * t_go) + 4.0 * yd / t_go)).clamp(-0.3 * a_t, 0.3 * a_t);
                    let a_h = (a_t * a_t - a_r * a_r - a_y * a_y).max(0.0).sqrt();
                    return (up * a_r + downrange * a_h + h * a_y) / a_t;
                }
                // Отклонение и гравитационный разворот в плоскости орбиты: тангаж — по проекции
                // скорости относительно воздуха на эту плоскость (нулевой угол атаки по тангажу),
                // рысканием ракета удерживается в плоскости.
                let tilt_cmd = self.kick * ((t - veh.vertical_rise_s) / 10.0).min(1.0);
                let (va_up, va_down) = (v_air.dot(&up), v_air.dot(&downrange));
                let tilt_v = va_down.atan2(va_up);
                let tilt = if t > veh.vertical_rise_s + 10.0 && tilt_v >= tilt_cmd { tilt_v } else { tilt_cmd };
                up * tilt.cos() + downrange * tilt.sin()
            };

            // Производные: ускорение (км/с²) с тягой, сопротивлением, притяжением с J2.
            let accel = |r: &Vector3, v: &Vector3, m: f64, t: f64| -> (Vector3, Vector3, f64) {
                let rn = r.norm();
                let z = r.z;
                let k = 1.5 * EARTH_J2 * mu * EARTH_EQUATORIAL_RADIUS_KM.powi(2) / rn.powi(5);
                let grav = -r * (mu / rn.powi(3))
                    + Vector3::new(r.x * (5.0 * z * z / (rn * rn) - 1.0), r.y * (5.0 * z * z / (rn * rn) - 1.0), z * (5.0 * z * z / (rn * rn) - 3.0)) * k;
                let alt = altitude_km(r);
                let (rho, p_amb, a_snd) = atmosphere(alt);
                let v_air = v - omega.cross(r);
                let va = v_air.norm() * 1000.0;
                let (f, _, _, ar) = propulsion(t, p_amb, burning, boosters_on, &stage_t0, first_stage);
                let cd = drag_coefficient(va / a_snd);
                let drag_n = 0.5 * rho * va * va * cd * ar;
                let drag = if va > 1e-6 { -v_air.normalize() * (drag_n / m / 1000.0) } else { Vector3::zeros() };
                let th = dir(r, v, t) * (f / m / 1000.0);
                (grav + drag + th, drag, f)
            };

            // Шаг RK4 с массой, меняющейся линейно.
            let step = |r: &Vector3, v: &Vector3, dt: f64| -> (Vector3, Vector3) {
                let total_mdot = mdot + mdot_b;
                let (a1, _, _) = accel(r, v, m, t);
                let (r2, v2) = (r + v * (dt / 2.0), v + a1 * (dt / 2.0));
                let (a2, _, _) = accel(&r2, &v2, m - total_mdot * dt / 2.0, t + dt / 2.0);
                let (r3, v3) = (r + v2 * (dt / 2.0), v + a2 * (dt / 2.0));
                let (a3, _, _) = accel(&r3, &v3, m - total_mdot * dt / 2.0, t + dt / 2.0);
                let (r4, v4) = (r + v3 * dt, v + a3 * dt);
                let (a4, _, _) = accel(&r4, &v4, m - total_mdot * dt, t + dt);
                (r + (v + v2 * 2.0 + v3 * 2.0 + v4) * (dt / 6.0), v + (a1 + a2 * 2.0 + a3 * 2.0 + a4) * (dt / 6.0))
            };

            // Диагностика: напор, перегрузка, потери.
            let alt = altitude_km(&r);
            let (rho, _, a_snd) = atmosphere(alt);
            let v_air = v - omega.cross(&r);
            let va = v_air.norm() * 1000.0;
            let q = 0.5 * rho * va * va;
            let (_, drag_vec, f_now) = accel(&r, &v, m, t);
            let thrust_dir = dir(&r, &v, t);
            if q > max_q {
                max_q = q;
                max_q_t = t;
            }
            let sensed = (thrust_dir * (f_now / m / 1000.0) + drag_vec).norm() * 1000.0 / G0;
            max_acc = max_acc.max(sensed);
            let up = r.normalize();
            let vn = v.norm().max(1e-9);
            let sin_gamma = v.dot(&up) / vn;
            g_loss += mu / r.norm_squared() * sin_gamma * DT * 1000.0;
            d_loss += drag_vec.norm() * DT * 1000.0;
            let cos_a = thrust_dir.dot(&(v / vn)).clamp(-1.0, 1.0);
            s_loss += f_now / m * (1.0 - cos_a) * DT;
            let f_vac = propulsion(t, 0.0, burning, boosters_on, &stage_t0, first_stage).0;
            p_loss += (f_vac - f_now) / m * DT;
            ideal += f_vac / m * DT;
            if self.record {
                if !max_q_logged && t > 5.0 && q < 0.97 * max_q && max_q > 1000.0 {
                    max_q_logged = true;
                    events.push(AscentEvent {
                        name: "max_q",
                        detail: format!("maximum dynamic pressure {:.1} kPa", max_q / 1000.0),
                        ..ev("max_q", String::new(), max_q_t, &r, &v, m)
                    });
                }
                if ((t / DT).round() as i64) % 5 == 0 {
                    let horiz = (v - up * v.dot(&up)).norm();
                    points.push(AscentPoint {
                        t_s: t,
                        r,
                        v,
                        altitude_km: alt,
                        mass_kg: m,
                        thrust_n: f_now,
                        dynamic_pressure_pa: q,
                        mach: va / a_snd,
                        pitch_deg: thrust_dir.dot(&up).atan2(thrust_dir.dot(&(v - up * v.dot(&up))) / horiz.max(1e-9)).to_degrees(),
                    });
                }
            }

            // Вход в режим наведения в пустоте: ускорители отработали, напор упал.
            if guidance.is_none() && !boosters_on && t > veh.vertical_rise_s + 10.0 && alt > 35.0 && q < 5000.0 {
                let up = r.normalize();
                let downrange = h.cross(&up).normalize();
                let tan_g = (v.dot(&up) / v.dot(&downrange).max(1e-6)).clamp(-1.0, 20.0);
                guidance = Some((t, tan_g, 0.0));
                if self.record {
                    events.push(ev("guidance", "closed-loop guidance starts (above the dense atmosphere)".into(), t, &r, &v, m));
                }
            }
            if !kick_done && t >= veh.vertical_rise_s {
                kick_done = true;
                if self.record {
                    events.push(ev("pitch_kick", format!("pitch-over {:.2} deg towards azimuth {:.1} deg", self.kick.to_degrees(), inp.azimuth_deg), t, &r, &v, m));
                }
            }

            // Шаг.
            let (mut r_new, mut v_new) = step(&r, &v, DT);
            let mut dt = DT;
            // Выключение на орбите: большая полуось достигла заданной (только когда двигатель работает).
            let energy = |r: &Vector3, v: &Vector3| v.norm_squared() / 2.0 - mu / r.norm();
            let thrusting = matches!(burning, Burning::Stage(_)) && thrust > 0.0;
            if thrusting && energy(&r_new, &v_new) >= e_target && guidance_now.is_some() {
                let (mut lo, mut hi) = (0.0, DT);
                for _ in 0..40 {
                    let mid = 0.5 * (lo + hi);
                    let (rm, vm) = step(&r, &v, mid);
                    if energy(&rm, &vm) >= e_target { hi = mid } else { lo = mid }
                }
                (r_new, v_new) = step(&r, &v, hi);
                dt = hi;
                let burned = (mdot + mdot_b) * dt;
                if let Burning::Stage(i) = burning {
                    stage_prop[i] = (stage_prop[i] - mdot * dt).max(0.0);
                    if i == ns {
                        upper_used = true;
                    }
                }
                let _ = burned;
                t += dt;
                let m_end = mass(booster_prop, boosters_on, &stage_prop, first_stage, fairing_on);
                if self.record {
                    events.push(ev("insertion", "engine cutoff: parking orbit reached".into(), t, &r_new, &v_new, m_end));
                    let up = r_new.normalize();
                    let horiz = (v_new - up * v_new.dot(&up)).norm();
                    points.push(AscentPoint {
                        t_s: t,
                        r: r_new,
                        v: v_new,
                        altitude_km: altitude_km(&r_new),
                        mass_kg: m_end,
                        thrust_n: 0.0,
                        dynamic_pressure_pa: 0.0,
                        mach: 0.0,
                        pitch_deg: v_new.dot(&up).atan2(horiz).to_degrees(),
                    });
                }
                // Масса на орбите: всё, что выше ещё работающей ступени.
                let in_orbit = inp.payload_kg + inp.upper.dry_mass_kg + stage_prop[ns];
                let total_at_cutoff = m_end;
                return Outcome {
                    inserted: true,
                    t,
                    r: r_new,
                    v: v_new,
                    upper_used,
                    upper_left: stage_prop[ns],
                    mass: in_orbit,
                    events,
                    points,
                    max_q: (max_q, max_q_t),
                    max_acc,
                    losses: (g_loss, d_loss, s_loss),
                    p_loss,
                    ideal,
                    v_start,
                    tan_g_switch: guidance_now.map(|g| g.1),
                    total_mass: total_at_cutoff,
                };
            }

            // Расход топлива и смена ступеней.
            if boosters_on {
                booster_prop -= mdot_b * dt;
            }
            match burning {
                Burning::Stage(i) => {
                    stage_prop[i] -= mdot * dt;
                    if i == ns {
                        upper_used = true;
                    }
                }
                Burning::Coast(i, until) => {
                    if t + dt >= until {
                        burning = Burning::Stage(i);
                        stage_t0[i] = t + dt;
                        if self.record {
                            let name = if i == ns { "upper_stage_ignition" } else { "stage_ignition" };
                            let label = if i == ns { "upper stage".to_string() } else { veh.stages[i].name.clone() };
                            events.push(ev(name, format!("{label} ignition"), t + dt, &r_new, &v_new, m));
                        }
                    }
                }
            }
            r = r_new;
            v = v_new;
            t += dt;
            if altitude_km(&r) < -0.1 && t > 5.0 {
                break;
            }
            if boosters_on && booster_prop <= 0.0 {
                boosters_on = false;
                booster_prop = 0.0;
                if self.record {
                    let b = veh.boosters.as_ref().expect("boosters");
                    let m_after = mass(booster_prop, boosters_on, &stage_prop, first_stage, fairing_on);
                    events.push(ev("booster_separation", format!("{} x {} separation", b_count as u32, b.name), t, &r, &v, m_after));
                }
            }
            if fairing_on && altitude_km(&r) >= veh.fairing_jettison_altitude_km {
                fairing_on = false;
                if self.record {
                    let m_after = mass(booster_prop, boosters_on, &stage_prop, first_stage, fairing_on);
                    events.push(ev("fairing_jettison", "payload fairing jettison".into(), t, &r, &v, m_after));
                }
            }
            if let Burning::Stage(i) = burning {
                if stage_prop[i] <= 0.0 {
                    stage_prop[i] = 0.0;
                    if i == ns {
                        break; // разгонная ступень исчерпана, орбита не достигнута
                    }
                    first_stage = i + 1;
                    let next = i + 1;
                    if self.record {
                        let m_after = mass(booster_prop, boosters_on, &stage_prop, first_stage, fairing_on);
                        events.push(ev("stage_cutoff", format!("{} cutoff and separation", veh.stages[i].name), t, &r, &v, m_after));
                    }
                    if next == ns && !veh.upper_stage_in_ascent {
                        break; // нижние ступени исчерпаны, а разгонная на подъёме не работает
                    }
                    burning = Burning::Coast(next, t + veh.stage_coast_s);
                }
            }
        }
        let m_end = mass(booster_prop, boosters_on, &stage_prop, first_stage.min(ns), fairing_on);
        Outcome {
            inserted: false,
            t,
            r,
            v,
            upper_used,
            upper_left: stage_prop[ns],
            mass: m_end,
            events,
            points,
            max_q: (max_q, max_q_t),
            max_acc,
            losses: (g_loss, d_loss, s_loss),
            p_loss,
            ideal,
            v_start,
            tan_g_switch: guidance.map(|g| g.1),
            total_mass: m_end,
        }
    }
}

/// Невязки: высота выключения (км) и вертикальная скорость (м/с, деленная на 10).
fn residuals(o: &Outcome, target: f64) -> Option<[f64; 2]> {
    if !o.inserted {
        return None;
    }
    let up = o.r.normalize();
    Some([o.r.norm() - target, o.v.dot(&up) * 100.0])
}

/// Исходный код модели подъёма — для ключа кэша (изменение модели делает кэш недействительным).
pub const MODEL_SOURCE: &str = include_str!("ascent.rs");

pub const REASON_ASCENT_ORBIT: &str =
    "the rocket cannot reach the parking orbit with this payload (ascent simulation runs out of propellant)";

/// Моделирует подъём и подбирает наведение. `Ok(Err(причина))` — орбита недостижима.
pub fn simulate(inp: &AscentInput) -> Result<std::result::Result<AscentResult, &'static str>> {
    let errors = inp.vehicle.validate();
    if !errors.is_empty() {
        return Err(Error::Config(format!("rocket ascent data: {}", errors.join("; "))));
    }
    let target = inp.target_radius_km;
    // Качество выведения: вся масса в момент выключения (чем больше осталось топлива, тем
    // экономичнее траектория), если орбита точная; с пределом скоростного напора или без.
    let run_score = |kick: f64, q_cap: f64| -> Option<f64> {
        let o = Sim { inp, kick, record: false }.run();
        let f = residuals(&o, target)?;
        let q_ok = o.max_q.0 <= q_cap;
        (f[0].abs() < 1.0 && f[1].abs() < 2.0 && q_ok).then_some(o.total_mass)
    };
    // Перебор угла отклонения по сетке (каждый вариант считается один раз), затем уточнение
    // золотым сечением вокруг лучшего. Если ни одна траектория не укладывается в предел напора —
    // допускается напор на 5% выше наименьшего достижимого (предупреждение в результате).
    let nominal = inp.vehicle.pitch_kick_deg;
    // Сетка с шагом 1° (золотое сечение затем уточняет в пределах ±1° до тысячных градуса).
    let mut grid: Vec<f64> = (1..=25).map(|i| i as f64).collect();
    grid.push(nominal);
    let mut iterations = 0u32;
    let q_max = inp.vehicle.max_q_limit_kpa * 1000.0;
    // Варианты сетки независимы — считаются параллельно (порядок результатов сохраняется).
    use rayon::prelude::*;
    iterations += grid.len() as u32;
    let tried: Vec<(f64, f64, f64)> = grid
        .par_iter()
        .filter_map(|&kd| {
            let o = Sim { inp, kick: kd.to_radians(), record: false }.run();
            let f = residuals(&o, target)?;
            (f[0].abs() < 1.0 && f[1].abs() < 2.0).then_some((kd, o.total_mass, o.max_q.0))
        })
        .collect();
    // Предел недостижим — он поднимается до наименьшего достижимого напора + 5%.
    let q_min = tried.iter().map(|t| t.2).fold(f64::INFINITY, f64::min);
    let q_cap = if q_min <= q_max { q_max } else { q_min * 1.05 };
    let best = tried.iter().filter(|t| t.2 <= q_cap).max_by(|a, b| a.1.total_cmp(&b.1)).map(|t| (t.0, t.1));
    let Some((k0, _)) = best else { return Ok(Err(REASON_ASCENT_ORBIT)) };
    let score = |kick: f64| run_score(kick, q_cap);
    let (mut lo, mut hi) = ((k0 - 1.0).max(0.05), k0 + 1.0);
    let phi = 0.618_033_988_75;
    let f = |kd: f64| score(kd.to_radians()).unwrap_or(f64::NEG_INFINITY);
    let (mut c, mut d) = (hi - phi * (hi - lo), lo + phi * (hi - lo));
    let (mut fc, mut fd) = (f(c), f(d));
    for _ in 0..12 {
        iterations += 1;
        if fc > fd {
            hi = d;
            (d, fd) = (c, fc);
            c = hi - phi * (hi - lo);
            fc = f(c);
        } else {
            lo = c;
            (c, fc) = (d, fd);
            d = lo + phi * (hi - lo);
            fd = f(d);
        }
    }
    let kick_deg = if fc.max(fd) > f64::NEG_INFINITY {
        if fc > fd { c } else { d }
    } else {
        k0
    };
    let kick = kick_deg.to_radians();
    let converged = score(kick).is_some();

    let o = Sim { inp, kick, record: true }.run();
    if !o.inserted {
        return Ok(Err(REASON_ASCENT_ORBIT));
    }
    let mu = inp.mu_km3_s2;
    let (lat_c, r_pad) = ellipsoid_point(inp.lat_deg, inp.alt_km);
    let site = Vector3::new(lat_c.cos(), 0.0, lat_c.sin()) * r_pad;
    let up0 = site.normalize();
    let az = inp.azimuth_deg.to_radians();
    let east0 = Vector3::new(0.0, 1.0, 0.0);
    let north0 = up0.cross(&east0);
    let h = up0.cross(&(north0 * az.cos() + east0 * az.sin())).normalize();
    let along = h.cross(&up0);
    let rp = o.r - h * o.r.dot(&h);
    let angle = rp.dot(&along).atan2(rp.dot(&up0));
    let energy = o.v.norm_squared() / 2.0 - mu / o.r.norm();
    let sma = -mu / (2.0 * energy);
    let hv = o.r.cross(&o.v);
    let ecc = (1.0 - hv.norm_squared() / (mu * sma)).max(0.0).sqrt();
    let ve_upper = G0 * inp.upper.isp_s;
    let m_dry_upper = inp.upper.dry_mass_kg + inp.payload_kg;
    let liftoff_mass = Sim { inp, kick, record: false }.run_mass0();
    let (lift_thrust, _) = liftoff_thrust(inp);
    Ok(Ok(AscentResult {
        converged,
        iterations,
        pitch_kick_deg: kick.to_degrees(),
        duration_s: o.t,
        downrange_angle_rad: angle,
        downrange_km: angle * EARTH_EQUATORIAL_RADIUS_KM,
        insertion_altitude_km: o.r.norm() - EARTH_EQUATORIAL_RADIUS_KM,
        insertion_speed_km_s: o.v.norm(),
        perigee_km: sma * (1.0 - ecc) - EARTH_EQUATORIAL_RADIUS_KM,
        apogee_km: sma * (1.0 + ecc) - EARTH_EQUATORIAL_RADIUS_KM,
        plane_error_deg: hv.normalize().dot(&h).clamp(-1.0, 1.0).acos().to_degrees(),
        max_q_pa: o.max_q.0,
        max_q_t_s: o.max_q.1,
        max_acceleration_g: o.max_acc,
        max_q_above_limit: o.max_q.0 > inp.vehicle.max_q_limit_kpa * 1000.0,
        gravity_loss_m_s: o.losses.0,
        drag_loss_m_s: o.losses.1,
        steering_loss_m_s: o.losses.2,
        pressure_loss_m_s: o.p_loss,
        ideal_delta_v_m_s: o.ideal,
        initial_speed_m_s: o.v_start * 1000.0,
        upper_stage_used: o.upper_used,
        upper_stage_propellant_left_kg: o.upper_left,
        upper_stage_delta_v_left_m_s: ve_upper * ((m_dry_upper + o.upper_left) / m_dry_upper).ln(),
        mass_in_orbit_kg: o.mass,
        liftoff_mass_kg: liftoff_mass,
        liftoff_thrust_to_weight: lift_thrust / (liftoff_mass * G0),
        events: o.events,
        points: o.points,
        insertion_r: o.r,
        insertion_v: o.v,
        site_r: site,
        plane_normal: h,
    }))
}

/// Тяга на старте (у земли), Н, и её доля от ускорителей.
fn liftoff_thrust(inp: &AscentInput) -> (f64, f64) {
    let veh = inp.vehicle;
    let mut f = 0.0;
    let mut fb = 0.0;
    if let Some(b) = &veh.boosters {
        let fv = b.thrust_vac_at(0.0);
        fb = b.count as f64 * (fv - b.exit_area(b.peak_thrust_vac()) * P0).max(0.0);
        f += fb;
    }
    if let Some(s) = veh.stages.first() {
        let fv = s.thrust_vac_at(0.0);
        f += s.count as f64 * (fv - s.exit_area(s.peak_thrust_vac()) * P0).max(0.0);
    }
    (f, fb)
}

impl Sim<'_> {
    fn run_mass0(&self) -> f64 {
        let veh = self.inp.vehicle;
        let mut m = self.inp.payload_kg + self.inp.upper.dry_mass_kg + self.inp.upper.propellant_kg + veh.fairing_mass_kg;
        for s in &veh.stages {
            m += (s.dry_mass_kg + s.propellant_kg) * s.count as f64;
        }
        if let Some(b) = &veh.boosters {
            m += (b.dry_mass_kg + b.propellant_kg) * b.count as f64;
        }
        m
    }
}

impl AscentResult {
    /// Матрица перехода из местной системы в J2000: космодром переходит в `site_dir`, а плоскость
    /// подъёма — точно в плоскость опорной орбиты с нормалью `plane_normal` (J2000).
    pub fn to_j2000(&self, site_dir: &Vector3, plane_normal: &Vector3) -> anise::math::Matrix3 {
        let basis = |s: &Vector3, h: &Vector3| {
            let h = h.normalize();
            let s = (s - h * h.dot(s)).normalize();
            anise::math::Matrix3::from_columns(&[s, h.cross(&s), h])
        };
        basis(site_dir, plane_normal) * basis(&self.site_r, &self.plane_normal).transpose()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_atmosphere_reference_values() {
        // US Standard Atmosphere 1976.
        let (rho0, p0, a0) = atmosphere(0.0);
        assert!((rho0 - 1.225).abs() < 1e-3 && (p0 - 101_325.0).abs() < 1.0 && (a0 - 340.3).abs() < 0.2);
        let (rho11, p11, _) = atmosphere(11.0);
        assert!((p11 - 22_700.0).abs() < 150.0, "{p11}");
        assert!((rho11 - 0.3639).abs() < 0.003, "{rho11}");
        let (rho50, _, _) = atmosphere(50.0);
        assert!((rho50 - 1.027e-3).abs() < 3e-5, "{rho50}");
        let (rho80, _, _) = atmosphere(80.0);
        assert!((rho80 - 1.846e-5).abs() < 1e-6, "{rho80}");
        let (rho200, _, _) = atmosphere(200.0);
        assert!((rho200 - 2.541e-10).abs() < 1e-12);
        // Монотонно убывает.
        let mut prev = f64::INFINITY;
        for h in (0..1000).step_by(5) {
            let r = atmosphere(h as f64).0;
            assert!(r < prev || r == 0.0, "{h}");
            prev = r;
        }
    }

    #[test]
    fn drag_coefficient_peaks_transonic() {
        assert!(drag_coefficient(1.2) > drag_coefficient(0.5));
        assert!(drag_coefficient(1.2) > drag_coefficient(4.0));
        assert_eq!(drag_coefficient(10.0), 0.28);
    }

    #[test]
    fn thrust_profile_burns_all_propellant() {
        let s = StageSpec {
            name: "srb".into(),
            count: 1,
            dry_mass_kg: 1.0,
            propellant_kg: 100_000.0,
            thrust_vac_n: None,
            thrust_sl_n: None,
            isp_vac_s: 270.0,
            isp_sl_s: 245.0,
            burn_time_s: Some(100.0),
            thrust_profile: Some(vec![[0.0, 1.0], [0.5, 0.6], [1.0, 0.0]]),
            diameter_m: 3.0,
        };
        let burned: f64 = (0..10000).map(|i| s.thrust_vac_at(i as f64 * 0.01) / (G0 * s.isp_vac_s) * 0.01).sum();
        assert!((burned - 100_000.0).abs() < 200.0, "{burned}");
    }

    fn simple_two_stage() -> AscentVehicle {
        AscentVehicle {
            boosters: None,
            stages: vec![
                StageSpec {
                    name: "first".into(),
                    count: 1,
                    dry_mass_kg: 20_000.0,
                    propellant_kg: 400_000.0,
                    thrust_vac_n: Some(7.6e6),
                    thrust_sl_n: Some(6.8e6),
                    isp_vac_s: 311.0,
                    isp_sl_s: 282.0,
                    burn_time_s: None,
                    thrust_profile: None,
                    diameter_m: 3.7,
                },
                StageSpec {
                    name: "second".into(),
                    count: 1,
                    dry_mass_kg: 4_000.0,
                    propellant_kg: 100_000.0,
                    thrust_vac_n: Some(980_000.0),
                    thrust_sl_n: None,
                    isp_vac_s: 348.0,
                    isp_sl_s: 200.0,
                    burn_time_s: None,
                    thrust_profile: None,
                    diameter_m: 3.7,
                },
            ],
            upper_stage_in_ascent: false,
            upper_stage_diameter_m: 3.7,
            fairing_mass_kg: 1900.0,
            fairing_jettison_altitude_km: 110.0,
            vertical_rise_s: 10.0,
            stage_coast_s: 5.0,
            pitch_kick_deg: 3.0,
            max_q_limit_kpa: 45.0,
            source: None,
        }
    }

    #[test]
    fn two_stage_rocket_reaches_circular_orbit() {
        let veh = simple_two_stage();
        let inp = AscentInput {
            vehicle: &veh,
            upper: UpperStage { dry_mass_kg: 1000.0, propellant_kg: 3000.0, thrust_n: 20_000.0, isp_s: 320.0 },
            payload_kg: 8000.0,
            lat_deg: 28.5,
            alt_km: 0.0,
            azimuth_deg: 90.0,
            target_radius_km: EARTH_EQUATORIAL_RADIUS_KM + 185.0,
            mu_km3_s2: 398_600.4418,
        };
        let r = simulate(&inp).unwrap().expect("reaches orbit");
        assert!(r.converged, "{r:?}");
        assert!((r.perigee_km - 185.0).abs() < 2.0 && (r.apogee_km - 185.0).abs() < 2.0, "{} x {}", r.perigee_km, r.apogee_km);
        assert!(r.plane_error_deg < 0.05, "{}", r.plane_error_deg);
        assert!((300.0..900.0).contains(&r.duration_s), "{}", r.duration_s);
        assert!((500.0..3000.0).contains(&r.downrange_km), "{}", r.downrange_km);
        // Потери: гравитационные 1-2 км/с, сопротивление ~100 м/с.
        assert!((700.0..2500.0).contains(&r.gravity_loss_m_s), "{}", r.gravity_loss_m_s);
        assert!((20.0..400.0).contains(&r.drag_loss_m_s), "{}", r.drag_loss_m_s);
        assert!(r.max_q_pa > 15_000.0 && r.max_q_pa < 60_000.0, "{}", r.max_q_pa);
        assert!(r.liftoff_thrust_to_weight > 1.0);
        // Баланс скорости: идеальный разгон − потери = прирост скорости (с точностью до шага и J2).
        let balance = r.ideal_delta_v_m_s - r.gravity_loss_m_s - r.drag_loss_m_s - r.steering_loss_m_s - r.pressure_loss_m_s;
        let gained = r.insertion_speed_km_s * 1000.0 - r.initial_speed_m_s;
        assert!((balance - gained).abs() < 0.01 * gained, "balance {balance} vs gained {gained}");
        let names: Vec<_> = r.events.iter().map(|e| e.name).collect();
        for n in ["liftoff", "pitch_kick", "max_q", "stage_cutoff", "stage_ignition", "fairing_jettison", "insertion"] {
            assert!(names.contains(&n), "{n}: {names:?}");
        }
        assert!(!r.upper_stage_used && r.upper_stage_propellant_left_kg == 3000.0);
        // Тяжёлый груз не выводится.
        let heavy = AscentInput { payload_kg: 80_000.0, ..inp.clone() };
        assert_eq!(simulate(&heavy).unwrap(), Err(REASON_ASCENT_ORBIT));
    }
}

#[cfg(test)]
mod real_rockets {
    use super::*;
    use crate::launch_vehicle::{LaunchVehicleLibrary, LaunchVehicleParams};

    #[test]
    fn catalog_rockets_reach_orbit() {
        let lib = LaunchVehicleLibrary::presets_only().unwrap();
        for (id, lat, az, payload) in [
            ("sls_block1", 28.6, 90.0, 26_000.0),
            ("long_march_5", 19.6, 100.0, 8_000.0),
            ("ariane_64", 5.2, 90.0, 8_000.0),
            ("ariane_62", 5.2, 90.0, 2_800.0),
            ("soyuz_21b_fregat", 51.9, 60.0, 1_800.0),
        ] {
            let rk = lib.resolve(id, &LaunchVehicleParams::default()).unwrap();
            let veh = rk.config.ascent.clone().expect("ascent data");
            let inp = AscentInput {
                vehicle: &veh,
                upper: rk.config.upper_stage(),
                payload_kg: payload,
                lat_deg: lat,
                alt_km: 0.0,
                azimuth_deg: az,
                target_radius_km: EARTH_EQUATORIAL_RADIUS_KM + 185.0,
                mu_km3_s2: 398_600.4418,
            };
            let r = simulate(&inp).unwrap();
            match &r {
                Ok(r) => println!(
                    "{id:18} conv {} it {:2} kick {:5.2} T {:6.1}s down {:5.0}km {:.1}x{:.1} q {:5.1}kPa@{:4.0}s acc {:4.1}g losses g{:4.0} d{:4.0} s{:4.0} p{:4.0} upper used {} left {:6.0}kg dv_left {:5.0} TW {:.2} m0 {:.0}t\n    {:?}",
                    r.converged, r.iterations, r.pitch_kick_deg, r.duration_s, r.downrange_km,
                    r.perigee_km, r.apogee_km, r.max_q_pa / 1000.0, r.max_q_t_s, r.max_acceleration_g,
                    r.gravity_loss_m_s, r.drag_loss_m_s, r.steering_loss_m_s, r.pressure_loss_m_s,
                    r.upper_stage_used, r.upper_stage_propellant_left_kg, r.upper_stage_delta_v_left_m_s,
                    r.liftoff_thrust_to_weight, r.liftoff_mass_kg / 1000.0,
                    r.events.iter().map(|e| format!("{}@{:.0}", e.name, e.t_s)).collect::<Vec<_>>()
                ),
                Err(e) => println!("{id}: {e}"),
            }
            let r = r.expect("orbit");
            assert!(r.converged, "{id}");
            // Орбита 185 км, плоскость — заданная.
            assert!((r.perigee_km - 185.0).abs() < 2.0 && (r.apogee_km - 185.0).abs() < 2.0, "{id}");
            assert!(r.plane_error_deg < 0.05, "{id}: {}", r.plane_error_deg);
            // Баланс скорости.
            let balance = r.ideal_delta_v_m_s - r.gravity_loss_m_s - r.drag_loss_m_s - r.steering_loss_m_s - r.pressure_loss_m_s;
            let gained = r.insertion_speed_km_s * 1000.0 - r.initial_speed_m_s;
            assert!((balance - gained).abs() < 0.01 * gained, "{id}: {balance} vs {gained}");
            // Разгонной ступени хватает на разгон к Луне с груза около паспортного.
            assert!(r.upper_stage_delta_v_left_m_s > 2900.0, "{id}: {}", r.upper_stage_delta_v_left_m_s);
            // Времена событий близки к опубликованным.
            let at = |name: &str| r.events.iter().find(|e| e.name == name).map(|e| e.t_s);
            let near = |name: &str, lo: f64, hi: f64| {
                let t = at(name).unwrap_or(f64::NAN);
                assert!(t >= lo && t <= hi, "{id}: {name} at {t}, expected {lo}..{hi}");
            };
            match id {
                // Artemis I: max-Q 1:10, отделение ускорителей 2:12, выключение центральной ступени 8:03.
                "sls_block1" => {
                    near("max_q", 50.0, 80.0);
                    near("booster_separation", 115.0, 140.0);
                    assert!((430.0..500.0).contains(&r.duration_s), "{}", r.duration_s);
                }
                // Soyuz-2.1b: боковые блоки ~118-120 с, центральный блок ~280-290 с.
                "soyuz_21b_fregat" => {
                    near("booster_separation", 110.0, 130.0);
                    near("stage_cutoff", 270.0, 300.0);
                }
                // Long March 5: ускорители 173 с, первая ступень 492 с.
                "long_march_5" => {
                    near("booster_separation", 165.0, 185.0);
                    near("stage_cutoff", 480.0, 510.0);
                }
                _ => {}
            }
        }
    }
}
