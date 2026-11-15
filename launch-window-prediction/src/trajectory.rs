//! Перелёт Земля -> Луна: связь даты старта и даты прилёта.
//!
//! «Обратное» направление означает обратное отображение дат (дата прилёта -> дата старта),
//! а не полёт назад к Земле. Новые типы траекторий добавляются реализацией трейта [`Trajectory`].
//!
//! # Модель Direct Transfer
//!
//! 1. Старт с космодрома на круговую опорную орбиту, плоскость которой проходит через
//!    положение Луны в момент прилёта. Наклонение — наименьшее, которое допускают широта
//!    космодрома, склонение Луны и ограничения космодрома (азимуты, коридоры наклонений).
//! 2. Старт возможен, когда космодром из-за вращения Земли проходит через эту плоскость
//!    (дважды в сутки), а азимут пуска в этот момент входит в разрешённый сектор.
//!    Вращение Земли берётся из ядер ANISE (ITRF93).
//! 3. Полёт по опорной орбите от `coast_min` до `coast_max`, затем разгон к Луне (TLI).
//! 4. Перелёт до Луны — задача Ламберта в поле Земли по реальному положению Луны (DE440).
//! 5. У Луны — торможение на круговую окололунную орбиту.
//! 6. Из всех вариантов выбирается минимум суммарного Δv (разгон + торможение).
//!
//! Упрощения (первое приближение, «patched conic»): притяжение Луны на перелёте и Солнца
//! не учитывается, участок выведения на опорную орбиту считается мгновенным, момент прилёта —
//! это момент встречи с центром Луны в этой модели (с точностью до нескольких часов совпадает
//! с выходом на окололунную орбиту).

use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

use anise::constants::frames::{EARTH_ITRF93, EARTH_J2000, IAU_EARTH_FRAME, MOON_J2000};
use anise::math::Vector3;
use anise::prelude::{Almanac, Frame};
use anise::time::{Duration, Epoch, Unit};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};
use crate::lambert;
use crate::lunar_landing::{DescentConfig, LandingInfo, LandingTarget, LandingWindowConfig, find_landing};
use crate::refine::{NumericalConfig, Refined};
use crate::runlog::{ser_duration_s, ser_epoch, ser_vec3};
use crate::sky::LunarSite;
use crate::spaceport::Spaceport;

/// Гравитационные параметры по умолчанию, км³/с² (если в ядрах их нет).
const MU_EARTH_DEFAULT: f64 = 398_600.435_436;
const MU_MOON_DEFAULT: f64 = 4_902.800_066;
/// Экваториальный радиус Земли и сжатие WGS84.
const EARTH_RADIUS_KM: f64 = 6_378.137;
const EARTH_FLATTENING: f64 = 1.0 / 298.257_223_563;
/// Радиус Луны, км (опорная сфера LOLA).
const MOON_RADIUS_KM: f64 = 1_737.4;
/// Скорость вращения Земли на экваторе, км/с.
const EARTH_EQUATOR_SPEED_KM_S: f64 = 0.465_1;
/// Наибольшее наклонение опорной орбиты (только прямые орбиты), градусы.
const MAX_INCLINATION_DEG: f64 = 89.0;
/// Шаг перебора наклонений, если космодром задаёт диапазон или ничего не задаёт, градусы.
const INCLINATION_STEP_DEG: f64 = 1.0;
/// Δv спуска с низкой окололунной орбиты на поверхность по умолчанию, м/с
/// (номинальная проектная величина, AIAA-2008-6939).
pub const DEFAULT_DESCENT_DELTA_V_M_S: f64 = 2053.0;
const REASON_BUDGET: &str = "the spacecraft's delta-v budget is not enough for any transfer arriving at this time";
const REASON_WINDOW: &str = "launch window does not overlap this arrival time";

// ---------------------------------------------------------------------------
// Общий интерфейс
// ---------------------------------------------------------------------------

/// Какая дата задана.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
#[serde(tag = "kind", content = "day_start_utc", rename_all = "snake_case")]
pub enum DateConstraint {
    /// Старт в эти сутки (UTC), начиная с указанного момента полуночи.
    LaunchDay(#[serde(serialize_with = "ser_epoch")] Epoch),
    /// Прилёт в эти сутки (UTC).
    ArrivalDay(#[serde(serialize_with = "ser_epoch")] Epoch),
}

/// Опорная орбита Земли.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ParkingOrbit {
    pub altitude_km: f64,
    pub radius_km: f64,
    pub circular_speed_km_s: f64,
    pub period_s: f64,
    pub inclination_deg: f64,
    /// Долгота восходящего узла в J2000, градусы.
    pub raan_deg: f64,
    /// Единичная нормаль плоскости в J2000.
    #[serde(serialize_with = "ser_vec3")]
    pub plane_normal: Vector3,
}

/// Старт с космодрома.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LaunchInfo {
    #[serde(serialize_with = "ser_epoch")]
    pub epoch: Epoch,
    pub azimuth_deg: f64,
    /// Космодром движется на север (восходящий проход плоскости) или на юг.
    pub ascending_pass: bool,
    /// Прибавка скорости от вращения Земли по направлению пуска, км/с (справочно, в Δv не входит).
    pub earth_rotation_gain_km_s: f64,
    /// Направление на космодром в момент старта, J2000.
    #[serde(serialize_with = "ser_vec3")]
    pub site_direction_j2000: Vector3,
    /// Земная система, по которой считалось вращение Земли.
    pub earth_frame: &'static str,
}

/// Подъём от стартового стола до выхода на опорную орбиту.
///
/// Если для ракеты есть данные ступеней — численное моделирование ([`crate::ascent`]).
/// Иначе упрощённая модель: траектория подъёма лежит в плоскости опорной орбиты. Высота растёт по закону
/// h = H·(6τ² − 8τ³ + 3τ⁴) (вертикальный старт, горизонтальный выход на орбиту), угол
/// от космодрома — кубически, от скорости вращения Земли на старте до орбитальной на выходе,
/// τ = t / длительность. Это кинематическая модель для траектории и дат, а не расчёт тяги.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct AscentInfo {
    /// "simulation" — моделирование ракеты по ступеням ([`crate::ascent`]); "kinematic" — упрощённо.
    pub model: &'static str,
    /// Ракета, по которой моделировался подъём.
    pub rocket: Option<String>,
    /// Почему подъём упрощённый (если так).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub model_note: Option<String>,
    pub duration_s: f64,
    pub downrange_km: f64,
    /// Угол от космодрома до точки выхода на орбиту (из центра Земли), градусы.
    pub downrange_angle_deg: f64,
    #[serde(serialize_with = "ser_epoch")]
    pub insertion_epoch: Epoch,
    /// Направление на космодром в момент старта, спроецированное на плоскость орбиты (J2000).
    #[serde(serialize_with = "ser_vec3")]
    pub launch_direction_in_plane: Vector3,
    #[serde(serialize_with = "ser_vec3")]
    pub insertion_position_km: Vector3,
    #[serde(serialize_with = "ser_vec3")]
    pub insertion_velocity_km_s: Vector3,
    /// Радиус точки старта (Земля + высота космодрома), км.
    pub pad_radius_km: f64,
    /// Скорость вращения Земли вдоль направления пуска, рад/с (угловая).
    pub earth_rotation_rate_along_rad_s: f64,
    /// Результат моделирования подъёма (если моделировался).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub simulation: Option<crate::ascent::AscentResult>,
    /// Переход из местной системы моделирования в J2000.
    #[serde(skip)]
    pub to_j2000: Option<anise::math::Matrix3>,
}

impl AscentInfo {
    /// Положение и скорость на подъёме через `t` секунд после старта (J2000, км, км/с).
    pub fn state_at(&self, t: f64, plane_normal: &Vector3, parking_radius_km: f64) -> (Vector3, Vector3) {
        let tt = self.duration_s.max(1e-9);
        let tau = (t / tt).clamp(0.0, 1.0);
        let h_top = parking_radius_km - self.pad_radius_km;
        let alt = h_top * tau * tau * (6.0 - 8.0 * tau + 3.0 * tau * tau);
        let alt_rate = h_top * 12.0 * tau * (1.0 - tau) * (1.0 - tau) / tt;
        let theta = self.downrange_angle_deg.to_radians();
        let w0 = self.earth_rotation_rate_along_rad_s;
        let n = self.insertion_velocity_km_s.norm() / parking_radius_km;
        let b = (3.0 * theta - (2.0 * w0 + n) * tt) / (tt * tt);
        let c = ((w0 + n) * tt - 2.0 * theta) / (tt * tt * tt);
        let tc = tau * tt;
        let phi = w0 * tc + b * tc * tc + c * tc * tc * tc;
        let phi_rate = w0 + 2.0 * b * tc + 3.0 * c * tc * tc;
        let l = self.launch_direction_in_plane;
        let along = plane_normal.cross(&l);
        let u = l * phi.cos() + along * phi.sin();
        let u_dot = (along * phi.cos() - l * phi.sin()) * phi_rate;
        let radius = self.pad_radius_km + alt;
        (u * radius, u * alt_rate + u_dot * radius)
    }
}

/// Разгон к Луне (TLI) и орбита перелёта.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct DepartureInfo {
    #[serde(serialize_with = "ser_epoch")]
    pub epoch: Epoch,
    /// Время на опорной орбите от старта до разгона.
    #[serde(serialize_with = "ser_duration_s")]
    pub coast_s: Duration,
    /// Угол, пройденный по опорной орбите до разгона, градусы.
    pub coast_angle_deg: f64,
    #[serde(serialize_with = "ser_vec3")]
    pub position_km: Vector3,
    /// Скорость на опорной орбите до разгона, км/с.
    #[serde(serialize_with = "ser_vec3")]
    pub velocity_parking_km_s: Vector3,
    /// Скорость после разгона (начало перелёта), км/с.
    #[serde(serialize_with = "ser_vec3")]
    pub velocity_transfer_km_s: Vector3,
    #[serde(serialize_with = "ser_vec3")]
    pub delta_v_km_s: Vector3,
    pub delta_v_magnitude_km_s: f64,
    /// Угол наклона траектории после разгона к местному горизонту, градусы.
    pub flight_path_angle_deg: f64,
    /// Характеристическая энергия C3 = v² - 2μ/r, км²/с².
    pub c3_km2_s2: f64,
    pub transfer_semi_major_axis_km: f64,
    pub transfer_eccentricity: f64,
    pub transfer_perigee_km: f64,
    /// `None` для незамкнутой орбиты.
    pub transfer_apogee_km: Option<f64>,
    pub transfer_angle_deg: f64,
    #[serde(serialize_with = "ser_duration_s")]
    pub time_of_flight_s: Duration,
}

/// Прилёт к Луне и торможение.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ArrivalInfo {
    #[serde(serialize_with = "ser_epoch")]
    pub epoch: Epoch,
    /// Положение и скорость Луны относительно Земли, J2000.
    #[serde(serialize_with = "ser_vec3")]
    pub moon_position_km: Vector3,
    #[serde(serialize_with = "ser_vec3")]
    pub moon_velocity_km_s: Vector3,
    pub moon_distance_km: f64,
    pub moon_declination_deg: f64,
    /// Скорость аппарата в конце перелёта относительно Земли, км/с.
    #[serde(serialize_with = "ser_vec3")]
    pub velocity_arrival_km_s: Vector3,
    /// Скорость подлёта относительно Луны (v∞), км/с.
    #[serde(serialize_with = "ser_vec3")]
    pub v_infinity_km_s: Vector3,
    pub v_infinity_magnitude_km_s: f64,
    pub lunar_orbit_altitude_km: f64,
    pub lunar_orbit_radius_km: f64,
    pub lunar_orbit_circular_speed_km_s: f64,
    /// Скорость в периселении гиперболы подлёта, км/с.
    pub periselene_speed_km_s: f64,
    /// Δv торможения на окололунную орбиту, км/с.
    pub delta_v_loi_km_s: f64,
}

/// Как выбиралось наклонение.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct InclinationChoice {
    /// Максимальное склонение космодрома за окно старта, градусы.
    pub site_max_declination_deg: f64,
    pub moon_declination_deg: f64,
    pub margin_deg: f64,
    /// Наименьшее наклонение по геометрии (без ограничений космодрома), градусы.
    pub required_min_deg: f64,
    /// Наклонения, которые перебирались для момента прилёта (по возрастанию).
    pub candidates_deg: Vec<f64>,
    pub chosen_deg: f64,
}

/// Лучший найденный вариант перелёта со всеми данными.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TransferSolution {
    pub trajectory: &'static str,
    pub spaceport_id: String,
    pub launch: LaunchInfo,
    /// Подъём до опорной орбиты.
    pub ascent: AscentInfo,
    pub parking_orbit: ParkingOrbit,
    pub departure: DepartureInfo,
    pub arrival: ArrivalInfo,
    pub inclination: InclinationChoice,
    pub delta_v_total_km_s: f64,
    pub mu_earth_km3_s2: f64,
    pub mu_moon_km3_s2: f64,
    /// Проверка по бюджету лунного аппарата (если он задан).
    pub budget: Option<BudgetCheck>,
    /// Посадка в заданную точку (если она задана).
    pub landing: Option<LandingInfo>,
    /// Уточнение численным интегрированием ([`crate::refine`]); заполняется после поиска.
    pub refined: Option<Refined>,
    /// Почему уточнить не удалось (если не удалось).
    pub refine_failure: Option<String>,
    /// Почему уточнение не запускалось (уточняется только лучший вариант).
    pub refine_skipped: Option<String>,
}

/// Лучший вариант для одного момента прилёта (строка перебора).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ScanPoint {
    #[serde(serialize_with = "ser_epoch")]
    pub arrival: Epoch,
    pub feasible: bool,
    /// Почему вариантов нет (если `feasible = false`).
    pub reason: Option<&'static str>,
    #[serde(serialize_with = "crate::runlog::ser_opt_epoch")]
    pub launch: Option<Epoch>,
    pub coast_min: Option<f64>,
    pub time_of_flight_h: Option<f64>,
    pub inclination_deg: Option<f64>,
    pub azimuth_deg: Option<f64>,
    pub delta_v_tli_km_s: Option<f64>,
    pub delta_v_loi_km_s: Option<f64>,
    pub delta_v_total_km_s: Option<f64>,
    #[serde(serialize_with = "crate::runlog::ser_opt_epoch")]
    pub landing: Option<Epoch>,
}

/// Счётчики перебора.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct SearchStats {
    pub arrival_epochs: u64,
    pub inclinations_tried: u64,
    pub planes_tried: u64,
    pub plane_crossings: u64,
    pub crossings_rejected_azimuth: u64,
    /// Пуски, отброшенные потому, что ракета не выводит груз на опорную орбиту (моделирование подъёма).
    pub crossings_rejected_ascent: u64,
    pub candidates_rejected_budget: u64,
    pub candidates_rejected_tli_limit: u64,
    pub arrivals_rejected_landing: u64,
    pub lambert_solves: u64,
    pub lambert_failures: u64,
    pub candidates_evaluated: u64,
}

/// Полный отчёт о поиске перелёта.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TransferReport {
    pub trajectory: &'static str,
    pub spaceport: Spaceport,
    pub constraint: DateConstraint,
    pub config: DirectTransferConfig,
    #[serde(serialize_with = "ser_epoch")]
    pub launch_window_start: Epoch,
    #[serde(serialize_with = "ser_epoch")]
    pub launch_window_end: Epoch,
    #[serde(serialize_with = "ser_epoch")]
    pub arrival_window_start: Epoch,
    #[serde(serialize_with = "ser_epoch")]
    pub arrival_window_end: Epoch,
    pub earth_frame: &'static str,
    pub site_max_declination_deg: f64,
    pub budget: Option<DeltaVBudget>,
    pub landing_site: Option<LunarSite>,
    pub landing_window: Option<LandingWindowConfig>,
    pub solution: Option<TransferSolution>,
    /// Главная причина, если вариантов нет (самая содержательная из причин перебора).
    pub no_solution_reason: Option<&'static str>,
    pub coarse_scan: Vec<ScanPoint>,
    pub fine_scan: Vec<ScanPoint>,
    pub stats: SearchStats,
    /// Варианты, которые проверялись численным уточнением по порядку (если первый не подошёл).
    pub refine_attempts: Vec<RefineAttempt>,
}

/// Один вариант, проверенный численным уточнением.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct RefineAttempt {
    /// Момент прилёта из поиска по формулам.
    #[serde(serialize_with = "ser_epoch")]
    pub arrival: Epoch,
    #[serde(serialize_with = "crate::runlog::ser_opt_epoch")]
    pub launch: Option<Epoch>,
    pub accepted: bool,
    /// Почему не подошёл (если не подошёл).
    pub reason: Option<String>,
}

/// Тип траектории. Реализации должны быть потокобезопасными: поиск идёт параллельно.
pub trait Trajectory: Send + Sync {
    fn name(&self) -> &'static str;

    /// Лучший вариант перелёта с космодрома при заданной дате старта или прилёта
    /// вместе со всеми данными поиска. `report.solution = None` — в эти сутки вариантов нет.
    ///
    /// С `budget` учитываются возможности лунного аппарата: отбрасываются варианты,
    /// на которые не хватает его Δv, и выбор идёт по `budget.objective`.
    ///
    /// С `landing` учитывается точка посадки: прилёт годится, только если после него
    /// аппарат может сесть в точку (проход орбиты над ней при подходящем освещении и связи).
    /// Тогда `DateConstraint::ArrivalDay` означает сутки посадки.
    fn solve(
        &self,
        almanac: &Almanac,
        spaceport: &Spaceport,
        constraint: DateConstraint,
        budget: Option<&DeltaVBudget>,
        landing: Option<&LandingTarget>,
    ) -> Result<TransferReport>;

    /// Лучший вариант с прилётом точно в момент `arrival` при тех же ограничениях, что у
    /// [`Trajectory::solve`]; `None` — в этот момент вариантов нет. Нужен, когда численное
    /// уточнение отвергло лучший вариант и надо проверить другие моменты прилёта из перебора.
    fn solve_at(
        &self,
        _almanac: &Almanac,
        _spaceport: &Spaceport,
        _constraint: DateConstraint,
        _budget: Option<&DeltaVBudget>,
        _landing: Option<&LandingTarget>,
        _arrival: Epoch,
    ) -> Result<Option<TransferSolution>> {
        Ok(None)
    }
}

// ---------------------------------------------------------------------------
// Параметры
// ---------------------------------------------------------------------------

/// Файл параметров перелёта (`data/transfer.toml`).
#[derive(Debug, Clone, Default, PartialEq, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct TransferConfig {
    pub direct_transfer: DirectTransferConfig,
    pub mission: MissionConfig,
    /// Условия посадки в заданную точку (если задана точка посадки).
    pub landing_window: LandingWindowConfig,
    /// Спуск с окололунной орбиты в точку посадки.
    pub descent: DescentConfig,
    /// Уточнение найденного перелёта численным интегрированием.
    pub numerical: NumericalConfig,
}

/// Манёвр, который может выполнять лунный аппарат собственным топливом.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Maneuver {
    /// Разгон к Луне с опорной орбиты.
    Tli,
    /// Торможение на окололунную орбиту.
    Loi,
    /// Спуск с окололунной орбиты на поверхность.
    Descent,
}

/// Что выбирать среди допустимых перелётов.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Objective {
    /// Наименьший суммарный Δv.
    MinDeltaV,
    /// Самый быстрый перелёт (с точностью до часа), на который хватает Δv аппарата
    /// и разгонной ступени; среди одинаково быстрых — самый экономичный.
    Fastest,
}

/// Как аппарат участвует в миссии (секция `[mission]`).
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct MissionConfig {
    /// Манёвры, которые выполняет лунный аппарат своим топливом.
    pub spacecraft_maneuvers: Vec<Maneuver>,
    /// Δv спуска с окололунной орбиты на поверхность, м/с. Не задан — считается моделированием
    /// торможения по параметрам аппарата ([`crate::lunar_landing::DescentModel::plan`]).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub descent_delta_v_m_s: Option<f64>,
    /// Запас Δv сверх требуемого, проценты.
    pub reserve_percent: f64,
    /// Критерий выбора перелёта, если задан аппарат.
    pub objective: Objective,
}

impl Default for MissionConfig {
    fn default() -> Self {
        Self {
            spacecraft_maneuvers: vec![Maneuver::Loi, Maneuver::Descent],
            descent_delta_v_m_s: None,
            reserve_percent: 5.0,
            objective: Objective::Fastest,
        }
    }
}

impl MissionConfig {
    pub fn validate(&self) -> Result<()> {
        let bad = |m: &str| Err(Error::Config(format!("mission: {m}")));
        if self.spacecraft_maneuvers.is_empty() {
            return bad("spacecraft_maneuvers must not be empty");
        }
        if self.descent_delta_v_m_s.is_some_and(|d| !(d >= 0.0 && d <= 10_000.0)) {
            return bad("descent_delta_v_m_s must be in [0, 10000]");
        }
        if !(self.reserve_percent >= 0.0 && self.reserve_percent <= 100.0) {
            return bad("reserve_percent must be in [0, 100]");
        }
        Ok(())
    }
}

/// Бюджет Δv лунного аппарата для поиска перелёта.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct DeltaVBudget {
    pub spacecraft_name: String,
    /// Характеристическая скорость аппарата (Циолковский), м/с.
    pub available_m_s: f64,
    pub maneuvers: Vec<Maneuver>,
    pub descent_delta_v_m_s: f64,
    pub reserve_percent: f64,
    pub objective: Objective,
}

impl DeltaVBudget {
    pub fn new(spacecraft_name: &str, available_m_s: f64, mission: &MissionConfig) -> Result<Self> {
        mission.validate()?;
        if !(available_m_s > 0.0 && available_m_s.is_finite()) {
            return Err(Error::Config(format!("spacecraft delta-v must be > 0 (got {available_m_s})")));
        }
        Ok(Self {
            spacecraft_name: spacecraft_name.to_string(),
            available_m_s,
            maneuvers: mission.spacecraft_maneuvers.clone(),
            // Без моделирования — номинальная проектная величина (AIAA-2008-6939, Georgia Tech).
            descent_delta_v_m_s: mission.descent_delta_v_m_s.unwrap_or(DEFAULT_DESCENT_DELTA_V_M_S),
            reserve_percent: mission.reserve_percent,
            objective: mission.objective,
        })
    }

    /// Требуемый Δv аппарата без запаса, м/с.
    fn required_m_s(&self, dv_tli_km_s: f64, dv_loi_km_s: f64) -> f64 {
        self.maneuvers
            .iter()
            .map(|m| match m {
                Maneuver::Tli => dv_tli_km_s * 1000.0,
                Maneuver::Loi => dv_loi_km_s * 1000.0,
                Maneuver::Descent => self.descent_delta_v_m_s,
            })
            .sum()
    }

    fn with_reserve(&self, required_m_s: f64) -> f64 {
        required_m_s * (1.0 + self.reserve_percent / 100.0)
    }

    /// Проверка бюджета для заданных манёвров. `descent_m_s` — уточнённый Δv спуска (если есть).
    pub fn check(
        &self,
        dv_tli_km_s: f64,
        dv_loi_km_s: f64,
        descent_m_s: Option<f64>,
        time_of_flight_at_lower_limit: bool,
    ) -> BudgetCheck {
        let has = |m: Maneuver| self.maneuvers.contains(&m);
        let descent = descent_m_s.unwrap_or(self.descent_delta_v_m_s);
        let required = self.required_m_s(dv_tli_km_s, dv_loi_km_s)
            + if has(Maneuver::Descent) { descent - self.descent_delta_v_m_s } else { 0.0 };
        let with_reserve = self.with_reserve(required);
        BudgetCheck {
            spacecraft_name: self.spacecraft_name.clone(),
            objective: self.objective,
            maneuvers: self.maneuvers.clone(),
            tli_m_s: has(Maneuver::Tli).then_some(dv_tli_km_s * 1000.0),
            loi_m_s: has(Maneuver::Loi).then_some(dv_loi_km_s * 1000.0),
            descent_m_s: has(Maneuver::Descent).then_some(descent),
            required_m_s: required,
            reserve_percent: self.reserve_percent,
            required_with_reserve_m_s: with_reserve,
            available_m_s: self.available_m_s,
            margin_m_s: self.available_m_s - with_reserve,
            time_of_flight_at_lower_limit,
        }
    }

    fn fits(&self, dv_tli_km_s: f64, dv_loi_km_s: f64) -> bool {
        self.with_reserve(self.required_m_s(dv_tli_km_s, dv_loi_km_s)) <= self.available_m_s
    }
}

/// Проверка бюджета для найденного перелёта.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct BudgetCheck {
    pub spacecraft_name: String,
    pub objective: Objective,
    pub maneuvers: Vec<Maneuver>,
    pub tli_m_s: Option<f64>,
    pub loi_m_s: Option<f64>,
    pub descent_m_s: Option<f64>,
    pub required_m_s: f64,
    pub reserve_percent: f64,
    pub required_with_reserve_m_s: f64,
    pub available_m_s: f64,
    pub margin_m_s: f64,
    /// Время перелёта в пределах часа от нижней границы `tof_min_days` (быстрее не позволяют настройки).
    pub time_of_flight_at_lower_limit: bool,
}

impl TransferConfig {
    pub fn parse(text: &str, origin: &Path) -> Result<Self> {
        let cfg: Self = toml::from_str(text)
            .map_err(|source| Error::Toml { path: origin.to_path_buf(), source })?;
        cfg.direct_transfer.validate()?;
        cfg.mission.validate()?;
        cfg.landing_window.validate()?;
        cfg.descent.validate()?;
        cfg.numerical.validate()?;
        Ok(cfg)
    }

    pub async fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = tokio::fs::read_to_string(path)
            .await
            .map_err(|source| Error::Io { path: path.to_path_buf(), source })?;
        Self::parse(&text, path)
    }
}

/// Параметры Direct Transfer.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct DirectTransferConfig {
    /// Высота круговой опорной орбиты над экватором Земли, км.
    pub parking_altitude_km: f64,
    /// Высота круговой окололунной орбиты после торможения, км.
    pub lunar_orbit_altitude_km: f64,
    /// Допустимое время перелёта от разгона до Луны, сутки.
    pub tof_min_days: f64,
    pub tof_max_days: f64,
    /// Допустимое время на опорной орбите от старта до разгона, минуты.
    pub coast_min_min: f64,
    pub coast_max_min: f64,
    /// Запас наклонения опорной орбиты над минимально возможным, градусы.
    pub inclination_margin_deg: f64,
    /// Шаг грубого перебора момента прилёта, минуты (потом уточняется до 1 минуты).
    pub search_step_min: f64,
    /// Наибольший Δv разгона к Луне, который может дать разгонная ступень, км/с.
    /// Варианты с большим разгоном отбрасываются.
    pub max_tli_delta_v_km_s: f64,
    /// Подъём от стартового стола до выхода на опорную орбиту: длительность, минуты.
    pub ascent_duration_min: f64,
    /// Дальность от космодрома до точки выхода на опорную орбиту по поверхности, км.
    pub ascent_downrange_km: f64,
}

impl Default for DirectTransferConfig {
    fn default() -> Self {
        Self {
            parking_altitude_km: 185.0,
            lunar_orbit_altitude_km: 100.0,
            tof_min_days: 2.5,
            tof_max_days: 5.5,
            coast_min_min: 20.0,
            coast_max_min: 180.0,
            inclination_margin_deg: 0.5,
            search_step_min: 30.0,
            max_tli_delta_v_km_s: 3.3,
            // Типичный подъём к Луне: Saturn V — 11 мин 39 с до выхода на орбиту,
            // SLS — около 8 мин до отделения центральной ступени. Дальность ~2000 км.
            ascent_duration_min: 10.0,
            ascent_downrange_km: 2000.0,
        }
    }
}

impl DirectTransferConfig {
    pub fn validate(&self) -> Result<()> {
        let bad = |m: &str| Err(Error::Config(format!("direct_transfer: {m}")));
        if !(self.parking_altitude_km >= 120.0 && self.parking_altitude_km <= 2000.0) {
            return bad("parking_altitude_km must be in [120, 2000]");
        }
        if !(self.lunar_orbit_altitude_km >= 10.0 && self.lunar_orbit_altitude_km <= 10_000.0) {
            return bad("lunar_orbit_altitude_km must be in [10, 10000]");
        }
        if !(self.tof_min_days >= 1.0 && self.tof_max_days > self.tof_min_days && self.tof_max_days <= 10.0) {
            return bad("need 1 <= tof_min_days < tof_max_days <= 10");
        }
        if !(self.coast_min_min >= 0.0 && self.coast_max_min > self.coast_min_min && self.coast_max_min <= 24.0 * 60.0) {
            return bad("need 0 <= coast_min_min < coast_max_min <= 1440");
        }
        if !(self.inclination_margin_deg >= 0.0 && self.inclination_margin_deg <= 10.0) {
            return bad("inclination_margin_deg must be in [0, 10]");
        }
        if !(self.search_step_min >= 1.0 && self.search_step_min <= 360.0) {
            return bad("search_step_min must be in [1, 360]");
        }
        if !(self.max_tli_delta_v_km_s >= 3.0 && self.max_tli_delta_v_km_s <= 15.0) {
            return bad("max_tli_delta_v_km_s must be in [3, 15]");
        }
        if !(self.ascent_duration_min >= 0.0 && self.ascent_duration_min <= 30.0) {
            return bad("ascent_duration_min must be in [0, 30]");
        }
        if !(self.ascent_downrange_km >= 0.0 && self.ascent_downrange_km <= 5000.0) {
            return bad("ascent_downrange_km must be in [0, 5000]");
        }
        if self.coast_min_min < self.ascent_duration_min {
            return bad("coast_min_min (time from launch to TLI) must be at least ascent_duration_min");
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Direct Transfer
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Default)]
pub struct DirectTransfer {
    pub config: DirectTransferConfig,
    /// Ракета для моделирования подъёма по космодромам (id космодрома -> ракета).
    /// Для космодрома без ракеты подъём упрощённый (длительность и дальность из конфигурации).
    pub ascent: std::collections::HashMap<String, AscentSetup>,
    ascent_cache: std::sync::Arc<std::sync::Mutex<std::collections::HashMap<(String, i64), AscentSlot>>>,
    /// Папка кэша узлов подъёма на диске (между запусками); `None` — без кэша.
    ascent_cache_dir: Option<std::path::PathBuf>,
}

/// Ракета и груз для моделирования подъёма.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct AscentSetup {
    pub rocket_id: String,
    pub rocket_name: String,
    pub vehicle: crate::ascent::AscentVehicle,
    pub upper: crate::ascent::UpperStage,
    /// Груз над разгонной ступенью (лунный аппарат), кг.
    pub payload_kg: f64,
}

/// Узел сетки азимутов подъёма: длительность и угол выхода на орбиту или причина отказа.
type AscentNode = std::result::Result<(f64, f64), &'static str>;
/// Ячейка кэша: узел сетки считается один раз, даже если его одновременно запросили несколько потоков.
type AscentSlot = std::sync::Arc<std::sync::OnceLock<std::result::Result<AscentNode, String>>>;

/// Файл кэша подъёма: ключ (входные данные) и узлы (индекс азимута -> длительность и угол или отказ).
#[derive(Serialize, serde::Deserialize)]
struct AscentCacheFile {
    key: String,
    nodes: std::collections::BTreeMap<i64, Option<(f64, f64)>>,
}

/// Хэш FNV-1a (стабильный между версиями Rust, в отличие от `DefaultHasher`).
fn fnv1a(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325u64, |h, &b| (h ^ b as u64).wrapping_mul(0x0000_0100_0000_01b3))
}

/// Шаг сетки азимутов для поиска (между узлами — линейная интерполяция), градусы.
const ASCENT_AZIMUTH_STEP_DEG: f64 = 5.0;

impl DirectTransfer {
    pub fn new(config: DirectTransferConfig) -> Result<Self> {
        config.validate()?;
        Ok(Self { config, ascent: Default::default(), ascent_cache: Default::default(), ascent_cache_dir: None })
    }

    /// Хранить узлы сетки подъёма в папке `dir`, чтобы не пересчитывать их при следующих запусках.
    pub fn with_ascent_cache_dir(mut self, dir: Option<std::path::PathBuf>) -> Self {
        self.ascent_cache_dir = dir;
        self
    }

    /// Файл кэша подъёма для космодрома. Ключ — всё, от чего зависит результат: ракета, груз,
    /// широта и высота космодрома, опорная орбита, GM Земли и исходный код модели подъёма
    /// (изменение `src/ascent.rs` делает старый кэш недействительным).
    fn ascent_cache_file(&self, port: &Spaceport, mu: f64) -> Option<(std::path::PathBuf, String)> {
        let dir = self.ascent_cache_dir.as_ref()?;
        let setup = self.ascent.get(&port.id)?;
        let key = serde_json::json!({
            "model": fnv1a(crate::ascent::MODEL_SOURCE.as_bytes()),
            "vehicle": setup.vehicle,
            "upper": setup.upper,
            "payload_kg": setup.payload_kg,
            "lat_deg": port.lat_deg,
            "alt_km": port.altitude_m.unwrap_or(0.0) / 1000.0,
            "target_radius_km": EARTH_RADIUS_KM + self.config.parking_altitude_km,
            "mu_km3_s2": mu,
            "azimuth_step_deg": ASCENT_AZIMUTH_STEP_DEG,
        })
        .to_string();
        Some((dir.join(format!("ascent-{:016x}.json", fnv1a(key.as_bytes()))), key))
    }

    /// Загружает сохранённые узлы подъёма космодрома (если есть и ключ совпадает).
    fn load_ascent_cache(&self, port: &Spaceport, mu: f64) -> usize {
        let Some((path, key)) = self.ascent_cache_file(port, mu) else { return 0 };
        let Ok(text) = std::fs::read_to_string(&path) else { return 0 };
        let Ok(file) = serde_json::from_str::<AscentCacheFile>(&text) else { return 0 };
        if file.key != key {
            return 0;
        }
        let mut cache = self.ascent_cache.lock().unwrap_or_else(|e| e.into_inner());
        for (k, v) in &file.nodes {
            let node: AscentNode = v.ok_or(crate::ascent::REASON_ASCENT_ORBIT);
            let _ = cache.entry((port.id.clone(), *k)).or_default().set(Ok(node));
        }
        file.nodes.len()
    }

    /// Сохраняет узлы подъёма космодрома, если появились новые. Ошибка записи не прерывает расчёт.
    fn save_ascent_cache(&self, port: &Spaceport, mu: f64, loaded: usize) {
        let Some((path, key)) = self.ascent_cache_file(port, mu) else { return };
        let nodes: std::collections::BTreeMap<i64, Option<(f64, f64)>> = self
            .ascent_cache
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .iter()
            .filter(|((id, _), _)| id == &port.id)
            .filter_map(|((_, k), slot)| match slot.get() {
                Some(Ok(node)) => Some((*k, node.ok())),
                _ => None,
            })
            .collect();
        if nodes.len() <= loaded {
            return;
        }
        let file = AscentCacheFile { key, nodes };
        let write = || -> std::io::Result<()> {
            std::fs::create_dir_all(path.parent().unwrap_or(std::path::Path::new(".")))?;
            // Сначала во временный файл, потом переименование: параллельные запуски не увидят половину файла.
            let tmp = path.with_extension(format!("tmp{}", std::process::id()));
            std::fs::write(&tmp, serde_json::to_string(&file).map_err(std::io::Error::other)?)?;
            std::fs::rename(&tmp, &path)
        };
        if let Err(e) = write() {
            tracing::warn!("could not save the ascent cache {}: {e}", path.display());
        }
    }

    /// Задаёт ракеты для моделирования подъёма.
    pub fn with_ascent(mut self, ascent: std::collections::HashMap<String, AscentSetup>) -> Self {
        self.ascent = ascent;
        self
    }

    fn ascent_input<'a>(&self, setup: &'a AscentSetup, port: &Spaceport, azimuth_deg: f64, mu: f64) -> crate::ascent::AscentInput<'a> {
        crate::ascent::AscentInput {
            vehicle: &setup.vehicle,
            upper: setup.upper,
            payload_kg: setup.payload_kg,
            lat_deg: port.lat_deg,
            alt_km: port.altitude_m.unwrap_or(0.0) / 1000.0,
            azimuth_deg,
            target_radius_km: EARTH_RADIUS_KM + self.config.parking_altitude_km,
            mu_km3_s2: mu,
        }
    }

    /// Длительность подъёма и угол выхода на орбиту для азимута (по сетке с интерполяцией).
    fn ascent_lookup(&self, port: &Spaceport, azimuth_deg: f64, mu: f64) -> Result<Option<AscentNode>> {
        let Some(setup) = self.ascent.get(&port.id) else { return Ok(None) };
        let node = |k: i64| -> Result<AscentNode> {
            let slot = self
                .ascent_cache
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .entry((port.id.clone(), k))
                .or_default()
                .clone();
            let v = slot.get_or_init(|| {
                let az = k as f64 * ASCENT_AZIMUTH_STEP_DEG;
                crate::ascent::simulate(&self.ascent_input(setup, port, az, mu))
                    .map(|r| r.map(|r| (r.duration_s, r.downrange_angle_rad)))
                    .map_err(|e| e.to_string())
            });
            v.clone().map_err(Error::Config)
        };
        let x = azimuth_deg / ASCENT_AZIMUTH_STEP_DEG;
        let (k0, k1) = (x.floor() as i64, x.floor() as i64 + 1);
        let w = x - k0 as f64;
        Ok(Some(match (node(k0)?, node(k1)?) {
            (Ok(a), Ok(b)) => Ok((a.0 + (b.0 - a.0) * w, a.1 + (b.1 - a.1) * w)),
            (Ok(a), Err(_)) | (Err(_), Ok(a)) => Ok(a),
            (Err(e), Err(_)) => Err(e),
        }))
    }
}

/// Направление на космодром в J2000 с шагом в 1 минуту.
struct SiteTrack {
    start: Epoch,
    step: Duration,
    dirs: Vec<Vector3>,
    /// Максимальный модуль склонения космодрома, радианы.
    max_abs_dec: f64,
    frame_name: &'static str,
}

impl SiteTrack {
    fn build(almanac: &Almanac, port: &Spaceport, from: Epoch, to: Epoch) -> Result<Self> {
        let ecef = geodetic_to_ecef(port.lat_deg, port.lon_deg, port.altitude_m.unwrap_or(0.0) / 1000.0);
        let step = 1 * Unit::Minute;
        let n = ((to - from).to_seconds() / step.to_seconds()).ceil() as usize + 2;
        let (frame, frame_name) = earth_fixed_frame(almanac, from, from + step * ((n - 1) as f64))?;
        let dirs: Vec<Vector3> = (0..n)
            .into_par_iter()
            .map(|i| {
                let t = from + step * (i as f64);
                let dcm = almanac
                    .rotate(frame, EARTH_J2000, t)
                    .map_err(|e| Error::Almanac(format!("Earth rotation at {t}: {e}")))?;
                Ok((dcm * ecef).normalize())
            })
            .collect::<Result<_>>()?;
        let max_abs_dec = dirs.iter().map(|d| d.z.abs().asin()).fold(0.0, f64::max);
        Ok(Self { start: from, step, dirs, max_abs_dec, frame_name })
    }

    fn epoch(&self, i: usize) -> Epoch {
        self.start + self.step * (i as f64)
    }

    fn index_range(&self, from: Epoch, to: Epoch) -> (usize, usize) {
        let s = self.step.to_seconds();
        let lo = ((from - self.start).to_seconds() / s).floor().max(0.0) as usize;
        let hi = (((to - self.start).to_seconds() / s).ceil().max(0.0) as usize).min(self.dirs.len() - 1);
        (lo, hi)
    }
}

/// Земная система на весь промежуток `from..=to`: точная ITRF93, если в ядрах есть ориентация
/// Земли на обоих концах (файл точной ориентации покрывает сплошной интервал), иначе IAU_EARTH —
/// одна система на весь промежуток, без скачка на границе данных.
fn earth_fixed_frame(almanac: &Almanac, from: Epoch, to: Epoch) -> Result<(Frame, &'static str)> {
    let precise = |t| almanac.rotate(EARTH_ITRF93, EARTH_J2000, t).is_ok();
    if precise(from) && precise(to) {
        return Ok((EARTH_ITRF93, "ITRF93"));
    }
    tracing::warn!("high-precision Earth orientation does not cover {from} .. {to}, using IAU_EARTH");
    for t in [from, to] {
        almanac
            .rotate(IAU_EARTH_FRAME, EARTH_J2000, t)
            .map_err(|e| Error::Almanac(format!("no Earth orientation for {t}: {e}")))?;
    }
    Ok((IAU_EARTH_FRAME, "IAU_EARTH"))
}

/// Геодезические координаты WGS84 -> вектор в земной системе, км.
/// Геодезические координаты (WGS84) → положение в связанной с Землёй системе, км.
pub fn geodetic_to_ecef(lat_deg: f64, lon_deg: f64, alt_km: f64) -> Vector3 {
    let (lat, lon) = (lat_deg.to_radians(), lon_deg.to_radians());
    let e2 = EARTH_FLATTENING * (2.0 - EARTH_FLATTENING);
    let n = EARTH_RADIUS_KM / (1.0 - e2 * lat.sin().powi(2)).sqrt();
    Vector3::new(
        (n + alt_km) * lat.cos() * lon.cos(),
        (n + alt_km) * lat.cos() * lon.sin(),
        (n * (1.0 - e2) + alt_km) * lat.sin(),
    )
}

/// Нормали плоскостей с наклонением `incl`, проходящих через вектор `r` (0, 1 или 2 решения).
fn planes_through(r: &Vector3, incl: f64) -> Vec<Vector3> {
    let rho = r.x.hypot(r.y);
    if rho < 1e-9 || incl.sin() < 1e-9 {
        return Vec::new();
    }
    let phi = r.y.atan2(r.x);
    let k = -(r.z / rho) * incl.cos() / incl.sin();
    if k.abs() > 1.0 {
        return Vec::new();
    }
    let a = k.asin();
    [phi + a, phi + std::f64::consts::PI - a]
        .iter()
        .map(|&raan| Vector3::new(incl.sin() * raan.sin(), -incl.sin() * raan.cos(), incl.cos()))
        .collect()
}

/// Азимут направления `dir` в точке с единичным вектором `up` (0 = север, 90 = восток).
fn azimuth_deg(up: &Vector3, dir: &Vector3) -> f64 {
    let z = Vector3::new(0.0, 0.0, 1.0);
    let east = z.cross(up).normalize();
    let north = up.cross(&east);
    dir.dot(&east).atan2(dir.dot(&north)).to_degrees().rem_euclid(360.0)
}

/// Все данные одного варианта (внутреннее представление).
#[derive(Debug, Clone, Copy)]
struct Candidate {
    launch: Epoch,
    tli: Epoch,
    arrival: Epoch,
    incl: f64,
    h: Vector3,
    site_dir: Vector3,
    azimuth_deg: f64,
    ascending: bool,
    coast_angle: f64,
    r1: Vector3,
    v_park: Vector3,
    v1: Vector3,
    v2: Vector3,
    r2: Vector3,
    v_moon: Vector3,
    angle: f64,
    dv_tli: f64,
    v_inf: f64,
    dv_loi: f64,
    moon_dec: f64,
    required_min_deg: f64,
    landing: Option<LandingInfo>,
}

impl Candidate {
    fn total(&self) -> f64 {
        self.dv_tli + self.dv_loi
    }

    fn tof_s(&self) -> f64 {
        (self.arrival - self.tli).to_seconds()
    }

    /// Конец миссии: посадка (если считается) или прилёт к Луне.
    fn end(&self) -> Epoch {
        self.landing.map_or(self.arrival, |l| l.epoch)
    }

    /// От старта до конца миссии, с.
    fn mission_s(&self) -> f64 {
        (self.end() - self.launch).to_seconds()
    }

    /// `true`, если `self` лучше `other` по критерию.
    fn better_than(&self, other: &Candidate, objective: Objective) -> bool {
        match objective {
            Objective::MinDeltaV => self.total() < other.total(),
            Objective::Fastest => {
                // Время перелёта сравнивается с точностью до часа: среди почти одинаково
                // быстрых вариантов берётся самый экономичный, а не случайный с огромным Δv.
                // С точкой посадки «быстрее» значит раньше сесть (от старта до посадки).
                let (a, b) = ((self.mission_s() / 3600.0).floor(), (other.mission_s() / 3600.0).floor());
                a < b || (a == b && self.total() < other.total())
            }
        }
    }
}

#[derive(Default)]
struct Counters {
    arrival_epochs: AtomicU64,
    inclinations_tried: AtomicU64,
    planes_tried: AtomicU64,
    plane_crossings: AtomicU64,
    crossings_rejected_azimuth: AtomicU64,
    crossings_rejected_ascent: AtomicU64,
    candidates_rejected_budget: AtomicU64,
    candidates_rejected_tli_limit: AtomicU64,
    arrivals_rejected_landing: AtomicU64,
    lambert_solves: AtomicU64,
    lambert_failures: AtomicU64,
    candidates_evaluated: AtomicU64,
}

impl Counters {
    fn add(c: &AtomicU64, n: u64) {
        c.fetch_add(n, Ordering::Relaxed);
    }

    fn snapshot(&self) -> SearchStats {
        let g = |c: &AtomicU64| c.load(Ordering::Relaxed);
        SearchStats {
            arrival_epochs: g(&self.arrival_epochs),
            inclinations_tried: g(&self.inclinations_tried),
            planes_tried: g(&self.planes_tried),
            plane_crossings: g(&self.plane_crossings),
            crossings_rejected_azimuth: g(&self.crossings_rejected_azimuth),
            crossings_rejected_ascent: g(&self.crossings_rejected_ascent),
            candidates_rejected_budget: g(&self.candidates_rejected_budget),
            candidates_rejected_tli_limit: g(&self.candidates_rejected_tli_limit),
            arrivals_rejected_landing: g(&self.arrivals_rejected_landing),
            lambert_solves: g(&self.lambert_solves),
            lambert_failures: g(&self.lambert_failures),
            candidates_evaluated: g(&self.candidates_evaluated),
        }
    }
}

struct Search<'a> {
    almanac: &'a Almanac,
    port: &'a Spaceport,
    cfg: &'a DirectTransferConfig,
    transfer: &'a DirectTransfer,
    site: &'a SiteTrack,
    mu_e: f64,
    mu_m: f64,
    counters: Counters,
    budget: Option<&'a DeltaVBudget>,
    landing: Option<&'a LandingTarget>,
    /// Допустимый интервал момента посадки (сутки посадки для ArrivalDay).
    landing_day: Option<(Epoch, Epoch)>,
}

impl Search<'_> {
    fn objective(&self) -> Objective {
        self.budget.map_or(Objective::MinDeltaV, |b| b.objective)
    }

    /// Лучший вариант для прилёта в момент `t_a`, либо причина, почему вариантов нет.
    fn at_arrival(
        &self,
        t_a: Epoch,
        coast_step: Duration,
        (launch_lo, launch_hi): (Epoch, Epoch),
    ) -> Result<std::result::Result<Candidate, &'static str>> {
        let cfg = self.cfg;
        Counters::add(&self.counters.arrival_epochs, 1);
        let moon = self
            .almanac
            .translate(MOON_J2000, EARTH_J2000, t_a, None)
            .map_err(|e| Error::Almanac(format!("Moon at {t_a}: {e}")))?;
        let (r2, v_moon) = (moon.radius_km, moon.velocity_km_s);
        let moon_dec = (r2.z / r2.norm()).asin();

        let tof_min = cfg.tof_min_days * Unit::Day;
        let tof_max = cfg.tof_max_days * Unit::Day;
        let coast_min = cfg.coast_min_min * Unit::Minute;
        let coast_max = cfg.coast_max_min * Unit::Minute;
        let w_lo = (t_a - tof_max - coast_max).max(launch_lo);
        let w_hi = (t_a - tof_min - coast_min).min(launch_hi);
        if w_hi <= w_lo {
            return Ok(Err(REASON_WINDOW));
        }

        let required_min_deg =
            self.site.max_abs_dec.max(moon_dec.abs()).to_degrees() + cfg.inclination_margin_deg;
        let inclinations =
            self.port.inclination_candidates(required_min_deg, INCLINATION_STEP_DEG, MAX_INCLINATION_DEG);
        if inclinations.is_empty() {
            return Ok(Err("no inclination allowed by the spaceport reaches the Moon's declination"));
        }

        let r_p = EARTH_RADIUS_KM + cfg.parking_altitude_km;
        let n_park = (self.mu_e / r_p.powi(3)).sqrt();
        // Подъём: через ascent_duration выход на орбиту на ascent_downrange от космодрома.
        // Положение на опорной орбите в момент t после старта — угол ascent_phase + n·t от космодрома.
        let kinematic_ascent = (cfg.ascent_duration_min * 60.0, cfg.ascent_downrange_km / EARTH_RADIUS_KM);
        let v_circ = (self.mu_e / r_p).sqrt();
        let r_llo = MOON_RADIUS_KM + cfg.lunar_orbit_altitude_km;
        let v_llo = (self.mu_m / r_llo).sqrt();
        let (i_lo, i_hi) = self.site.index_range(w_lo, w_hi);
        let mut rejected_by_azimuth = false;
        let mut rejected_by_ascent = false;
        let mut rejected_by_budget = false;
        let mut rejected_by_tli = false;
        let objective = self.objective();

        // Наклонения по возрастанию: берём наименьшее, при котором есть хоть один вариант
        // (меньшее наклонение — пуск ближе к востоку, больше выигрыш от вращения Земли).
        for &incl_deg in &inclinations {
            Counters::add(&self.counters.inclinations_tried, 1);
            let incl = incl_deg.to_radians();
            let mut best: Option<Candidate> = None;
            for h in planes_through(&r2, incl) {
                Counters::add(&self.counters.planes_tried, 1);
                for i in i_lo..i_hi {
                    let (g0, g1) = (h.dot(&self.site.dirs[i]), h.dot(&self.site.dirs[i + 1]));
                    if g0 == 0.0 || g0.signum() == g1.signum() {
                        continue;
                    }
                    // Момент и место пересечения плоскости (линейная интерполяция внутри минуты).
                    let f = g0 / (g0 - g1);
                    let t_l = self.site.epoch(i) + self.site.step * f;
                    if t_l < w_lo || t_l >= w_hi {
                        continue;
                    }
                    Counters::add(&self.counters.plane_crossings, 1);
                    let s = (self.site.dirs[i] * (1.0 - f) + self.site.dirs[i + 1] * f).normalize();
                    let l_hat = (s - h * h.dot(&s)).normalize();
                    let along = h.cross(&l_hat);
                    let az = azimuth_deg(&s, &along);
                    if !self.port.azimuth_allowed(az) {
                        Counters::add(&self.counters.crossings_rejected_azimuth, 1);
                        rejected_by_azimuth = true;
                        continue;
                    }

                    // Подъём: моделирование ракеты (если задана) или упрощённая модель.
                    let (asc_t, asc_angle) = match self.transfer.ascent_lookup(self.port, az, self.mu_e)? {
                        None => kinematic_ascent,
                        Some(Ok(v)) => v,
                        Some(Err(_)) => {
                            Counters::add(&self.counters.crossings_rejected_ascent, 1);
                            rejected_by_ascent = true;
                            continue;
                        }
                    };
                    let ascent_phase = asc_angle - n_park * asc_t;
                    let k_min = (coast_min.to_seconds() / coast_step.to_seconds()).ceil() as i64;
                    let k_max = (coast_max.to_seconds() / coast_step.to_seconds()).floor() as i64;
                    // Внутренний цикл — в секундах (f64) и с локальными счётчиками: арифметика дат
                    // и общие атомарные счётчики здесь были дороже самой задачи Ламберта.
                    let step_s = coast_step.to_seconds();
                    let launch_to_arrival_s = (t_a - t_l).to_seconds();
                    let (tof_min_s, tof_max_s) = (tof_min.to_seconds(), tof_max.to_seconds());
                    let (mut n_solves, mut n_fail, mut n_tli, mut n_budget) = (0u64, 0u64, 0u64, 0u64);
                    for k in k_min..=k_max {
                        let coast_s = step_s * k as f64;
                        if coast_s < asc_t + 60.0 {
                            continue;
                        }
                        let tof_s = launch_to_arrival_s - coast_s;
                        if tof_s < tof_min_s || tof_s > tof_max_s {
                            continue;
                        }
                        let ang = ascent_phase + n_park * coast_s;
                        let r1_hat = l_hat * ang.cos() + along * ang.sin();
                        let r1 = r1_hat * r_p;
                        n_solves += 1;
                        let Ok(sol) = lambert::solve(&r1, &r2, tof_s, self.mu_e, &h) else {
                            n_fail += 1;
                            continue;
                        };
                        let v_park = h.cross(&r1_hat) * v_circ;
                        let dv_tli = (sol.v1 - v_park).norm();
                        if dv_tli > cfg.max_tli_delta_v_km_s {
                            n_tli += 1;
                            continue;
                        }
                        let v_inf = (sol.v2 - v_moon).norm();
                        let dv_loi = (v_inf * v_inf + 2.0 * self.mu_m / r_llo).sqrt() - v_llo;
                        if let Some(b) = self.budget
                            && !b.fits(dv_tli, dv_loi)
                        {
                            n_budget += 1;
                            continue;
                        }
                        let t_d = t_l + coast_step * (k as f64);
                        let cand = Candidate {
                            launch: t_l,
                            tli: t_d,
                            arrival: t_a,
                            incl,
                            h,
                            site_dir: s,
                            azimuth_deg: az,
                            ascending: along.z > 0.0,
                            coast_angle: ang,
                            r1,
                            v_park,
                            v1: sol.v1,
                            v2: sol.v2,
                            r2,
                            v_moon,
                            angle: sol.transfer_angle_rad,
                            dv_tli,
                            v_inf,
                            dv_loi,
                            moon_dec,
                            required_min_deg,
                            landing: None,
                        };
                        if best.is_none_or(|b| cand.better_than(&b, objective)) {
                            best = Some(cand);
                        }
                    }
                    Counters::add(&self.counters.lambert_solves, n_solves);
                    Counters::add(&self.counters.lambert_failures, n_fail);
                    Counters::add(&self.counters.candidates_evaluated, n_solves - n_fail);
                    Counters::add(&self.counters.candidates_rejected_tli_limit, n_tli);
                    Counters::add(&self.counters.candidates_rejected_budget, n_budget);
                    rejected_by_tli |= n_tli > 0;
                    rejected_by_budget |= n_budget > 0;
                }
            }
            if let Some(mut b) = best {
                if let Some(target) = self.landing {
                    // Δv, израсходованный аппаратом до спуска, — для его массы в начале торможения.
                    let prior_m_s: f64 = self.budget.map_or(0.0, |bg| {
                        bg.maneuvers
                            .iter()
                            .map(|m| match m {
                                Maneuver::Tli => b.dv_tli * 1000.0,
                                Maneuver::Loi => b.dv_loi * 1000.0,
                                Maneuver::Descent => 0.0,
                            })
                            .sum()
                    });
                    match find_landing(
                        self.almanac,
                        target,
                        b.arrival,
                        b.v2 - b.v_moon,
                        self.mu_m,
                        r_llo,
                        prior_m_s,
                        self.landing_day,
                    )? {
                        Ok(info) => b.landing = Some(info),
                        Err(reason) => {
                            Counters::add(&self.counters.arrivals_rejected_landing, 1);
                            return Ok(Err(reason));
                        }
                    }
                }
                return Ok(Ok(b));
            }
        }
        Ok(Err(if rejected_by_budget {
            REASON_BUDGET
        } else if rejected_by_tli {
            "every transfer arriving at this time needs more TLI delta-v than max_tli_delta_v_km_s"
        } else if rejected_by_ascent {
            crate::ascent::REASON_ASCENT_ORBIT
        } else if rejected_by_azimuth {
            "all launch opportunities fall outside the spaceport's azimuth sector"
        } else {
            "no launch opportunity within the coast and flight-time limits"
        }))
    }

    fn scan(
        &self,
        ta_lo: Epoch,
        ta_hi: Epoch,
        ta_step: Duration,
        coast_step: Duration,
        launch: (Epoch, Epoch),
    ) -> Result<(Option<Candidate>, Vec<ScanPoint>)> {
        let n = ((ta_hi - ta_lo).to_seconds() / ta_step.to_seconds()).floor() as usize;
        let results: Vec<(Epoch, std::result::Result<Candidate, &'static str>)> = (0..=n)
            .into_par_iter()
            .map(|i| {
                let t_a = ta_lo + ta_step * (i as f64);
                Ok((t_a, self.at_arrival(t_a, coast_step, launch)?))
            })
            .collect::<Result<_>>()?;
        let points = results
            .iter()
            .map(|(t_a, r)| match r {
                Ok(c) => ScanPoint {
                    arrival: *t_a,
                    feasible: true,
                    reason: None,
                    launch: Some(c.launch),
                    coast_min: Some((c.tli - c.launch).to_seconds() / 60.0),
                    time_of_flight_h: Some((c.arrival - c.tli).to_seconds() / 3600.0),
                    inclination_deg: Some(c.incl.to_degrees()),
                    azimuth_deg: Some(c.azimuth_deg),
                    delta_v_tli_km_s: Some(c.dv_tli),
                    delta_v_loi_km_s: Some(c.dv_loi),
                    delta_v_total_km_s: Some(c.total()),
                    landing: c.landing.map(|l| l.epoch),
                },
                Err(reason) => ScanPoint {
                    arrival: *t_a,
                    feasible: false,
                    reason: Some(reason),
                    launch: None,
                    coast_min: None,
                    time_of_flight_h: None,
                    inclination_deg: None,
                    azimuth_deg: None,
                    delta_v_tli_km_s: None,
                    delta_v_loi_km_s: None,
                    delta_v_total_km_s: None,
                    landing: None,
                },
            })
            .collect();
        let objective = self.objective();
        let best = results.into_iter().filter_map(|(_, r)| r.ok()).fold(None, |acc: Option<Candidate>, c| {
            match acc {
                Some(a) if !c.better_than(&a, objective) => Some(a),
                _ => Some(c),
            }
        });
        Ok((best, points))
    }

    fn solution(&self, c: &Candidate, name: &'static str) -> Result<TransferSolution> {
        let cfg = self.cfg;
        let mu = self.mu_e;
        // Подъём для найденного азимута: точное моделирование ракеты (если задана).
        let setup = self.transfer.ascent.get(&self.port.id);
        let (sim, sim_failure) = match setup {
            Some(st) => match crate::ascent::simulate(&self.transfer.ascent_input(st, self.port, c.azimuth_deg, mu))? {
                Ok(r) => (Some(r), None),
                Err(reason) => (None, Some(reason)),
            },
            None => (None, None),
        };
        let to_j2000 = sim.as_ref().map(|r| r.to_j2000(&c.site_dir, &c.h));
        let r_p = EARTH_RADIUS_KM + cfg.parking_altitude_km;
        let r1n = c.r1.norm();
        // Элементы орбиты перелёта по состоянию сразу после разгона.
        let energy = c.v1.norm_squared() / 2.0 - mu / r1n;
        let h_vec = c.r1.cross(&c.v1);
        let e_vec = c.v1.cross(&h_vec) / mu - c.r1 / r1n;
        let ecc = e_vec.norm();
        let sma = -mu / (2.0 * energy);
        let p = h_vec.norm_squared() / mu;
        let perigee = p / (1.0 + ecc);
        let apogee = (ecc < 1.0).then(|| p / (1.0 - ecc));
        let fpa = (c.r1.dot(&c.v1) / (r1n * c.v1.norm())).clamp(-1.0, 1.0).asin();
        let r_llo = MOON_RADIUS_KM + cfg.lunar_orbit_altitude_km;
        let v_llo = (self.mu_m / r_llo).sqrt();
        let v_peri = (c.v_inf * c.v_inf + 2.0 * self.mu_m / r_llo).sqrt();
        let port_lat = self.port.lat_deg.to_radians();
        let available = self.port.inclination_candidates(
            c.required_min_deg,
            INCLINATION_STEP_DEG,
            MAX_INCLINATION_DEG,
        );

        Ok(TransferSolution {
            trajectory: name,
            spaceport_id: self.port.id.clone(),
            launch: LaunchInfo {
                epoch: c.launch,
                azimuth_deg: c.azimuth_deg,
                ascending_pass: c.ascending,
                earth_rotation_gain_km_s: EARTH_EQUATOR_SPEED_KM_S
                    * port_lat.cos()
                    * c.azimuth_deg.to_radians().sin(),
                site_direction_j2000: c.site_dir,
                earth_frame: self.site.frame_name,
            },
            ascent: {
                let l = (c.site_dir - c.h * c.h.dot(&c.site_dir)).normalize();
                let along = c.h.cross(&l);
                let pad = EARTH_RADIUS_KM + self.port.altitude_m.unwrap_or(0.0) / 1000.0;
                let rot_rate = EARTH_EQUATOR_SPEED_KM_S * port_lat.cos() * c.azimuth_deg.to_radians().sin() / pad;
                match (&sim, to_j2000) {
                    (Some(r), Some(m)) => AscentInfo {
                        model: "simulation",
                        rocket: setup.map(|st| st.rocket_name.clone()),
                        model_note: None,
                        duration_s: r.duration_s,
                        downrange_km: r.downrange_km,
                        downrange_angle_deg: r.downrange_angle_rad.to_degrees(),
                        insertion_epoch: c.launch + r.duration_s * Unit::Second,
                        launch_direction_in_plane: l,
                        insertion_position_km: m * r.insertion_r,
                        insertion_velocity_km_s: m * r.insertion_v,
                        pad_radius_km: r.site_r.norm(),
                        earth_rotation_rate_along_rad_s: rot_rate,
                        simulation: Some(r.clone()),
                        to_j2000: Some(m),
                    },
                    _ => {
                        let theta = cfg.ascent_downrange_km / EARTH_RADIUS_KM;
                        let u = l * theta.cos() + along * theta.sin();
                        AscentInfo {
                            model: "kinematic",
                            rocket: setup.map(|st| st.rocket_name.clone()),
                            model_note: Some(match sim_failure {
                                Some(reason) => format!("ascent simulation failed ({reason}); simplified ascent profile used"),
                                None => "no rocket with stage data for this spaceport; simplified ascent profile (duration and downrange from transfer.toml)".to_string(),
                            }),
                            duration_s: cfg.ascent_duration_min * 60.0,
                            downrange_km: cfg.ascent_downrange_km,
                            downrange_angle_deg: theta.to_degrees(),
                            insertion_epoch: c.launch + cfg.ascent_duration_min * Unit::Minute,
                            launch_direction_in_plane: l,
                            insertion_position_km: u * r_p,
                            insertion_velocity_km_s: c.h.cross(&u) * (mu / r_p).sqrt(),
                            pad_radius_km: pad,
                            earth_rotation_rate_along_rad_s: rot_rate,
                            simulation: None,
                            to_j2000: None,
                        }
                    }
                }
            },
            parking_orbit: ParkingOrbit {
                altitude_km: cfg.parking_altitude_km,
                radius_km: r_p,
                circular_speed_km_s: (mu / r_p).sqrt(),
                period_s: std::f64::consts::TAU * (r_p.powi(3) / mu).sqrt(),
                inclination_deg: c.incl.to_degrees(),
                raan_deg: c.h.x.atan2(-c.h.y).to_degrees().rem_euclid(360.0),
                plane_normal: c.h,
            },
            departure: DepartureInfo {
                epoch: c.tli,
                coast_s: c.tli - c.launch,
                coast_angle_deg: c.coast_angle.to_degrees(),
                position_km: c.r1,
                velocity_parking_km_s: c.v_park,
                velocity_transfer_km_s: c.v1,
                delta_v_km_s: c.v1 - c.v_park,
                delta_v_magnitude_km_s: c.dv_tli,
                flight_path_angle_deg: fpa.to_degrees(),
                c3_km2_s2: 2.0 * energy,
                transfer_semi_major_axis_km: sma,
                transfer_eccentricity: ecc,
                transfer_perigee_km: perigee,
                transfer_apogee_km: apogee,
                transfer_angle_deg: c.angle.to_degrees(),
                time_of_flight_s: c.arrival - c.tli,
            },
            arrival: ArrivalInfo {
                epoch: c.arrival,
                moon_position_km: c.r2,
                moon_velocity_km_s: c.v_moon,
                moon_distance_km: c.r2.norm(),
                moon_declination_deg: c.moon_dec.to_degrees(),
                velocity_arrival_km_s: c.v2,
                v_infinity_km_s: c.v2 - c.v_moon,
                v_infinity_magnitude_km_s: c.v_inf,
                lunar_orbit_altitude_km: cfg.lunar_orbit_altitude_km,
                lunar_orbit_radius_km: r_llo,
                lunar_orbit_circular_speed_km_s: v_llo,
                periselene_speed_km_s: v_peri,
                delta_v_loi_km_s: c.dv_loi,
            },
            inclination: InclinationChoice {
                site_max_declination_deg: self.site.max_abs_dec.to_degrees(),
                moon_declination_deg: c.moon_dec.to_degrees(),
                margin_deg: cfg.inclination_margin_deg,
                required_min_deg: c.required_min_deg,
                candidates_deg: available,
                chosen_deg: c.incl.to_degrees(),
            },
            delta_v_total_km_s: c.total(),
            mu_earth_km3_s2: self.mu_e,
            mu_moon_km3_s2: self.mu_m,
            landing: c.landing,
            refined: None,
            refine_failure: None,
            refine_skipped: None,
            budget: self.budget.map(|b| {
                let tof_min_s = cfg.tof_min_days * 86_400.0;
                b.check(c.dv_tli, c.dv_loi, None, c.tof_s() < tof_min_s + 3600.0)
            }),
        })
    }
}

/// Самая содержательная причина отсутствия вариантов: бюджет аппарата важнее остальных,
/// «окно не перекрывается» (края диапазона перебора) — наименее содержательна.
fn main_reason(scan: &[ScanPoint]) -> Option<&'static str> {
    let reasons: Vec<&'static str> = scan.iter().filter_map(|p| p.reason).collect();
    if reasons.contains(&REASON_BUDGET) {
        return Some(REASON_BUDGET);
    }
    let mut counts: Vec<(&'static str, usize)> = Vec::new();
    for r in reasons.iter().copied().filter(|r| *r != REASON_WINDOW) {
        match counts.iter_mut().find(|(x, _)| *x == r) {
            Some((_, n)) => *n += 1,
            None => counts.push((r, 1)),
        }
    }
    counts.into_iter().max_by_key(|(_, n)| *n).map(|(r, _)| r).or_else(|| reasons.first().copied())
}

impl Trajectory for DirectTransfer {
    fn name(&self) -> &'static str {
        "direct_transfer"
    }

    fn solve(
        &self,
        almanac: &Almanac,
        spaceport: &Spaceport,
        constraint: DateConstraint,
        budget: Option<&DeltaVBudget>,
        landing: Option<&LandingTarget>,
    ) -> Result<TransferReport> {
        let cfg = &self.config;
        cfg.validate()?;
        let w = self.windows(constraint, landing);
        let (launch_lo, launch_hi, ta_lo, ta_hi) = (w.launch_lo, w.launch_hi, w.ta_lo, w.ta_hi);
        let site = SiteTrack::build(almanac, spaceport, launch_lo, launch_hi)?;
        let mu_e = gm(almanac, EARTH_J2000, MU_EARTH_DEFAULT);
        let ascent_loaded = self.load_ascent_cache(spaceport, mu_e);
        let search = Search {
            almanac,
            port: spaceport,
            cfg,
            transfer: self,
            site: &site,
            mu_e,
            mu_m: gm(almanac, MOON_J2000, MU_MOON_DEFAULT),
            counters: Counters::default(),
            budget,
            landing,
            landing_day: w.landing_day,
        };

        // Грубый перебор момента прилёта, затем уточнение вокруг лучшего.
        let coarse_step = cfg.search_step_min * Unit::Minute;
        let (coarse, coarse_scan) = search.scan(ta_lo, ta_hi, coarse_step, 2 * Unit::Minute, (launch_lo, launch_hi))?;
        let (best, fine_scan) = match coarse {
            Some(c) => {
                let fine_lo = (c.arrival - coarse_step).max(ta_lo);
                let fine_hi = (c.arrival + coarse_step).min(ta_hi);
                // Окно пуска не сужается: поле Δv пологое, и близкие по Δv варианты бывают в
                // других окнах пуска (сужение меняло ответ при разнице ~1 м/с).
                let (fine, scan) =
                    search.scan(fine_lo, fine_hi, 1 * Unit::Minute, 30 * Unit::Second, (launch_lo, launch_hi))?;
                (Some(fine.unwrap_or(c)), scan)
            }
            None => (None, Vec::new()),
        };
        self.save_ascent_cache(spaceport, mu_e, ascent_loaded);

        Ok(TransferReport {
            trajectory: self.name(),
            spaceport: spaceport.clone(),
            constraint,
            config: cfg.clone(),
            launch_window_start: launch_lo,
            launch_window_end: launch_hi,
            arrival_window_start: ta_lo,
            arrival_window_end: ta_hi,
            earth_frame: site.frame_name,
            site_max_declination_deg: site.max_abs_dec.to_degrees(),
            budget: budget.cloned(),
            landing_site: landing.map(|l| l.site.site),
            landing_window: landing.map(|l| l.window.clone()),
            no_solution_reason: if best.is_some() { None } else { main_reason(&coarse_scan) },
            solution: best.map(|c| search.solution(&c, self.name())).transpose()?,
            coarse_scan,
            fine_scan,
            stats: search.counters.snapshot(),
            refine_attempts: Vec::new(),
        })
    }

    fn solve_at(
        &self,
        almanac: &Almanac,
        spaceport: &Spaceport,
        constraint: DateConstraint,
        budget: Option<&DeltaVBudget>,
        landing: Option<&LandingTarget>,
        arrival: Epoch,
    ) -> Result<Option<TransferSolution>> {
        let cfg = &self.config;
        cfg.validate()?;
        let w = self.windows(constraint, landing);
        let site = SiteTrack::build(almanac, spaceport, w.launch_lo, w.launch_hi)?;
        let search = Search {
            almanac,
            port: spaceport,
            cfg,
            transfer: self,
            site: &site,
            mu_e: gm(almanac, EARTH_J2000, MU_EARTH_DEFAULT),
            mu_m: gm(almanac, MOON_J2000, MU_MOON_DEFAULT),
            counters: Counters::default(),
            budget,
            landing,
            landing_day: w.landing_day,
        };
        // Тот же шаг по времени на опорной орбите, что и в уточняющем переборе.
        match search.at_arrival(arrival, 30 * Unit::Second, (w.launch_lo, w.launch_hi))? {
            Ok(c) => Ok(Some(search.solution(&c, self.name())?)),
            Err(_) => Ok(None),
        }
    }
}

/// Гравитационный параметр тела из ядер (или значение по умолчанию), км³/с².
fn gm(almanac: &Almanac, frame: Frame, default: f64) -> f64 {
    almanac.frame_info(frame).ok().and_then(|f| f.mu_km3_s2().ok()).unwrap_or(default)
}

/// Окна поиска для ограничения по дате.
struct Windows {
    launch_lo: Epoch,
    launch_hi: Epoch,
    ta_lo: Epoch,
    ta_hi: Epoch,
    /// Допустимый интервал момента посадки (сутки посадки для ArrivalDay с точкой посадки).
    landing_day: Option<(Epoch, Epoch)>,
}

impl DirectTransfer {
    fn windows(&self, constraint: DateConstraint, landing: Option<&LandingTarget>) -> Windows {
        let cfg = &self.config;
        let day = 1 * Unit::Day;
        let tof_min = cfg.tof_min_days * Unit::Day;
        let tof_max = cfg.tof_max_days * Unit::Day;
        let coast_min = cfg.coast_min_min * Unit::Minute;
        let coast_max = cfg.coast_max_min * Unit::Minute;
        // С точкой посадки «сутки прилёта» — это сутки посадки: прилёт к Луне может быть
        // раньше на время ожидания на орбите и спуска (полвитка эллипса + торможение, < 2 ч).
        let wait = landing.map_or(Duration::ZERO, |l| l.window.max_wait_in_orbit_h * Unit::Hour + 2 * Unit::Hour);
        let (launch_lo, launch_hi, ta_lo, ta_hi) = match constraint {
            DateConstraint::LaunchDay(d) => (d, d + day, d + coast_min + tof_min, d + day + coast_max + tof_max),
            DateConstraint::ArrivalDay(d) => {
                (d - wait - tof_max - coast_max, d + day - tof_min - coast_min, d - wait, d + day)
            }
        };
        let landing_day = match (constraint, landing) {
            (DateConstraint::ArrivalDay(d), Some(_)) => Some((d, d + day)),
            _ => None,
        };
        Windows { launch_lo, launch_hi, ta_lo, ta_hi, landing_day }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ephemeris::{AlmanacSource, load_almanac_blocking};
    use crate::spaceport::SpaceportCatalog;
    use approx::assert_abs_diff_eq;

    #[test]
    fn planes_contain_the_vector_and_have_the_inclination() {
        let r = Vector3::new(300_000.0, 200_000.0, -120_000.0);
        let dec = (r.z / r.norm()).asin();
        let incl = dec.abs() + 0.1;
        let planes = planes_through(&r, incl);
        assert_eq!(planes.len(), 2);
        for h in planes {
            assert_abs_diff_eq!(h.norm(), 1.0, epsilon = 1e-12);
            assert_abs_diff_eq!(h.dot(&r) / r.norm(), 0.0, epsilon = 1e-12);
            assert_abs_diff_eq!(h.z.acos(), incl, epsilon = 1e-12);
        }
        // Наклонение меньше склонения — плоскости нет.
        assert!(planes_through(&r, dec.abs() - 0.05).is_empty());
    }

    #[test]
    fn azimuth_east_and_north() {
        let up = Vector3::new(1.0, 0.0, 0.0);
        assert_abs_diff_eq!(azimuth_deg(&up, &Vector3::new(0.0, 1.0, 0.0)), 90.0, epsilon = 1e-9);
        assert_abs_diff_eq!(azimuth_deg(&up, &Vector3::new(0.0, 0.0, 1.0)), 0.0, epsilon = 1e-9);
    }

    #[test]
    fn geodetic_conversion_on_axes() {
        let eq = geodetic_to_ecef(0.0, 90.0, 0.0);
        assert_abs_diff_eq!(eq.y, EARTH_RADIUS_KM, epsilon = 1e-6);
        let pole = geodetic_to_ecef(90.0, 0.0, 0.0);
        assert_abs_diff_eq!(pole.z, EARTH_RADIUS_KM * (1.0 - EARTH_FLATTENING), epsilon = 1e-6);
    }

    #[test]
    fn config_defaults_and_validation() {
        let c = TransferConfig::parse("", Path::new("t.toml")).unwrap();
        assert_eq!(c, TransferConfig::default());
        assert!(TransferConfig::parse("[direct_transfer]\ntof_min_days = 6.0\n", Path::new("t.toml")).is_err());
        assert!(TransferConfig::parse("[direct_transfer]\nparking_alt = 200\n", Path::new("t.toml")).is_err());
    }

    #[test]
    fn budget_arithmetic() {
        let b = DeltaVBudget::new("x", 3000.0, &MissionConfig::default()).unwrap();
        // LOI 0.8 км/с + спуск 2053 м/с = 2853; с запасом 5 % = 2995.65 <= 3000.
        assert_abs_diff_eq!(b.required_m_s(3.1, 0.8), 2853.0, epsilon = 1e-9);
        assert!(b.fits(3.1, 0.8));
        assert!(!b.fits(3.1, 1.0)); // 3053 * 1.05 = 3205.65 > 3000
        let with_tli = MissionConfig { spacecraft_maneuvers: vec![Maneuver::Tli, Maneuver::Loi], ..Default::default() };
        let b = DeltaVBudget::new("x", 5000.0, &with_tli).unwrap();
        assert_abs_diff_eq!(b.required_m_s(3.1, 0.8), 3900.0, epsilon = 1e-9);
        assert!(DeltaVBudget::new("x", 0.0, &MissionConfig::default()).is_err());
        let empty = MissionConfig { spacecraft_maneuvers: vec![], ..Default::default() };
        assert!(empty.validate().is_err());
    }

    #[test]
    fn bundled_config_file_is_valid() {
        let text = include_str!("../data/transfer.toml");
        TransferConfig::parse(text, Path::new("data/transfer.toml")).unwrap();
    }

    // ---- С настоящими ядрами из data/kernels: `cargo test -- --ignored` ----

    fn real() -> (Almanac, SpaceportCatalog) {
        let almanac = load_almanac_blocking(AlmanacSource::LocalDir("data/kernels".into())).unwrap();
        let text = std::fs::read_to_string("data/spaceports.toml").unwrap();
        let catalog = SpaceportCatalog::parse(&text, Path::new("data/spaceports.toml")).unwrap();
        (almanac, catalog)
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_earth_rotation_matches_sidereal_time() {
        // Прямое восхождение космодрома ~ GMST + долгота (с точностью до прецессии с 2000 г., < 0.5 градуса).
        let (almanac, catalog) = real();
        let ksc = catalog.get("ksc").unwrap();
        let t = Epoch::from_gregorian_utc_at_midnight(2026, 11, 14);
        let track = SiteTrack::build(&almanac, ksc, t, t + 1 * Unit::Minute).unwrap();
        let ra = track.dirs[0].y.atan2(track.dirs[0].x).to_degrees().rem_euclid(360.0);
        let d = (t - Epoch::from_gregorian_utc_hms(2000, 1, 1, 12, 0, 0)).to_unit(Unit::Day);
        let gmst = (280.460_618_37 + 360.985_647_366_29 * d).rem_euclid(360.0);
        let expected = (gmst + ksc.lon_deg).rem_euclid(360.0);
        let diff = ((ra - expected + 540.0).rem_euclid(360.0) - 180.0).abs();
        assert!(diff < 0.5, "RA {ra:.3}, expected {expected:.3}");
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_solve_at_reproduces_the_search_result() {
        let (almanac, catalog) = real();
        let ksc = catalog.get("ksc").unwrap();
        let dt = DirectTransfer::default();
        let c = DateConstraint::LaunchDay(Epoch::from_gregorian_utc_at_midnight(2026, 12, 1));
        let r = dt.solve(&almanac, ksc, c, None, None).unwrap();
        let s = r.solution.expect("solution");
        // Тот же момент прилёта даёт тот же вариант (шаг по опорной орбите тот же, 30 с).
        let again = dt.solve_at(&almanac, ksc, c, None, None, s.arrival.epoch).unwrap().expect("solution at arrival");
        assert_eq!(again.arrival.epoch, s.arrival.epoch);
        assert_eq!(again.launch.epoch, s.launch.epoch);
        assert!((again.delta_v_total_km_s - s.delta_v_total_km_s).abs() < 1e-9);
        // Другой момент прилёта из перебора — свой вариант с этим прилётом.
        let other = r.coarse_scan.iter().find(|p| p.feasible && (p.arrival - s.arrival.epoch).abs() > 2 * Unit::Hour).unwrap();
        let alt = dt.solve_at(&almanac, ksc, c, None, None, other.arrival).unwrap().expect("alternative");
        assert_eq!(alt.arrival.epoch, other.arrival);
        assert!(alt.delta_v_total_km_s <= other.delta_v_total_km_s.unwrap() + 1e-9);
    }

    #[test]
    #[ignore = "needs real kernels in data/kernels"]
    fn real_transfer_is_physical_and_consistent_both_ways() {
        let (almanac, catalog) = real();
        let ksc = catalog.get("ksc").unwrap();
        let dt = DirectTransfer::default();
        let day = Epoch::from_gregorian_utc_at_midnight(2026, 12, 1);
        let fwd = dt.solve(&almanac, ksc, DateConstraint::LaunchDay(day), None, None).unwrap();
        let s = fwd.solution.clone().expect("solution");

        // Старт в заданные сутки, цепочка времён согласована.
        assert!(s.launch.epoch >= day && s.launch.epoch < day + 1 * Unit::Day);
        assert_eq!(s.launch.epoch + s.departure.coast_s, s.departure.epoch);
        // Подъём: старт с поверхности, выход на опорную орбиту с её скоростью, высота растёт монотонно.
        let a = &s.ascent;
        let (h, rp) = (s.parking_orbit.plane_normal, s.parking_orbit.radius_km);
        let (r0, v0) = a.state_at(0.0, &h, rp);
        assert!((r0.norm() - a.pad_radius_km).abs() < 1e-9 && v0.norm() < 0.5);
        let (r1, v1) = a.state_at(a.duration_s, &h, rp);
        assert!((r1 - a.insertion_position_km).norm() < 1e-6, "{}", (r1 - a.insertion_position_km).norm());
        assert!((v1 - a.insertion_velocity_km_s).norm() < 1e-9);
        let alts: Vec<f64> = (0..=60).map(|i| a.state_at(a.duration_s * i as f64 / 60.0, &h, rp).0.norm()).collect();
        assert!(alts.windows(2).all(|w| w[1] >= w[0] - 1e-9));
        assert!(((r0.normalize().dot(&r1.normalize())).acos() - a.downrange_angle_deg.to_radians()).abs() < 1e-3);
        assert_eq!(a.insertion_epoch, s.launch.epoch + a.duration_s * Unit::Second);
        assert_eq!(s.departure.epoch + s.departure.time_of_flight_s, s.arrival.epoch);
        // Типичные значения для перелёта с опорной орбиты 185 км на окололунную 100 км.
        assert!(s.departure.delta_v_magnitude_km_s > 3.0 && s.departure.delta_v_magnitude_km_s < 3.3);
        assert!(s.arrival.delta_v_loi_km_s > 0.6 && s.arrival.delta_v_loi_km_s < 1.2);
        // Азимут в секторе Восточного полигона, перигей орбиты перелёта — на опорной орбите.
        assert!(ksc.azimuth_allowed(s.launch.azimuth_deg), "az {}", s.launch.azimuth_deg);
        assert_abs_diff_eq!(s.departure.transfer_perigee_km, s.parking_orbit.radius_km, epsilon = 50.0);
        assert!(s.departure.c3_km2_s2 < 0.0, "transfer to the Moon is a bound ellipse");
        // Отчёт содержит перебор и счётчики.
        assert!(!fwd.coarse_scan.is_empty() && !fwd.fine_scan.is_empty());
        assert!(fwd.stats.lambert_solves > 0);

        // Обратная задача: по дате прилёта найденного варианта он же (или не хуже) должен находиться.
        let (y, m, d, ..) = s.arrival.epoch.to_gregorian_utc();
        let arr_day = Epoch::from_gregorian_utc_at_midnight(y, m, d);
        let inv = dt.solve(&almanac, ksc, DateConstraint::ArrivalDay(arr_day), None, None).unwrap();
        let si = inv.solution.expect("solution");
        assert!(si.arrival.epoch >= arr_day && si.arrival.epoch < arr_day + 1 * Unit::Day);
        assert!(si.delta_v_total_km_s <= s.delta_v_total_km_s + 0.005);
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_budget_controls_speed_and_feasibility() {
        let (almanac, catalog) = real();
        let ksc = catalog.get("ksc").unwrap();
        let dt = DirectTransfer::default();
        let day = Epoch::from_gregorian_utc_at_midnight(2026, 12, 1);
        let econ = dt.solve(&almanac, ksc, DateConstraint::LaunchDay(day), None, None).unwrap().solution.unwrap();

        // Большой запас Δv: самый быстрый перелёт, упирается в нижнюю границу времени.
        let big = DeltaVBudget::new("big", 5000.0, &MissionConfig::default()).unwrap();
        let fast = dt.solve(&almanac, ksc, DateConstraint::LaunchDay(day), Some(&big), None).unwrap().solution.unwrap();
        let b = fast.budget.clone().unwrap();
        assert!(fast.departure.time_of_flight_s < econ.departure.time_of_flight_s);
        assert!(b.time_of_flight_at_lower_limit);
        assert!(b.margin_m_s >= 0.0);
        assert!(fast.arrival.delta_v_loi_km_s > econ.arrival.delta_v_loi_km_s);
        assert!(fast.departure.delta_v_magnitude_km_s <= DirectTransferConfig::default().max_tli_delta_v_km_s);

        // Бюджет между экономичным и быстрым: перелёт быстрее экономичного, но в пределах запаса.
        let mid_need = (econ.arrival.delta_v_loi_km_s * 1000.0 + DEFAULT_DESCENT_DELTA_V_M_S) * 1.05 + 60.0;
        let mid = DeltaVBudget::new("mid", mid_need, &MissionConfig::default()).unwrap();
        let m = dt.solve(&almanac, ksc, DateConstraint::LaunchDay(day), Some(&mid), None).unwrap().solution.unwrap();
        assert!(m.budget.as_ref().unwrap().margin_m_s >= 0.0);
        assert!(m.departure.time_of_flight_s <= econ.departure.time_of_flight_s + 1 * Unit::Second);

        // Слишком маленький бюджет: перелёта нет, причина названа.
        let tiny = DeltaVBudget::new("tiny", 1500.0, &MissionConfig::default()).unwrap();
        let none = dt.solve(&almanac, ksc, DateConstraint::LaunchDay(day), Some(&tiny), None).unwrap();
        assert!(none.solution.is_none());
        assert_eq!(none.no_solution_reason, Some(REASON_BUDGET));
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_landing_site_drives_the_dates() {
        use crate::lunar_landing::{LandingTarget, LandingWindowConfig};
        use crate::visibility::{SiteVisibility, VisibilityConfig};
        let (almanac, catalog) = real();
        let ksc = catalog.get("ksc").unwrap();
        let dt = DirectTransfer::default();
        let target = |lat: f64, lon: f64| LandingTarget {
            site: SiteVisibility::new(None, lat, lon, 0.0, VisibilityConfig::default()).unwrap(),
            window: LandingWindowConfig::default(),
            descent: crate::lunar_landing::DescentModel::default(),
        };
        let day = Epoch::from_gregorian_utc_at_midnight(2027, 1, 1);

        // Запад видимой стороны у линии заката: посадка в заданные сутки, условия выполнены.
        let west = target(5.0, -10.0);
        let rep = dt.solve(&almanac, ksc, DateConstraint::ArrivalDay(day), None, Some(&west)).unwrap();
        let s = rep.solution.expect("landing on New Year's day in the west");
        let l = s.landing.expect("landing info");
        assert!(l.epoch >= day && l.epoch < day + 1 * Unit::Day, "landing {}", l.epoch);
        assert!(l.epoch >= s.arrival.epoch);
        assert!(l.sun.visible && l.earth.visible && l.sun.elevation_deg <= 20.0);

        // Море Спокойствия в эти сутки за линией заката: посадки нет, и причина — Солнце.
        let east = target(0.674, 23.473);
        let rep = dt.solve(&almanac, ksc, DateConstraint::ArrivalDay(day), None, Some(&east)).unwrap();
        assert!(rep.solution.is_none());
        assert!(rep.no_solution_reason.is_some_and(|r| r.contains("Sun")), "{:?}", rep.no_solution_reason);

        // Разные точки — разные даты старта при одной дате посадки (координаты влияют на даты).
        let west2 = target(0.0, -30.0);
        let s2 = dt.solve(&almanac, ksc, DateConstraint::ArrivalDay(day), None, Some(&west2)).unwrap().solution.unwrap();
        assert_ne!(s.launch.epoch, s2.launch.epoch);
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_spaceport_limits_are_respected() {
        let (almanac, catalog) = real();
        let dt = DirectTransfer::default();
        let day = Epoch::from_gregorian_utc_at_midnight(2026, 12, 1);
        for port in catalog.iter() {
            let rep = dt.solve(&almanac, port, DateConstraint::LaunchDay(day), None, None).unwrap();
            let Some(s) = rep.solution else { continue };
            assert!(port.azimuth_allowed(s.launch.azimuth_deg), "{}: az {}", port.id, s.launch.azimuth_deg);
            let i = s.parking_orbit.inclination_deg;
            if let (Some(lo), Some(hi)) = (port.inclination_min_deg, port.inclination_max_deg) {
                assert!(i >= lo - 1e-9 && i <= hi + 1e-9, "{}: incl {i}", port.id);
            }
            if let Some(list) = &port.allowed_inclinations_deg {
                assert!(list.iter().any(|c| (i - c).abs() <= 1.0 + 1e-9), "{}: incl {i}", port.id);
            }
            assert!(i + 1e-9 >= port.lat_deg.abs(), "{}: incl {i} below latitude", port.id);
        }
    }
}
