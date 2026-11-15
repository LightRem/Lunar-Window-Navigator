//! Перелёт Земля -> Луна: связь даты старта и даты прилёта.
//!
//! «Обратное» направление означает обратное отображение дат (дата прилёта -> дата старта),
//! а не полёт назад к Земле.
//!
//! # Типы траекторий
//!
//! Тип выбирается при запуске (`--trajectory`, поле `trajectory` в API), см. [`TrajectoryKind`]:
//!
//! - [`DirectTransfer`] (`direct`) — прямой перелёт: задача Ламберта при любом угле перелёта и
//!   времени полёта в диапазоне `tof_min_days..tof_max_days` ([`direct`]);
//! - [`HohmannTransfer`] (`hohmann`) — перелёт Гомана: разгон по касательной в перигее,
//!   апогей на расстоянии Луны, угол перелёта 180°, время полёта — половина периода ([`hohmann`]).
//!
//! Общая часть поиска — в [`common`]:
//!
//! 1. Старт с космодрома на круговую опорную орбиту, плоскость которой проходит через
//!    положение Луны в момент прилёта. Наклонение — наименьшее, которое допускают широта
//!    космодрома, склонение Луны и ограничения космодрома (азимуты, коридоры наклонений).
//! 2. Старт возможен, когда космодром из-за вращения Земли проходит через эту плоскость
//!    (дважды в сутки), а азимут пуска в этот момент входит в разрешённый сектор.
//!    Вращение Земли берётся из ядер ANISE (ITRF93).
//! 3. Подъём на опорную орбиту, полёт по ней от `coast_min` до `coast_max`, затем разгон к Луне (TLI).
//! 4. Перелёт до Луны — своя модель у каждого типа траектории ([`common::LegModel`]).
//! 5. У Луны — торможение на круговую окололунную орбиту, посадка в точку (если задана).
//! 6. Из всех вариантов выбирается лучший по критерию (минимум Δv или самый быстрый).
//!
//! Новый тип траектории — это новая реализация [`common::LegModel`] (расчёт разгона для прохода
//! космодрома через плоскость), вариант в [`TrajectoryKind`] и секция в `transfer.toml`.
//!
//! Упрощения (первое приближение, «patched conic»): притяжение Луны на перелёте и Солнца
//! не учитывается, момент прилёта — это момент встречи с центром Луны в этой модели.
//! Лучший вариант затем уточняется численным интегрированием ([`crate::refine`]).

pub mod common;
pub mod direct;
pub mod hohmann;

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anise::math::Vector3;
use anise::prelude::Almanac;
use anise::time::{Duration, Epoch};
use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};
use crate::lunar_landing::{DescentConfig, LandingInfo, LandingTarget, LandingWindowConfig};
use crate::refine::{NumericalConfig, Refined};
use crate::runlog::{ser_duration_s, ser_epoch, ser_vec3};
use crate::sky::LunarSite;
use crate::spaceport::Spaceport;

pub use common::{AscentSetup, LegModel, Transfer, TransferCommonConfig, geodetic_to_ecef};
pub use direct::{DirectTransfer, DirectTransferConfig};
pub use hohmann::{HohmannConfig, HohmannInfo, HohmannTransfer};

/// Δv спуска с низкой окололунной орбиты на поверхность по умолчанию, м/с
/// (номинальная проектная величина, AIAA-2008-6939).
pub const DEFAULT_DESCENT_DELTA_V_M_S: f64 = 2053.0;
/// Причина «не хватает Δv аппарата» (самая содержательная из причин перебора).
pub const REASON_BUDGET: &str = "the spacecraft's delta-v budget is not enough for any transfer arriving at this time";

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
    /// Тип траектории: "direct_transfer" или "hohmann".
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
    /// Данные перелёта Гомана (только для траектории "hohmann").
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hohmann: Option<HohmannInfo>,
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
    /// Варианты, отброшенные моделью траектории (у Гомана — разгон не в тот момент:
    /// время полёта по эллипсу Гомана не совпадает с моментом прилёта в пределах допуска).
    pub candidates_rejected_timing: u64,
    pub arrivals_rejected_landing: u64,
    pub lambert_solves: u64,
    pub lambert_failures: u64,
    pub candidates_evaluated: u64,
}

/// Полный отчёт о поиске перелёта.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TransferReport {
    /// Тип траектории: "direct_transfer" или "hohmann".
    pub trajectory: &'static str,
    pub spaceport: Spaceport,
    pub constraint: DateConstraint,
    /// Параметры поиска: `{ "transfer": общие, "<тип траектории>": свои }`.
    pub config: serde_json::Value,
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
    /// Общие параметры всех траекторий (секция `[transfer]`).
    pub transfer: TransferCommonConfig,
    /// Свои параметры прямого перелёта (секция `[direct_transfer]`).
    pub direct_transfer: DirectTransferConfig,
    /// Свои параметры перелёта Гомана (секция `[hohmann]`).
    pub hohmann: HohmannConfig,
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
    pub(crate) fn required_m_s(&self, dv_tli_km_s: f64, dv_loi_km_s: f64) -> f64 {
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

    pub(crate) fn fits(&self, dv_tli_km_s: f64, dv_loi_km_s: f64) -> bool {
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
        check_moved_parameters(text, origin)?;
        let cfg: Self = toml::from_str(text)
            .map_err(|source| Error::Toml { path: origin.to_path_buf(), source })?;
        cfg.transfer.validate()?;
        cfg.direct_transfer.validate(&cfg.transfer)?;
        cfg.hohmann.validate(&cfg.transfer)?;
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

/// Параметры, которые раньше были в `[direct_transfer]`, а теперь общие (`[transfer]`):
/// понятная ошибка вместо «unknown field» для старого файла.
fn check_moved_parameters(text: &str, origin: &Path) -> Result<()> {
    let Ok(doc) = text.parse::<toml::Table>() else { return Ok(()) };
    let Some(old) = doc.get("direct_transfer").and_then(|v| v.as_table()) else { return Ok(()) };
    let moved: Vec<&str> = old.keys().map(String::as_str).filter(|k| TransferCommonConfig::FIELDS.contains(k)).collect();
    if moved.is_empty() {
        return Ok(());
    }
    Err(Error::Config(format!(
        "{}: parameters {} moved from [direct_transfer] to the common section [transfer] (shared by all trajectories); \
         [direct_transfer] keeps only tof_min_days and tof_max_days",
        origin.display(),
        moved.join(", ")
    )))
}

// ---------------------------------------------------------------------------
// Выбор траектории
// ---------------------------------------------------------------------------

/// Тип траектории перелёта.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TrajectoryKind {
    /// Прямой перелёт (задача Ламберта), см. [`DirectTransfer`].
    #[default]
    #[serde(alias = "direct_transfer")]
    Direct,
    /// Перелёт Гомана, см. [`HohmannTransfer`].
    Hohmann,
}

impl TrajectoryKind {
    pub const ALL: [TrajectoryKind; 2] = [TrajectoryKind::Direct, TrajectoryKind::Hohmann];

    /// Имя для командной строки и API.
    pub fn as_str(self) -> &'static str {
        match self {
            TrajectoryKind::Direct => "direct",
            TrajectoryKind::Hohmann => "hohmann",
        }
    }

    /// Имя траектории в отчётах и журналах (`TransferReport::trajectory`).
    pub fn report_name(self) -> &'static str {
        match self {
            TrajectoryKind::Direct => <DirectTransferConfig as LegModel>::NAME,
            TrajectoryKind::Hohmann => <HohmannConfig as LegModel>::NAME,
        }
    }

    /// Человекочитаемое описание модели.
    pub fn description(self) -> &'static str {
        match self {
            TrajectoryKind::Direct => <DirectTransferConfig as LegModel>::DESCRIPTION,
            TrajectoryKind::Hohmann => <HohmannConfig as LegModel>::DESCRIPTION,
        }
    }
}

impl std::fmt::Display for TrajectoryKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl std::str::FromStr for TrajectoryKind {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, String> {
        match s.trim().to_ascii_lowercase().replace('-', "_").as_str() {
            "direct" | "direct_transfer" => Ok(TrajectoryKind::Direct),
            "hohmann" | "hohmann_transfer" | "homan" | "goman" => Ok(TrajectoryKind::Hohmann),
            other => Err(format!(
                "unknown trajectory '{other}' (known: {})",
                TrajectoryKind::ALL.map(|k| k.as_str()).join(", ")
            )),
        }
    }
}

/// Создаёт траекторию выбранного типа по параметрам из `transfer.toml`.
///
/// `ascent` — ракеты для моделирования подъёма по космодромам, `ascent_cache_dir` — папка кэша
/// подъёма (`None` — без кэша).
pub fn build(
    kind: TrajectoryKind,
    cfg: &TransferConfig,
    ascent: HashMap<String, AscentSetup>,
    ascent_cache_dir: Option<PathBuf>,
) -> Result<Box<dyn Trajectory>> {
    Ok(match kind {
        TrajectoryKind::Direct => Box::new(
            DirectTransfer::new(cfg.transfer.clone(), cfg.direct_transfer.clone())?
                .with_ascent(ascent)
                .with_ascent_cache_dir(ascent_cache_dir),
        ),
        TrajectoryKind::Hohmann => Box::new(
            HohmannTransfer::new(cfg.transfer.clone(), cfg.hohmann.clone())?
                .with_ascent(ascent)
                .with_ascent_cache_dir(ascent_cache_dir),
        ),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn config_defaults_and_validation() {
        let c = TransferConfig::parse("", Path::new("t.toml")).unwrap();
        assert_eq!(c, TransferConfig::default());
        assert!(TransferConfig::parse("[direct_transfer]\ntof_min_days = 6.0\n", Path::new("t.toml")).is_err());
        assert!(TransferConfig::parse("[transfer]\nparking_alt = 200\n", Path::new("t.toml")).is_err());
        assert!(TransferConfig::parse("[hohmann]\ntiming_tolerance_min = 0.0\n", Path::new("t.toml")).is_err());
        let ok = TransferConfig::parse("[transfer]\nparking_altitude_km = 200\n[hohmann]\ntiming_tolerance_min = 5\n", Path::new("t.toml")).unwrap();
        assert_eq!(ok.transfer.parking_altitude_km, 200.0);
        assert_eq!(ok.hohmann.timing_tolerance_min, 5.0);
    }

    #[test]
    fn old_layout_gives_a_clear_error() {
        let e = TransferConfig::parse("[direct_transfer]\nparking_altitude_km = 185\ntof_min_days = 2.5\n", Path::new("t.toml"))
            .unwrap_err()
            .to_string();
        assert!(e.contains("parking_altitude_km") && e.contains("[transfer]"), "{e}");
    }

    #[test]
    fn trajectory_kind_names() {
        for k in TrajectoryKind::ALL {
            assert_eq!(k.as_str().parse::<TrajectoryKind>().unwrap(), k);
            assert_eq!(k.report_name().parse::<TrajectoryKind>().unwrap(), k);
            let json = serde_json::to_string(&k).unwrap();
            assert_eq!(serde_json::from_str::<TrajectoryKind>(&json).unwrap(), k);
        }
        assert_eq!("Direct-Transfer".parse::<TrajectoryKind>().unwrap(), TrajectoryKind::Direct);
        assert!("lunar_flyby".parse::<TrajectoryKind>().is_err());
        assert_eq!(TrajectoryKind::default(), TrajectoryKind::Direct);
        let cfg = TransferConfig::default();
        for k in TrajectoryKind::ALL {
            assert_eq!(build(k, &cfg, HashMap::new(), None).unwrap().name(), k.report_name());
        }
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
        let text = include_str!("../../data/transfer.toml");
        TransferConfig::parse(text, Path::new("data/transfer.toml")).unwrap();
    }
}
