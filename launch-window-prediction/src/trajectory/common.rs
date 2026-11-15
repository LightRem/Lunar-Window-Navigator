//! Общая часть поиска перелёта для всех типов траекторий.
//!
//! Здесь всё, что не зависит от формы перелёта: окна поиска по дате, трек космодрома
//! (вращение Земли), плоскости опорной орбиты через Луну, проходы космодрома через плоскость,
//! азимуты, подъём ракеты (и его кэш), перебор моментов прилёта, торможение у Луны, посадка
//! в точку, бюджет Δv аппарата, причины отказа и сборка [`TransferSolution`].
//!
//! Форма перелёта задаётся реализацией [`LegModel`]: по проходу космодрома через плоскость
//! она даёт варианты разгона к Луне (момент, положение и скорости до и после разгона, скорость
//! у Луны). [`Transfer<M>`] реализует [`Trajectory`] для любой такой модели.

use std::sync::atomic::{AtomicU64, Ordering};

use anise::constants::frames::{EARTH_ITRF93, EARTH_J2000, IAU_EARTH_FRAME, MOON_J2000};
use anise::math::Vector3;
use anise::prelude::{Almanac, Frame};
use anise::time::{Duration, Epoch, Unit};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use super::{
    AscentInfo, ArrivalInfo, DateConstraint, DeltaVBudget, DepartureInfo, HohmannInfo, InclinationChoice,
    LaunchInfo, Maneuver, Objective, ParkingOrbit, REASON_BUDGET, ScanPoint, SearchStats, Trajectory, TransferReport,
    TransferSolution,
};
use crate::error::{Error, Result};
use crate::lunar_landing::{LandingInfo, LandingTarget, find_landing};
use crate::spaceport::Spaceport;

/// Гравитационные параметры по умолчанию, км³/с² (если в ядрах их нет).
pub(crate) const MU_EARTH_DEFAULT: f64 = 398_600.435_436;
pub(crate) const MU_MOON_DEFAULT: f64 = 4_902.800_066;
/// Экваториальный радиус Земли и сжатие WGS84.
pub(crate) const EARTH_RADIUS_KM: f64 = 6_378.137;
const EARTH_FLATTENING: f64 = 1.0 / 298.257_223_563;
/// Радиус Луны, км (опорная сфера LOLA).
const MOON_RADIUS_KM: f64 = 1_737.4;
/// Скорость вращения Земли на экваторе, км/с.
const EARTH_EQUATOR_SPEED_KM_S: f64 = 0.465_1;
/// Наибольшее наклонение опорной орбиты (только прямые орбиты), градусы.
const MAX_INCLINATION_DEG: f64 = 89.0;
/// Шаг перебора наклонений, если космодром задаёт диапазон или ничего не задаёт, градусы.
const INCLINATION_STEP_DEG: f64 = 1.0;
const REASON_WINDOW: &str = "launch window does not overlap this arrival time";
/// Шаг уточняющего перебора момента прилёта, с.
const FINE_ARRIVAL_STEP_S: f64 = 60.0;

// ---------------------------------------------------------------------------
// Общие параметры (секция [transfer])
// ---------------------------------------------------------------------------

/// Общие параметры перелёта для всех траекторий (секция `[transfer]` в `transfer.toml`).
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct TransferCommonConfig {
    /// Высота круговой опорной орбиты над экватором Земли, км.
    pub parking_altitude_km: f64,
    /// Высота круговой окололунной орбиты после торможения, км.
    pub lunar_orbit_altitude_km: f64,
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
    /// Часть суток старта (UTC), в которой ищется старт при `--launch-date`: от и до, часы 0..24.
    /// По умолчанию вся сутки. Не влияет на поиск по дате прилёта.
    pub launch_from_utc_h: f64,
    pub launch_to_utc_h: f64,
}

impl Default for TransferCommonConfig {
    fn default() -> Self {
        Self {
            parking_altitude_km: 185.0,
            lunar_orbit_altitude_km: 100.0,
            coast_min_min: 20.0,
            coast_max_min: 180.0,
            inclination_margin_deg: 0.5,
            search_step_min: 30.0,
            max_tli_delta_v_km_s: 3.3,
            // Типичный подъём к Луне: Saturn V — 11 мин 39 с до выхода на орбиту,
            // SLS — около 8 мин до отделения центральной ступени. Дальность ~2000 км.
            ascent_duration_min: 10.0,
            ascent_downrange_km: 2000.0,
            launch_from_utc_h: 0.0,
            launch_to_utc_h: 24.0,
        }
    }
}

impl TransferCommonConfig {
    /// Имена параметров секции (для подсказки при старом формате файла).
    pub const FIELDS: [&'static str; 11] = [
        "parking_altitude_km",
        "lunar_orbit_altitude_km",
        "coast_min_min",
        "coast_max_min",
        "inclination_margin_deg",
        "search_step_min",
        "max_tli_delta_v_km_s",
        "ascent_duration_min",
        "ascent_downrange_km",
        "launch_from_utc_h",
        "launch_to_utc_h",
    ];

    pub fn validate(&self) -> Result<()> {
        let bad = |m: &str| Err(Error::Config(format!("transfer: {m}")));
        if !(self.parking_altitude_km >= 120.0 && self.parking_altitude_km <= 2000.0) {
            return bad("parking_altitude_km must be in [120, 2000]");
        }
        if !(self.lunar_orbit_altitude_km >= 10.0 && self.lunar_orbit_altitude_km <= 10_000.0) {
            return bad("lunar_orbit_altitude_km must be in [10, 10000]");
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
        if !(self.launch_from_utc_h >= 0.0 && self.launch_to_utc_h > self.launch_from_utc_h && self.launch_to_utc_h <= 24.0) {
            return bad("need 0 <= launch_from_utc_h < launch_to_utc_h <= 24");
        }
        if self.coast_min_min < self.ascent_duration_min {
            return bad("coast_min_min (time from launch to TLI) must be at least ascent_duration_min");
        }
        Ok(())
    }

    /// Радиус опорной орбиты, км.
    pub fn parking_radius_km(&self) -> f64 {
        EARTH_RADIUS_KM + self.parking_altitude_km
    }
}

// ---------------------------------------------------------------------------
// Модель перелёта (что отличает типы траекторий)
// ---------------------------------------------------------------------------

/// Проход космодрома через плоскость опорной орбиты для заданного момента прилёта:
/// всё, что нужно модели перелёта, чтобы предложить варианты разгона к Луне.
///
/// Векторы — в J2000 (относительно центра Земли), км и км/с; время — в секундах.
#[derive(Debug, Clone, Copy)]
pub struct CrossingGeometry {
    /// Момент прилёта к Луне.
    pub arrival: Epoch,
    /// Момент старта (проход космодрома через плоскость).
    pub launch: Epoch,
    /// Положение и скорость Луны в момент прилёта.
    pub moon_position_km: Vector3,
    pub moon_velocity_km_s: Vector3,
    /// Единичная нормаль плоскости опорной орбиты (движение против часовой стрелки с её конца).
    pub plane_normal: Vector3,
    /// Направление на космодром в момент старта, спроецированное на плоскость (единичное).
    pub l_hat: Vector3,
    /// Направление движения в плоскости: `plane_normal × l_hat`.
    pub along: Vector3,
    /// Положение на опорной орбите через `t` с после старта — угол `ascent_phase + n_park·t`
    /// от `l_hat` к `along` (после выхода на орбиту).
    pub ascent_phase_rad: f64,
    /// Длительность подъёма до выхода на опорную орбиту, с.
    pub ascent_s: f64,
    /// Среднее движение, радиус и круговая скорость опорной орбиты.
    pub n_park_rad_s: f64,
    pub parking_radius_km: f64,
    pub circular_speed_km_s: f64,
    pub mu_earth_km3_s2: f64,
    /// Допустимое время от старта до разгона, с.
    pub coast_min_s: f64,
    pub coast_max_s: f64,
    /// Шаг перебора времени на опорной орбите, с (для моделей, которые его перебирают).
    pub coast_step_s: f64,
    /// Шаг перебора момента прилёта, с: момент прилёта известен с этой точностью
    /// (модели с жёсткой связью времён используют его для допуска).
    pub arrival_step_s: f64,
}

impl CrossingGeometry {
    /// От старта до прилёта, с.
    pub fn launch_to_arrival_s(&self) -> f64 {
        (self.arrival - self.launch).to_seconds()
    }

    /// Единичный вектор положения на опорной орбите при угле `ang` от `l_hat`.
    pub fn parking_direction(&self, ang: f64) -> Vector3 {
        self.l_hat * ang.cos() + self.along * ang.sin()
    }
}

/// Вариант разгона к Луне, предложенный моделью перелёта.
#[derive(Debug, Clone, Copy)]
pub struct Leg {
    /// Время от старта до разгона, с.
    pub coast_s: f64,
    /// Угол на опорной орбите от `l_hat` в момент разгона, рад.
    pub coast_angle_rad: f64,
    /// Положение в момент разгона и скорость на опорной орбите до него.
    pub r1: Vector3,
    pub v_park: Vector3,
    /// Скорость сразу после разгона.
    pub v1: Vector3,
    /// Скорость в момент прилёта (относительно Земли).
    pub v2: Vector3,
    /// Угол перелёта, рад.
    pub transfer_angle_rad: f64,
    /// Данные перелёта Гомана (если это он).
    pub hohmann: Option<HohmannInfo>,
}

/// Счётчики модели перелёта для одного прохода.
#[derive(Debug, Clone, Copy, Default)]
pub struct LegCounters {
    /// Решения задачи Ламберта (или другого расчёта перелёта) и неудачи.
    pub solves: u64,
    pub failures: u64,
    /// Варианты, отброшенные моделью по времени (у Гомана — разгон не в тот момент).
    pub rejected_timing: u64,
}

/// Форма перелёта от опорной орбиты до Луны: то, чем отличаются типы траекторий.
///
/// Реализация — это параметры траектории (секция в `transfer.toml`), например
/// [`super::DirectTransferConfig`] или [`super::HohmannConfig`].
pub trait LegModel: std::fmt::Debug + Clone + Default + Serialize + Send + Sync + 'static {
    /// Имя траектории в отчётах (`TransferReport::trajectory`).
    const NAME: &'static str;
    /// Имя секции параметров в `transfer.toml`.
    const SECTION: &'static str;
    /// Описание модели для вывода.
    const DESCRIPTION: &'static str;

    fn validate(&self, common: &TransferCommonConfig) -> Result<()>;

    /// Допустимое время перелёта от разгона до прилёта, с: от него зависят окна поиска.
    fn time_of_flight_bounds_s(&self, common: &TransferCommonConfig) -> (f64, f64);

    /// Варианты разгона для одного прохода космодрома через плоскость: каждый передаётся в `out`.
    fn legs(&self, g: &CrossingGeometry, counters: &mut LegCounters, out: &mut dyn FnMut(Leg));

    /// Причина, если все варианты для момента прилёта отброшены моделью по времени.
    fn timing_reason(&self) -> &'static str {
        "no launch opportunity within the coast and flight-time limits"
    }

    /// Насколько сдвинуть момент прилёта, чтобы вариант точнее подходил модели, с
    /// (`None` — сдвигать не нужно). Используется для доводки лучшего варианта.
    fn arrival_correction_s(&self, _leg: &Leg) -> Option<f64> {
        None
    }
}

// ---------------------------------------------------------------------------
// Поиск перелёта для любой модели
// ---------------------------------------------------------------------------

/// Перелёт Земля -> Луна с формой перелёта `M` (см. [`super::DirectTransfer`],
/// [`super::HohmannTransfer`]).
#[derive(Debug, Clone, Default)]
pub struct Transfer<M: LegModel> {
    /// Общие параметры (секция `[transfer]`).
    pub config: TransferCommonConfig,
    /// Параметры формы перелёта (своя секция).
    pub model: M,
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

impl<M: LegModel> Transfer<M> {
    pub fn new(config: TransferCommonConfig, model: M) -> Result<Self> {
        config.validate()?;
        model.validate(&config)?;
        Ok(Self { config, model, ascent: Default::default(), ascent_cache: Default::default(), ascent_cache_dir: None })
    }

    /// Параметры поиска для отчёта: общие и свои.
    fn config_json(&self) -> serde_json::Value {
        let mut m = serde_json::Map::new();
        m.insert("transfer".into(), serde_json::to_value(&self.config).unwrap_or_default());
        m.insert(M::SECTION.into(), serde_json::to_value(&self.model).unwrap_or_default());
        serde_json::Value::Object(m)
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
    hohmann: Option<HohmannInfo>,
}

impl Candidate {
    fn total(&self) -> f64 {
        self.dv_tli + self.dv_loi
    }

    fn leg(&self) -> Leg {
        Leg {
            coast_s: (self.tli - self.launch).to_seconds(),
            coast_angle_rad: self.coast_angle,
            r1: self.r1,
            v_park: self.v_park,
            v1: self.v1,
            v2: self.v2,
            transfer_angle_rad: self.angle,
            hohmann: self.hohmann,
        }
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
    candidates_rejected_timing: AtomicU64,
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
            candidates_rejected_timing: g(&self.candidates_rejected_timing),
            arrivals_rejected_landing: g(&self.arrivals_rejected_landing),
            lambert_solves: g(&self.lambert_solves),
            lambert_failures: g(&self.lambert_failures),
            candidates_evaluated: g(&self.candidates_evaluated),
        }
    }
}

struct Search<'a, M: LegModel> {
    almanac: &'a Almanac,
    port: &'a Spaceport,
    cfg: &'a TransferCommonConfig,
    transfer: &'a Transfer<M>,
    site: &'a SiteTrack,
    mu_e: f64,
    mu_m: f64,
    counters: Counters,
    budget: Option<&'a DeltaVBudget>,
    landing: Option<&'a LandingTarget>,
    /// Допустимый интервал момента посадки (сутки посадки для ArrivalDay).
    landing_day: Option<(Epoch, Epoch)>,
}

impl<M: LegModel> Search<'_, M> {
    fn objective(&self) -> Objective {
        self.budget.map_or(Objective::MinDeltaV, |b| b.objective)
    }

    /// Лучший вариант для прилёта в момент `t_a`, либо причина, почему вариантов нет.
    ///
    /// `arrival_step` — с какой точностью известен момент прилёта (шаг перебора); модели с жёсткой
    /// связью времён (Гоман) используют его для допуска.
    fn at_arrival(
        &self,
        t_a: Epoch,
        coast_step: Duration,
        arrival_step: Duration,
        (launch_lo, launch_hi): (Epoch, Epoch),
    ) -> Result<std::result::Result<Candidate, &'static str>> {
        let cfg = self.cfg;
        let model = &self.transfer.model;
        Counters::add(&self.counters.arrival_epochs, 1);
        let moon = self
            .almanac
            .translate(MOON_J2000, EARTH_J2000, t_a, None)
            .map_err(|e| Error::Almanac(format!("Moon at {t_a}: {e}")))?;
        let (r2, v_moon) = (moon.radius_km, moon.velocity_km_s);
        let moon_dec = (r2.z / r2.norm()).asin();

        let (tof_min_s, tof_max_s) = model.time_of_flight_bounds_s(cfg);
        let tof_min = tof_min_s * Unit::Second;
        let tof_max = tof_max_s * Unit::Second;
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

        let r_p = cfg.parking_radius_km();
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
        let mut rejected_by_timing = false;
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
                    let geometry = CrossingGeometry {
                        arrival: t_a,
                        launch: t_l,
                        moon_position_km: r2,
                        moon_velocity_km_s: v_moon,
                        plane_normal: h,
                        l_hat,
                        along,
                        ascent_phase_rad: asc_angle - n_park * asc_t,
                        ascent_s: asc_t,
                        n_park_rad_s: n_park,
                        parking_radius_km: r_p,
                        circular_speed_km_s: v_circ,
                        mu_earth_km3_s2: self.mu_e,
                        coast_min_s: coast_min.to_seconds(),
                        coast_max_s: coast_max.to_seconds(),
                        coast_step_s: coast_step.to_seconds(),
                        arrival_step_s: arrival_step.to_seconds(),
                    };
                    // Внутренний цикл — с локальными счётчиками: общие атомарные счётчики здесь
                    // были дороже самой задачи Ламберта.
                    let mut legc = LegCounters::default();
                    let (mut n_tli, mut n_budget) = (0u64, 0u64);
                    model.legs(&geometry, &mut legc, &mut |leg: Leg| {
                        let dv_tli = (leg.v1 - leg.v_park).norm();
                        if dv_tli > cfg.max_tli_delta_v_km_s {
                            n_tli += 1;
                            return;
                        }
                        let v_inf = (leg.v2 - v_moon).norm();
                        let dv_loi = (v_inf * v_inf + 2.0 * self.mu_m / r_llo).sqrt() - v_llo;
                        if let Some(b) = self.budget
                            && !b.fits(dv_tli, dv_loi)
                        {
                            n_budget += 1;
                            return;
                        }
                        let cand = Candidate {
                            launch: t_l,
                            tli: t_l + leg.coast_s * Unit::Second,
                            arrival: t_a,
                            incl,
                            h,
                            site_dir: s,
                            azimuth_deg: az,
                            ascending: along.z > 0.0,
                            coast_angle: leg.coast_angle_rad,
                            r1: leg.r1,
                            v_park: leg.v_park,
                            v1: leg.v1,
                            v2: leg.v2,
                            r2,
                            v_moon,
                            angle: leg.transfer_angle_rad,
                            dv_tli,
                            v_inf,
                            dv_loi,
                            moon_dec,
                            required_min_deg,
                            landing: None,
                            hohmann: leg.hohmann,
                        };
                        if best.is_none_or(|b| cand.better_than(&b, objective)) {
                            best = Some(cand);
                        }
                    });
                    Counters::add(&self.counters.lambert_solves, legc.solves);
                    Counters::add(&self.counters.lambert_failures, legc.failures);
                    Counters::add(&self.counters.candidates_evaluated, legc.solves - legc.failures);
                    Counters::add(&self.counters.candidates_rejected_tli_limit, n_tli);
                    Counters::add(&self.counters.candidates_rejected_budget, n_budget);
                    Counters::add(&self.counters.candidates_rejected_timing, legc.rejected_timing);
                    rejected_by_tli |= n_tli > 0;
                    rejected_by_budget |= n_budget > 0;
                    rejected_by_timing |= legc.rejected_timing > 0;
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
        } else if rejected_by_timing {
            model.timing_reason()
        } else {
            "no launch opportunity within the coast and flight-time limits"
        }))
    }

    /// Доводка варианта для моделей с жёсткой связью времён (Гоман): момент прилёта
    /// сдвигается, пока вариант не подойдёт модели точно (несколько шагов, каждый —
    /// пересчёт для нового момента прилёта). Возвращает лучший из найденных.
    fn polish(&self, mut best: Candidate, coast_step: Duration, launch: (Epoch, Epoch), (ta_lo, ta_hi): (Epoch, Epoch)) -> Result<Candidate> {
        let model = &self.transfer.model;
        for _ in 0..6 {
            let Some(shift_s) = model.arrival_correction_s(&best.leg()) else { break };
            let t = best.arrival + shift_s * Unit::Second;
            if t < ta_lo || t > ta_hi {
                break;
            }
            let step = FINE_ARRIVAL_STEP_S * Unit::Second;
            match self.at_arrival(t, coast_step, step, launch)? {
                Ok(c) => {
                    let err = |c: &Candidate| model.arrival_correction_s(&c.leg()).map_or(0.0, f64::abs);
                    if err(&c) < err(&best) {
                        best = c;
                    } else {
                        break;
                    }
                }
                Err(_) => break,
            }
        }
        Ok(best)
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
                Ok((t_a, self.at_arrival(t_a, coast_step, ta_step, launch)?))
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

    fn solution(&self, c: &Candidate) -> Result<TransferSolution> {
        let name = M::NAME;
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
        let r_p = cfg.parking_radius_km();
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
            hohmann: c.hohmann,
            mu_earth_km3_s2: self.mu_e,
            mu_moon_km3_s2: self.mu_m,
            landing: c.landing,
            refined: None,
            refine_failure: None,
            refine_skipped: None,
            budget: self.budget.map(|b| {
                let (tof_min_s, _) = self.transfer.model.time_of_flight_bounds_s(cfg);
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

impl<M: LegModel> Transfer<M> {
    fn search<'a>(
        &'a self,
        almanac: &'a Almanac,
        spaceport: &'a Spaceport,
        site: &'a SiteTrack,
        budget: Option<&'a DeltaVBudget>,
        landing: Option<&'a LandingTarget>,
        landing_day: Option<(Epoch, Epoch)>,
    ) -> Search<'a, M> {
        Search {
            almanac,
            port: spaceport,
            cfg: &self.config,
            transfer: self,
            site,
            mu_e: gm(almanac, EARTH_J2000, MU_EARTH_DEFAULT),
            mu_m: gm(almanac, MOON_J2000, MU_MOON_DEFAULT),
            counters: Counters::default(),
            budget,
            landing,
            landing_day,
        }
    }
}

impl<M: LegModel> Trajectory for Transfer<M> {
    fn name(&self) -> &'static str {
        M::NAME
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
        self.model.validate(cfg)?;
        let w = self.windows(constraint, landing);
        let (launch_lo, launch_hi, ta_lo, ta_hi) = (w.launch_lo, w.launch_hi, w.ta_lo, w.ta_hi);
        let site = SiteTrack::build(almanac, spaceport, launch_lo, launch_hi)?;
        let search = self.search(almanac, spaceport, &site, budget, landing, w.landing_day);
        let ascent_loaded = self.load_ascent_cache(spaceport, search.mu_e);

        // Грубый перебор момента прилёта, затем уточнение вокруг лучшего.
        let coarse_step = cfg.search_step_min * Unit::Minute;
        let (coarse, coarse_scan) = search.scan(ta_lo, ta_hi, coarse_step, 2 * Unit::Minute, (launch_lo, launch_hi))?;
        let (best, fine_scan) = match coarse {
            Some(c) => {
                let fine_lo = (c.arrival - coarse_step).max(ta_lo);
                let fine_hi = (c.arrival + coarse_step).min(ta_hi);
                // Окно пуска не сужается: поле Δv пологое, и близкие по Δv варианты бывают в
                // других окнах пуска (сужение меняло ответ при разнице ~1 м/с).
                let (fine, scan) = search.scan(
                    fine_lo,
                    fine_hi,
                    FINE_ARRIVAL_STEP_S * Unit::Second,
                    30 * Unit::Second,
                    (launch_lo, launch_hi),
                )?;
                let best = fine.unwrap_or(c);
                (Some(search.polish(best, 30 * Unit::Second, (launch_lo, launch_hi), (ta_lo, ta_hi))?), scan)
            }
            None => (None, Vec::new()),
        };
        self.save_ascent_cache(spaceport, search.mu_e, ascent_loaded);

        Ok(TransferReport {
            trajectory: self.name(),
            spaceport: spaceport.clone(),
            constraint,
            config: self.config_json(),
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
            solution: best.map(|c| search.solution(&c)).transpose()?,
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
        self.model.validate(cfg)?;
        let w = self.windows(constraint, landing);
        let site = SiteTrack::build(almanac, spaceport, w.launch_lo, w.launch_hi)?;
        let search = self.search(almanac, spaceport, &site, budget, landing, w.landing_day);
        let launch = (w.launch_lo, w.launch_hi);
        // Тот же шаг по времени на опорной орбите, что и в уточняющем переборе. Момент прилёта
        // обычно взят из грубого перебора — с его точностью.
        let arrival_step = cfg.search_step_min * Unit::Minute;
        match search.at_arrival(arrival, 30 * Unit::Second, arrival_step, launch)? {
            Ok(c) => {
                let c = search.polish(c, 30 * Unit::Second, launch, (w.ta_lo, w.ta_hi))?;
                Ok(Some(search.solution(&c)?))
            }
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

impl<M: LegModel> Transfer<M> {
    fn windows(&self, constraint: DateConstraint, landing: Option<&LandingTarget>) -> Windows {
        let cfg = &self.config;
        let day = 1 * Unit::Day;
        let (tof_min_s, tof_max_s) = self.model.time_of_flight_bounds_s(cfg);
        let tof_min = tof_min_s * Unit::Second;
        let tof_max = tof_max_s * Unit::Second;
        let coast_min = cfg.coast_min_min * Unit::Minute;
        let coast_max = cfg.coast_max_min * Unit::Minute;
        // С точкой посадки «сутки прилёта» — это сутки посадки: прилёт к Луне может быть
        // раньше на время ожидания на орбите и спуска (полвитка эллипса + торможение, < 2 ч).
        let wait = landing.map_or(Duration::ZERO, |l| l.window.max_wait_in_orbit_h * Unit::Hour + 2 * Unit::Hour);
        let (launch_lo, launch_hi, ta_lo, ta_hi) = match constraint {
            DateConstraint::LaunchDay(d) => {
                let lo = d + cfg.launch_from_utc_h * Unit::Hour;
                let hi = d + cfg.launch_to_utc_h * Unit::Hour;
                (lo, hi, lo + coast_min + tof_min, hi + coast_max + tof_max)
            }
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
    use std::path::Path;

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
    fn common_config_validation() {
        assert!(TransferCommonConfig::default().validate().is_ok());
        let bad = TransferCommonConfig { coast_min_min: 5.0, ascent_duration_min: 10.0, ..Default::default() };
        assert!(bad.validate().is_err());
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_earth_rotation_matches_sidereal_time() {
        // Прямое восхождение космодрома ~ GMST + долгота (с точностью до прецессии с 2000 г., < 0.5 градуса).
        let almanac = load_almanac_blocking(AlmanacSource::LocalDir("data/kernels".into())).unwrap();
        let text = std::fs::read_to_string("data/spaceports.toml").unwrap();
        let catalog = SpaceportCatalog::parse(&text, Path::new("data/spaceports.toml")).unwrap();
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
}
