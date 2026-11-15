//! Видимость Солнца и Земли из точки на Луне с учётом рельефа.
//!
//! Это вход для подсистемы ресурсного баланса: освещённость (`h_S(t)`, затенение)
//! и прямая видимость Земли (`LOS_Earth(t)`).
//!
//! Все настраиваемые параметры собраны в [`VisibilityConfig`] и читаются из TOML-файла
//! (по умолчанию `data/visibility.toml`), поэтому менять правила можно без правки кода.
//!
//! Как устроено:
//! 1. [`SiteVisibility::new`] берёт высоту точки из карт рельефа и строит маску горизонта
//!    (угол возвышения рельефа по азимутам). Без карт горизонт считается ровным (0°).
//! 2. [`crate::sky`] даёт азимут, угол возвышения и дальность Солнца/Земли на каждый момент.
//! 3. [`SiteVisibility::evaluate`] сравнивает их с горизонтом по правилу [`BodyRule`].
//! 4. [`VisibilitySummary::from_samples`] считает доли времени и самые длинные перерывы.

use std::path::Path;

use anise::time::{Duration, Epoch};
use serde::{Deserialize, Serialize};

use crate::runlog::{ser_duration_s, ser_epoch};

use crate::error::{Error, Result};
use crate::sky::{AzEl, LunarSite, SkySample};
use crate::terrain::{HorizonMask, HorizonParams, TerrainSet};

/// Радиус Солнца, км (номинальное значение IAU 2015).
pub const SUN_RADIUS_KM: f64 = 695_700.0;
/// Экваториальный радиус Земли, км (WGS84).
pub const EARTH_RADIUS_KM: f64 = 6_378.137;

// ---------------------------------------------------------------------------
// Параметры
// ---------------------------------------------------------------------------

/// Все параметры видимости. Каждое поле имеет значение по умолчанию,
/// поэтому в TOML-файле достаточно указать только то, что меняется.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct VisibilityConfig {
    pub horizon: HorizonConfig,
    pub sun: BodyRule,
    pub earth: BodyRule,
}

impl Default for VisibilityConfig {
    fn default() -> Self {
        Self { horizon: HorizonConfig::default(), sun: BodyRule::default(), earth: BodyRule::default() }
    }
}

/// Как строится маска горизонта.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct HorizonConfig {
    /// Высота наблюдателя (солнечной панели, антенны) над поверхностью, м.
    pub observer_height_m: f64,
    /// Насколько далеко учитывать рельеф, км. Высокие горы видны издалека:
    /// пик 8 км над наблюдателем на высоте 5 км выглядывает из-за горизонта примерно с 300 км.
    pub max_range_km: f64,
    /// Число направлений по азимуту (720 = шаг 0.5 градуса).
    pub azimuths: usize,
}

impl Default for HorizonConfig {
    fn default() -> Self {
        Self { observer_height_m: 2.0, max_range_km: 300.0, azimuths: 720 }
    }
}

/// Какая часть диска должна быть над горизонтом, чтобы тело считалось видимым.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum DiscCriterion {
    /// Центр диска над горизонтом.
    Centre,
    /// Хотя бы верхний край диска над горизонтом.
    AnyPart,
    /// Весь диск над горизонтом.
    FullDisc,
}

/// Правило видимости для одного тела.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct BodyRule {
    pub criterion: DiscCriterion,
    /// Запас над горизонтом, градусы (например, минимальный угол места для антенны).
    /// Отсчитывается от точки диска, заданной `criterion`.
    pub min_clearance_deg: f64,
}

impl Default for BodyRule {
    fn default() -> Self {
        Self { criterion: DiscCriterion::Centre, min_clearance_deg: 0.0 }
    }
}

impl VisibilityConfig {
    pub fn parse(text: &str, origin: &Path) -> Result<Self> {
        let cfg: Self = toml::from_str(text)
            .map_err(|source| Error::Toml { path: origin.to_path_buf(), source })?;
        cfg.validate()?;
        Ok(cfg)
    }

    pub async fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = tokio::fs::read_to_string(path)
            .await
            .map_err(|source| Error::Io { path: path.to_path_buf(), source })?;
        Self::parse(&text, path)
    }

    pub fn validate(&self) -> Result<()> {
        let h = &self.horizon;
        if !(h.observer_height_m >= 0.0 && h.observer_height_m.is_finite()) {
            return Err(Error::Config("horizon.observer_height_m must be >= 0".into()));
        }
        if !(h.max_range_km > 0.0 && h.max_range_km.is_finite()) {
            return Err(Error::Config("horizon.max_range_km must be > 0".into()));
        }
        if h.azimuths < 36 {
            return Err(Error::Config("horizon.azimuths must be >= 36".into()));
        }
        for (name, rule) in [("sun", &self.sun), ("earth", &self.earth)] {
            if !(rule.min_clearance_deg.abs() <= 45.0) {
                return Err(Error::Config(format!("{name}.min_clearance_deg must be within ±45")));
            }
        }
        Ok(())
    }

    fn horizon_params(&self) -> HorizonParams {
        HorizonParams {
            azimuths: self.horizon.azimuths,
            max_range_km: self.horizon.max_range_km,
            observer_height_m: self.horizon.observer_height_m,
        }
    }
}

// ---------------------------------------------------------------------------
// Точка и расчёт
// ---------------------------------------------------------------------------

/// Сведения о рельефе в точке.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SiteTerrain {
    /// Высота поверхности над опорной сферой, м.
    pub surface_height_m: f64,
    /// Уклон, градусы (на базе одного пикселя карты), если его удалось посчитать.
    pub slope_deg: Option<f64>,
    /// Карта, по которой взята точка.
    pub map_name: String,
    /// Размер пикселя этой карты в точке, м.
    pub resolution_m: f64,
    /// Лучей горизонта, которые вышли за край карт раньше `max_range_km`.
    pub truncated_rays: usize,
}

/// Точка на Луне вместе с её горизонтом и правилами видимости.
#[derive(Debug, Clone)]
pub struct SiteVisibility {
    /// Положение наблюдателя (высота уже включает высоту над поверхностью).
    pub site: LunarSite,
    pub config: VisibilityConfig,
    /// `None` — рельеф не учитывается (ровный горизонт).
    pub terrain: Option<SiteTerrain>,
    mask: Option<HorizonMask>,
}

impl SiteVisibility {
    /// Готовит точку. С `terrain` высота берётся из карты и строится маска горизонта;
    /// без него используется `flat_height_m` и ровный горизонт.
    /// CPU-нагруженный вызов: из async — через `spawn_blocking`.
    pub fn new(
        terrain: Option<&TerrainSet>,
        lat_deg: f64,
        lon_deg: f64,
        flat_height_m: f64,
        config: VisibilityConfig,
    ) -> Result<Self> {
        config.validate()?;
        let observer_m = config.horizon.observer_height_m;
        match terrain {
            Some(set) => {
                let here = set.sample(lat_deg, lon_deg).ok_or_else(|| {
                    Error::Config(format!("site ({lat_deg}, {lon_deg}) is not covered by any terrain map"))
                })?;
                let mask = HorizonMask::build(set, lat_deg, lon_deg, config.horizon_params())?;
                let info = SiteTerrain {
                    surface_height_m: here.height_m,
                    slope_deg: set.slope_deg(lat_deg, lon_deg),
                    map_name: here.map_name.to_string(),
                    resolution_m: here.resolution_m,
                    truncated_rays: mask.truncated_rays,
                };
                let site = LunarSite::new(lat_deg, lon_deg, (here.height_m + observer_m) / 1000.0)?;
                Ok(Self { site, config, terrain: Some(info), mask: Some(mask) })
            }
            None => {
                let site = LunarSite::new(lat_deg, lon_deg, (flat_height_m + observer_m) / 1000.0)?;
                Ok(Self { site, config, terrain: None, mask: None })
            }
        }
    }

    /// Маска горизонта: угол возвышения рельефа по азимутам (индекс i — азимут i * 360 / N).
    /// `None` без рельефа.
    pub fn horizon_mask_deg(&self) -> Option<&[f64]> {
        self.mask.as_ref().map(|m| m.elevation_deg.as_slice())
    }

    /// Угол возвышения горизонта в направлении `azimuth_deg`, градусы (0 без рельефа).
    pub fn horizon_deg(&self, azimuth_deg: f64) -> f64 {
        self.mask.as_ref().map_or(0.0, |m| m.at(azimuth_deg))
    }

    /// Видимость Солнца и Земли в один момент.
    pub fn evaluate(&self, sky: &SkySample) -> VisibilitySample {
        VisibilitySample {
            epoch: sky.epoch,
            sun: self.body(&sky.sun, SUN_RADIUS_KM, &self.config.sun),
            earth: self.body(&sky.earth, EARTH_RADIUS_KM, &self.config.earth),
        }
    }

    pub fn evaluate_all(&self, sky: &[SkySample]) -> Vec<VisibilitySample> {
        sky.iter().map(|s| self.evaluate(s)).collect()
    }

    fn body(&self, pos: &AzEl, radius_km: f64, rule: &BodyRule) -> BodyVisibility {
        let horizon_deg = self.horizon_deg(pos.azimuth_deg);
        let angular_radius_deg = (radius_km / pos.range_km).clamp(-1.0, 1.0).asin().to_degrees();
        let above = pos.elevation_deg - horizon_deg;
        let reference = match rule.criterion {
            DiscCriterion::Centre => above,
            DiscCriterion::AnyPart => above + angular_radius_deg,
            DiscCriterion::FullDisc => above - angular_radius_deg,
        };
        BodyVisibility {
            azimuth_deg: pos.azimuth_deg,
            elevation_deg: pos.elevation_deg,
            horizon_deg,
            angular_radius_deg,
            visible_fraction: disc_fraction_above(above, angular_radius_deg),
            visible: reference > rule.min_clearance_deg,
        }
    }
}

/// Доля площади диска с угловым радиусом `radius_deg`, чей центр на `centre_above_deg`
/// выше линии горизонта (горизонт в пределах диска считается прямой линией).
pub fn disc_fraction_above(centre_above_deg: f64, radius_deg: f64) -> f64 {
    if radius_deg <= 0.0 {
        return if centre_above_deg > 0.0 { 1.0 } else { 0.0 };
    }
    let h = (centre_above_deg / radius_deg).clamp(-1.0, 1.0);
    1.0 - (h.acos() - h * (1.0 - h * h).sqrt()) / std::f64::consts::PI
}

// ---------------------------------------------------------------------------
// Результаты
// ---------------------------------------------------------------------------

/// Видимость одного тела в один момент.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct BodyVisibility {
    pub azimuth_deg: f64,
    /// Угол возвышения над ровным горизонтом, градусы.
    pub elevation_deg: f64,
    /// Угол возвышения рельефа в этом направлении, градусы.
    pub horizon_deg: f64,
    /// Видимый угловой радиус диска, градусы.
    pub angular_radius_deg: f64,
    /// Доля диска над рельефом, 0..=1 (для Солнца — доля доступного света).
    pub visible_fraction: f64,
    /// Видимо ли тело по правилу из конфигурации.
    pub visible: bool,
}

/// Видимость Солнца и Земли в один момент.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct VisibilitySample {
    #[serde(serialize_with = "ser_epoch")]
    pub epoch: Epoch,
    pub sun: BodyVisibility,
    pub earth: BodyVisibility,
}

impl VisibilitySample {
    pub fn both(&self) -> bool {
        self.sun.visible && self.earth.visible
    }
}

/// Сводка по ряду с постоянным шагом.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct VisibilitySummary {
    pub samples: usize,
    pub sun_pct: f64,
    pub earth_pct: f64,
    pub both_pct: f64,
    /// Самые длинные перерывы (подряд идущие моменты без видимости) с точностью до шага.
    /// Перерыв на краю окна расчёта может продолжаться за его пределами.
    #[serde(rename = "longest_sun_gap_s", serialize_with = "ser_duration_s")]
    pub longest_sun_gap: Duration,
    #[serde(rename = "longest_earth_gap_s", serialize_with = "ser_duration_s")]
    pub longest_earth_gap: Duration,
    #[serde(rename = "longest_both_gap_s", serialize_with = "ser_duration_s")]
    pub longest_both_gap: Duration,
}

impl VisibilitySummary {
    pub fn from_samples(samples: &[VisibilitySample], step: Duration) -> Self {
        let n = samples.len();
        let pct = |f: &dyn Fn(&VisibilitySample) -> bool| {
            if n == 0 { 0.0 } else { 100.0 * samples.iter().filter(|s| f(s)).count() as f64 / n as f64 }
        };
        let gap = |f: &dyn Fn(&VisibilitySample) -> bool| {
            let (mut best, mut run) = (0usize, 0usize);
            for s in samples {
                run = if f(s) { 0 } else { run + 1 };
                best = best.max(run);
            }
            step * (best as f64)
        };
        Self {
            samples: n,
            sun_pct: pct(&|s| s.sun.visible),
            earth_pct: pct(&|s| s.earth.visible),
            both_pct: pct(&|s| s.both()),
            longest_sun_gap: gap(&|s| s.sun.visible),
            longest_earth_gap: gap(&|s| s.earth.visible),
            longest_both_gap: gap(&|s| s.both()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anise::time::Unit;
    use approx::assert_abs_diff_eq;

    fn azel(az: f64, el: f64, range_km: f64) -> AzEl {
        AzEl { azimuth_deg: az, elevation_deg: el, range_km }
    }

    fn flat(config: VisibilityConfig) -> SiteVisibility {
        SiteVisibility::new(None, 0.0, 0.0, 0.0, config).unwrap()
    }

    fn sample(sun_el: f64, earth_el: f64) -> SkySample {
        SkySample {
            epoch: Epoch::from_gregorian_utc_at_midnight(2026, 11, 14),
            sun: azel(90.0, sun_el, 1.496e8),
            earth: azel(0.0, earth_el, 384_400.0),
        }
    }

    #[test]
    fn disc_fraction_limits_and_half() {
        assert_abs_diff_eq!(disc_fraction_above(0.0, 0.25), 0.5, epsilon = 1e-12);
        assert_abs_diff_eq!(disc_fraction_above(0.25, 0.25), 1.0, epsilon = 1e-12);
        assert_abs_diff_eq!(disc_fraction_above(-0.25, 0.25), 0.0, epsilon = 1e-12);
        assert_abs_diff_eq!(disc_fraction_above(5.0, 0.25), 1.0, epsilon = 1e-12);
        // Симметрия: доля над + доля под = 1.
        let f = disc_fraction_above(0.1, 0.25);
        assert_abs_diff_eq!(f + disc_fraction_above(-0.1, 0.25), 1.0, epsilon = 1e-12);
        assert!(f > 0.5 && f < 1.0);
    }

    #[test]
    fn angular_radii_are_physical() {
        let v = flat(VisibilityConfig::default()).evaluate(&sample(10.0, 10.0));
        assert_abs_diff_eq!(v.sun.angular_radius_deg, 0.2664, epsilon = 0.001);
        assert_abs_diff_eq!(v.earth.angular_radius_deg, 0.9507, epsilon = 0.001);
    }

    #[test]
    fn criteria_differ_near_horizon() {
        // Центр Земли на 0.5 градуса ниже горизонта, радиус диска ~0.95: верхний край виден.
        let s = sample(10.0, -0.5);
        let mut cfg = VisibilityConfig::default();
        cfg.earth.criterion = DiscCriterion::Centre;
        assert!(!flat(cfg.clone()).evaluate(&s).earth.visible);
        cfg.earth.criterion = DiscCriterion::AnyPart;
        assert!(flat(cfg.clone()).evaluate(&s).earth.visible);
        // Центр на 0.5 выше: целиком диск ещё не поднялся.
        let s = sample(10.0, 0.5);
        cfg.earth.criterion = DiscCriterion::FullDisc;
        assert!(!flat(cfg.clone()).evaluate(&s).earth.visible);
        cfg.earth.criterion = DiscCriterion::Centre;
        assert!(flat(cfg).evaluate(&s).earth.visible);
    }

    #[test]
    fn clearance_margin_applies() {
        let mut cfg = VisibilityConfig::default();
        cfg.earth.min_clearance_deg = 5.0;
        let site = flat(cfg);
        assert!(!site.evaluate(&sample(10.0, 4.0)).earth.visible);
        assert!(site.evaluate(&sample(10.0, 6.0)).earth.visible);
    }

    #[test]
    fn observer_height_is_added_to_site() {
        let mut cfg = VisibilityConfig::default();
        cfg.horizon.observer_height_m = 10.0;
        let site = SiteVisibility::new(None, 0.0, 0.0, 500.0, cfg).unwrap();
        assert_abs_diff_eq!(site.site.height_km, 0.51, epsilon = 1e-12);
    }

    #[test]
    fn summary_percentages_and_gaps() {
        let site = flat(VisibilityConfig::default());
        // Солнце: видно, нет, нет, нет, видно; Земля видна всегда.
        let sky: Vec<SkySample> =
            [5.0, -1.0, -1.0, -1.0, 5.0].iter().map(|&e| sample(e, 3.0)).collect();
        let vis = site.evaluate_all(&sky);
        let sum = VisibilitySummary::from_samples(&vis, 1 * Unit::Hour);
        assert_eq!(sum.samples, 5);
        assert_abs_diff_eq!(sum.sun_pct, 40.0, epsilon = 1e-9);
        assert_abs_diff_eq!(sum.earth_pct, 100.0, epsilon = 1e-9);
        assert_abs_diff_eq!(sum.both_pct, 40.0, epsilon = 1e-9);
        assert_eq!(sum.longest_sun_gap, 3 * Unit::Hour);
        assert_eq!(sum.longest_earth_gap, Duration::ZERO);
    }

    #[test]
    fn config_defaults_partial_file_and_unknown_keys() {
        let cfg = VisibilityConfig::parse("", Path::new("v.toml")).unwrap();
        assert_eq!(cfg, VisibilityConfig::default());
        let cfg = VisibilityConfig::parse(
            "[earth]\ncriterion = \"any_part\"\nmin_clearance_deg = 3.0\n",
            Path::new("v.toml"),
        )
        .unwrap();
        assert_eq!(cfg.earth.criterion, DiscCriterion::AnyPart);
        assert_eq!(cfg.sun, BodyRule::default());
        assert_eq!(cfg.horizon, HorizonConfig::default());
        // Опечатка в имени параметра — ошибка, а не тихое значение по умолчанию.
        assert!(VisibilityConfig::parse("[horizon]\nobserver_heigth_m = 3\n", Path::new("v.toml")).is_err());
        assert!(VisibilityConfig::parse("[horizon]\nmax_range_km = -1\n", Path::new("v.toml")).is_err());
    }

    #[test]
    fn bundled_config_file_is_valid() {
        let text = include_str!("../data/visibility.toml");
        VisibilityConfig::parse(text, Path::new("data/visibility.toml")).unwrap();
    }
}
