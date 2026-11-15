//! Проверка пригодности площадки для посадки по карте рельефа.
//!
//! Уклон, который видит аппарат, зависит от масштаба: на грубой карте склон небольшого
//! кратера сглаживается. Поэтому вердикт честно учитывает детализацию карты:
//!
//! - [`Verdict::Unsuitable`] — уклон превышает допустимый уже на этой карте. На более
//!   подробной карте он обычно только круче, поэтому «нельзя» можно утверждать и по грубой.
//! - [`Verdict::NotVerified`] — на этой карте нарушений не видно, но её пиксель крупнее
//!   [`LandingConfig::max_verified_resolution_m`], на масштабе опор аппарата не проверено.
//! - [`Verdict::Suitable`] — нарушений нет, и карта достаточно подробная.
//! - [`Verdict::NotChecked`] — проверить нельзя (часть области не покрыта картами).
//!
//! Валуны и кратеры меньше пикселя карта высот не показывает: их здесь нет ни в одном вердикте.
//!
//! Параметры — в [`LandingConfig`], читаются из TOML (по умолчанию `data/landing.toml`).

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};
use crate::sky::{LunarSite, local_frame};
use crate::terrain::{TerrainSet, offset_lat_lon};

// ---------------------------------------------------------------------------
// Параметры
// ---------------------------------------------------------------------------

/// Требования аппарата к площадке.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct LandingConfig {
    /// Максимально допустимый уклон, градусы.
    pub max_slope_deg: f64,
    /// Размах посадочных опор аппарата, м. На этой базе (но не мельче пикселя карты) меряется уклон.
    pub footprint_m: f64,
    /// Радиус проверяемой области вокруг точки, м (учитывает точность посадки).
    pub check_radius_m: f64,
    /// Карта с пикселем не крупнее этого значения считается достаточной для вердикта «пригодно», м.
    pub max_verified_resolution_m: f64,
    /// Допустимая доля области с уклоном выше `max_slope_deg`, 0..1.
    /// 0 — вся область должна быть пологой (строго); например, 0.05 — допускается 5 % крутых точек
    /// (имеет смысл, если аппарат может уйти от препятствия при посадке).
    pub max_steep_fraction: f64,
}

impl Default for LandingConfig {
    fn default() -> Self {
        Self {
            // Требование к площадкам Artemis в материалах NASA: уклон меньше 10 градусов.
            max_slope_deg: 10.0,
            // Допущение: типичный размах опор посадочного аппарата.
            footprint_m: 5.0,
            // Точность посадки пилотируемого аппарата Artemis: 100 x 100 м.
            check_radius_m: 50.0,
            max_verified_resolution_m: 10.0,
            max_steep_fraction: 0.0,
        }
    }
}

impl LandingConfig {
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
        let positive = |v: f64| v > 0.0 && v.is_finite();
        if !(self.max_slope_deg > 0.0 && self.max_slope_deg < 90.0) {
            return Err(Error::Config("max_slope_deg must be in (0, 90)".into()));
        }
        if !positive(self.footprint_m) {
            return Err(Error::Config("footprint_m must be > 0".into()));
        }
        if !(self.check_radius_m >= 0.0 && self.check_radius_m.is_finite()) {
            return Err(Error::Config("check_radius_m must be >= 0".into()));
        }
        if self.check_radius_m > 5_000.0 {
            return Err(Error::Config("check_radius_m must be <= 5000".into()));
        }
        if !positive(self.max_verified_resolution_m) {
            return Err(Error::Config("max_verified_resolution_m must be > 0".into()));
        }
        if !(0.0..1.0).contains(&self.max_steep_fraction) {
            return Err(Error::Config("max_steep_fraction must be in [0, 1)".into()));
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Результат
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(tag = "verdict", rename_all = "snake_case")]
pub enum Verdict {
    /// Сесть нельзя: уклон больше допустимого уже на доступной карте.
    Unsuitable,
    /// Нарушений не видно, но карта слишком грубая, чтобы утверждать «можно».
    NotVerified,
    /// Пригодно по уклону на масштабе опор (валуны не проверялись).
    Suitable,
    /// Проверить нельзя.
    NotChecked { reason: String },
}

/// Одна проверенная точка площадки.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LandingPoint {
    /// Смещение от центра площадки на восток и на север, м.
    pub east_m: f64,
    pub north_m: f64,
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub height_m: f64,
    pub slope_deg: f64,
    pub baseline_m: f64,
    pub resolution_m: f64,
    pub map_name: String,
    pub steep: bool,
}

/// Результат проверки площадки.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LandingAssessment {
    pub verdict: Verdict,
    /// Карта в центре площадки (`None`, если точка не покрыта).
    pub map_name: Option<String>,
    /// Самый крупный пиксель среди использованных карт, м.
    pub resolution_m: f64,
    /// База, на которой мерился уклон (самая крупная в области), м.
    pub baseline_m: f64,
    /// Число проверенных точек в области.
    pub points: usize,
    pub max_slope_deg: f64,
    pub mean_slope_deg: f64,
    /// Доля точек области с уклоном выше допустимого, 0..1.
    pub steep_fraction: f64,
    /// Самая крутая точка области (широта, долгота).
    pub steepest_at: Option<(f64, f64)>,
    /// Перепад высот в области, м.
    pub relief_m: f64,
    /// Все проверенные точки.
    pub samples: Vec<LandingPoint>,
}

impl LandingAssessment {
    /// Короткое объяснение вердикта для человека.
    pub fn explain(&self, cfg: &LandingConfig) -> String {
        match &self.verdict {
            Verdict::Unsuitable => format!(
                "UNSUITABLE: {:.0}% of the area is steeper than {:.1} deg (allowed {:.0}%), \
                 max {:.1} deg, already at {:.0} m/pix",
                self.steep_fraction * 100.0,
                cfg.max_slope_deg,
                cfg.max_steep_fraction * 100.0,
                self.max_slope_deg,
                self.resolution_m
            ),
            Verdict::NotVerified => format!(
                "NOT VERIFIED: slope within limits on this map ({:.0}% above {:.1} deg), but the map \
                 is {:.0} m/pix (need <= {:.0} m/pix to judge at lander scale)",
                self.steep_fraction * 100.0,
                cfg.max_slope_deg,
                self.resolution_m,
                cfg.max_verified_resolution_m
            ),
            Verdict::Suitable => format!(
                "SUITABLE by slope: {:.0}% above {:.1} deg (allowed {:.0}%), max {:.1} deg, \
                 at {:.0} m/pix (boulders not checked)",
                self.steep_fraction * 100.0,
                cfg.max_slope_deg,
                cfg.max_steep_fraction * 100.0,
                self.max_slope_deg,
                self.resolution_m
            ),
            Verdict::NotChecked { reason } => format!("NOT CHECKED: {reason}"),
        }
    }

    fn not_checked(reason: String) -> Self {
        Self {
            verdict: Verdict::NotChecked { reason },
            map_name: None,
            resolution_m: f64::NAN,
            baseline_m: f64::NAN,
            points: 0,
            max_slope_deg: f64::NAN,
            mean_slope_deg: f64::NAN,
            steep_fraction: f64::NAN,
            steepest_at: None,
            relief_m: f64::NAN,
            samples: Vec::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// Расчёт
// ---------------------------------------------------------------------------

/// Проверяет площадку радиусом `cfg.check_radius_m` вокруг точки.
/// CPU-нагруженный вызов: из async — через `spawn_blocking`.
pub fn assess(terrain: &TerrainSet, lat_deg: f64, lon_deg: f64, cfg: &LandingConfig) -> Result<LandingAssessment> {
    cfg.validate()?;
    let centre = LunarSite::new(lat_deg, lon_deg, 0.0)?;
    let Some(here) = terrain.sample(lat_deg, lon_deg) else {
        return Ok(LandingAssessment::not_checked(format!(
            "point ({lat_deg}, {lon_deg}) is not covered by any terrain map"
        )));
    };
    let map_name = here.map_name.to_string();

    // Сетка точек в круге: шаг не мельче половины размаха опор и не мельче пикселя.
    let half_baseline = cfg.footprint_m / 2.0;
    let spacing = half_baseline.max(here.resolution_m);
    let n = (cfg.check_radius_m / spacing).floor() as i64;
    let (up, north, east) = local_frame(&centre.up());

    let mut slopes = Vec::new();
    let mut samples = Vec::new();
    let (mut worst_res, mut worst_base) = (0.0_f64, 0.0_f64);
    let (mut h_min, mut h_max) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut steepest: Option<(f64, (f64, f64))> = None;

    for i in -n..=n {
        for j in -n..=n {
            let (x, y) = (i as f64 * spacing, j as f64 * spacing);
            let d = x.hypot(y);
            if d > cfg.check_radius_m + 1e-9 {
                continue;
            }
            let (lat, lon) = if d == 0.0 {
                (lat_deg, lon_deg)
            } else {
                let dir = (north * y + east * x) / d;
                offset_lat_lon(&up, &dir, d)
            };
            let Some(s) = terrain.slope_over(lat, lon, half_baseline) else {
                return Ok(LandingAssessment::not_checked(format!(
                    "part of the {:.0} m area around the point is not covered by the terrain maps",
                    cfg.check_radius_m
                )));
            };
            worst_res = worst_res.max(s.resolution_m);
            worst_base = worst_base.max(s.baseline_m);
            h_min = h_min.min(s.height_m);
            h_max = h_max.max(s.height_m);
            if steepest.is_none_or(|(best, _)| s.slope_deg > best) {
                steepest = Some((s.slope_deg, (lat, lon)));
            }
            slopes.push(s.slope_deg);
            samples.push(LandingPoint {
                east_m: x,
                north_m: y,
                lat_deg: lat,
                lon_deg: lon,
                height_m: s.height_m,
                slope_deg: s.slope_deg,
                baseline_m: s.baseline_m,
                resolution_m: s.resolution_m,
                map_name: s.map_name.to_string(),
                steep: s.slope_deg > cfg.max_slope_deg,
            });
        }
    }

    let max_slope_deg = steepest.map_or(f64::NAN, |s| s.0);
    let steep = slopes.iter().filter(|&&s| s > cfg.max_slope_deg).count();
    let steep_fraction = steep as f64 / slopes.len() as f64;
    let verdict = if steep_fraction > cfg.max_steep_fraction {
        Verdict::Unsuitable
    } else if worst_res > cfg.max_verified_resolution_m {
        Verdict::NotVerified
    } else {
        Verdict::Suitable
    };

    Ok(LandingAssessment {
        verdict,
        map_name: Some(map_name),
        resolution_m: worst_res,
        baseline_m: worst_base,
        points: slopes.len(),
        max_slope_deg,
        mean_slope_deg: slopes.iter().sum::<f64>() / slopes.len() as f64,
        steep_fraction,
        steepest_at: steepest.map(|s| s.1),
        relief_m: h_max - h_min,
        samples,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::terrain::test_util::{cyl_label, polar_scaled, write_map};
    use approx::assert_abs_diff_eq;

    fn cfg() -> LandingConfig {
        LandingConfig::default()
    }

    #[test]
    fn flat_detailed_map_is_suitable() {
        let set = TerrainSet::new(vec![polar_scaled("land-flat", 200, 5.0, |_, _| 100)]);
        let a = assess(&set, -90.0, 0.0, &cfg()).unwrap();
        assert_eq!(a.verdict, Verdict::Suitable);
        assert!(a.points > 100, "points {}", a.points);
        assert_abs_diff_eq!(a.max_slope_deg, 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(a.relief_m, 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(a.baseline_m, 10.0, epsilon = 1e-9);
    }

    #[test]
    fn steep_ramp_is_unsuitable_even_on_detailed_map() {
        // Подъём 1 м на пиксель 5 м вдоль sample: уклон atan(0.2) = 11.3 градуса > 10.
        let set = TerrainSet::new(vec![polar_scaled("land-ramp", 200, 5.0, |s, _| (s as i16) * 2)]);
        let a = assess(&set, -90.0, 0.0, &cfg()).unwrap();
        assert_eq!(a.verdict, Verdict::Unsuitable);
        assert_abs_diff_eq!(a.max_slope_deg, 11.31, epsilon = 0.01);
    }

    #[test]
    fn gentle_ramp_is_suitable() {
        // atan(0.1) = 5.7 градуса.
        let set = TerrainSet::new(vec![polar_scaled("land-gentle", 200, 5.0, |s, _| s as i16)]);
        let a = assess(&set, -90.0, 0.0, &cfg()).unwrap();
        assert_eq!(a.verdict, Verdict::Suitable);
        assert_abs_diff_eq!(a.max_slope_deg, 5.71, epsilon = 0.01);
    }

    #[test]
    fn steep_fraction_tolerance() {
        // Ровная площадка с одной крутой «стенкой» в нескольких столбцах карты.
        let set = TerrainSet::new(vec![polar_scaled("land-wall", 200, 5.0, |s, _| {
            if s >= 106 { 40 } else { 0 }
        })]);
        let strict = assess(&set, -90.0, 0.0, &cfg()).unwrap();
        assert_eq!(strict.verdict, Verdict::Unsuitable);
        assert!(strict.steep_fraction > 0.0 && strict.steep_fraction < 0.3, "{}", strict.steep_fraction);
        let tolerant = LandingConfig { max_steep_fraction: 0.3, ..cfg() };
        assert_eq!(assess(&set, -90.0, 0.0, &tolerant).unwrap().verdict, Verdict::Suitable);
    }

    #[test]
    fn flat_coarse_map_is_not_verified() {
        let set = TerrainSet::new(vec![write_map("land-coarse", &cyl_label(2), 360, 720, |_, _| 0)]);
        let a = assess(&set, 10.0, 20.0, &cfg()).unwrap();
        assert_eq!(a.verdict, Verdict::NotVerified);
        assert!(a.resolution_m > 1000.0);
        // Пиксель больше области — проверяется одна точка.
        assert_eq!(a.points, 1);
    }

    #[test]
    fn steep_coarse_map_is_unsuitable() {
        // 20 м на пиксель 80 м: atan(0.25) = 14 градусов уже на грубой карте.
        let set = TerrainSet::new(vec![polar_scaled("land-coarse-ramp", 200, 80.0, |s, _| (s as i16) * 40)]);
        let a = assess(&set, -90.0, 0.0, &cfg()).unwrap();
        assert_eq!(a.verdict, Verdict::Unsuitable);
    }

    #[test]
    fn area_leaving_the_map_is_not_checked() {
        // Карта 20 x 20 пикселей по 5 м (100 м): круг радиусом 50 м у полюса выходит за край.
        let set = TerrainSet::new(vec![polar_scaled("land-small", 20, 5.0, |_, _| 0)]);
        let a = assess(&set, -90.0, 0.0, &cfg()).unwrap();
        assert!(matches!(a.verdict, Verdict::NotChecked { .. }));
        let b = assess(&set, 0.0, 0.0, &cfg()).unwrap();
        assert!(matches!(b.verdict, Verdict::NotChecked { .. }));
    }

    #[test]
    fn config_defaults_and_validation() {
        assert_eq!(LandingConfig::parse("", Path::new("l.toml")).unwrap(), LandingConfig::default());
        let c = LandingConfig::parse("max_slope_deg = 15.0\n", Path::new("l.toml")).unwrap();
        assert_eq!(c.max_slope_deg, 15.0);
        assert_eq!(c.footprint_m, 5.0);
        assert!(LandingConfig::parse("max_slope = 15.0\n", Path::new("l.toml")).is_err());
        assert!(LandingConfig::parse("footprint_m = 0\n", Path::new("l.toml")).is_err());
        assert!(LandingConfig::parse("max_steep_fraction = 1.0\n", Path::new("l.toml")).is_err());
    }

    #[test]
    fn bundled_config_file_is_valid() {
        let text = include_str!("../data/landing.toml");
        LandingConfig::parse(text, Path::new("data/landing.toml")).unwrap();
    }

    // ---- Настоящие карты из data/terrain: `cargo test -- --ignored` ----

    #[test]
    #[ignore = "needs LOLA maps in data/terrain"]
    fn real_sites() {
        let t = TerrainSet::open_dir("data/terrain").unwrap();
        // Центральная горка кратера Тихо — крутая уже на глобальной карте.
        let tycho = assess(&t, -43.3, -11.36, &cfg()).unwrap();
        assert_eq!(tycho.verdict, Verdict::Unsuitable, "{tycho:?}");
        // Море Спокойствия (Apollo 11): ровно, но глобальная карта слишком грубая.
        let a11 = assess(&t, 0.674, 23.473, &cfg()).unwrap();
        assert_eq!(a11.verdict, Verdict::NotVerified, "{a11:?}");
        // У южного полюса есть карта 10 м: вердикт выносится на масштабе опор.
        let pole = assess(&t, -89.45, 222.8, &cfg()).unwrap();
        assert!(pole.resolution_m <= 10.0, "{pole:?}");
        assert!(matches!(pole.verdict, Verdict::Suitable | Verdict::Unsuitable), "{pole:?}");
    }
}
