//! Положение Солнца и Земли на небе точки на поверхности Луны.
//!
//! Здесь только геометрия: азимут, угол возвышения над ровным горизонтом и дальность.
//! Видна ли цель с учётом рельефа и размера диска, решает [`crate::visibility`].
//!
//! Луна считается сферой радиуса [`MOON_RADIUS_KM`] (опорная сфера LOLA); высота точки
//! задаётся в [`LunarSite::height_km`].

use anise::constants::frames::{EARTH_J2000, MOON_ME_DE440_ME421_FRAME, SUN_J2000};
use anise::math::Vector3;
use anise::prelude::{Aberration, Almanac};
use anise::time::{Duration, Epoch};
use rayon::prelude::*;
use serde::Serialize;

use crate::error::{Error, Result};

/// Радиус опорной сферы Луны, км (тот же, что у высот LOLA).
pub const MOON_RADIUS_KM: f64 = 1737.4;

/// Точка на Луне в системе Mean Earth / Polar Axis (в ней заданы карты LOLA).
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct LunarSite {
    /// Широта, градусы (+ север).
    pub lat_deg: f64,
    /// Долгота, градусы (+ восток).
    pub lon_deg: f64,
    /// Высота над опорной сферой, км.
    pub height_km: f64,
}

impl LunarSite {
    pub fn new(lat_deg: f64, lon_deg: f64, height_km: f64) -> Result<Self> {
        if !(-90.0..=90.0).contains(&lat_deg) {
            return Err(Error::Config(format!("site latitude {lat_deg} out of range [-90, 90]")));
        }
        if !(-180.0..=360.0).contains(&lon_deg) {
            return Err(Error::Config(format!("site longitude {lon_deg} out of range [-180, 360]")));
        }
        if !height_km.is_finite() || height_km <= -MOON_RADIUS_KM {
            return Err(Error::Config(format!("site height {height_km} km is invalid")));
        }
        Ok(Self { lat_deg, lon_deg, height_km })
    }

    /// Единичный вектор «вверх» (от центра Луны к точке).
    pub fn up(&self) -> Vector3 {
        let (lat, lon) = (self.lat_deg.to_radians(), self.lon_deg.to_radians());
        Vector3::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin())
    }

    /// Радиус-вектор точки в связанной с Луной системе, км.
    pub fn position_km(&self) -> Vector3 {
        self.up() * (MOON_RADIUS_KM + self.height_km)
    }
}

/// Азимут, угол возвышения и дальность цели из точки на поверхности.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct AzEl {
    /// Азимут, градусы [0, 360): 0 = север, 90 = восток.
    /// Точно на полюсе «севера» нет: там азимут отсчитывается от направления меридиана 0°.
    pub azimuth_deg: f64,
    /// Угол возвышения над ровным горизонтом, градусы [-90, 90].
    pub elevation_deg: f64,
    pub range_km: f64,
}

/// Локальные оси (вверх, север, восток) для единичного вектора `up`.
/// Север — проекция оси вращения на плоскость горизонта; точно на полюсе
/// вместо него берётся направление меридиана 0°.
pub fn local_frame(up: &Vector3) -> (Vector3, Vector3, Vector3) {
    let z = Vector3::new(0.0, 0.0, 1.0);
    let mut north = z - up * z.dot(up);
    if north.norm() < 1e-9 {
        let x = Vector3::new(1.0, 0.0, 0.0);
        north = x - up * x.dot(up);
    }
    let north = north.normalize();
    let east = north.cross(up);
    (*up, north, east)
}

/// Чистая геометрия: `target_km` — положение цели относительно центра Луны
/// в той же связанной системе, что и точка.
pub fn az_el(site: &LunarSite, target_km: Vector3) -> AzEl {
    let (up, north, east) = local_frame(&site.up());
    let los = target_km - site.position_km();
    let range_km = los.norm();
    let dir = los / range_km;

    let elevation_deg = dir.dot(&up).clamp(-1.0, 1.0).asin().to_degrees();

    let azimuth_deg = dir.dot(&east).atan2(dir.dot(&north)).to_degrees().rem_euclid(360.0);

    AzEl { azimuth_deg, elevation_deg, range_km }
}

/// Небо над точкой в один момент времени.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SkySample {
    pub epoch: Epoch,
    pub sun: AzEl,
    pub earth: AzEl,
}

/// Положение Солнца и Земли над точкой в момент `epoch`.
/// Используется видимое положение (поправка за световое время и аберрацию, LT+S).
pub fn sky_at(almanac: &Almanac, site: &LunarSite, epoch: Epoch) -> Result<SkySample> {
    let body = |target, name: &str| -> Result<AzEl> {
        let state = almanac
            .transform(target, MOON_ME_DE440_ME421_FRAME, epoch, Aberration::LT_S)
            .map_err(|e| Error::Almanac(format!("{name} at {epoch}: {e}")))?;
        Ok(az_el(site, state.radius_km))
    };
    Ok(SkySample {
        epoch,
        sun: body(SUN_J2000, "Sun")?,
        earth: body(EARTH_J2000, "Earth")?,
    })
}

/// Ряд значений от `start` до `end` включительно с шагом `step` (считается параллельно).
pub fn sky_series(
    almanac: &Almanac,
    site: &LunarSite,
    start: Epoch,
    end: Epoch,
    step: Duration,
) -> Result<Vec<SkySample>> {
    if step <= Duration::ZERO {
        return Err(Error::Config("time step must be positive".into()));
    }
    if end < start {
        return Err(Error::Config("end epoch is before start epoch".into()));
    }
    let steps = ((end - start).to_seconds() / step.to_seconds()).floor() as u64;
    (0..=steps)
        .into_par_iter()
        .map(|i| sky_at(almanac, site, start + step * (i as f64)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    const FAR: f64 = 1.0e9;

    fn site(lat: f64, lon: f64) -> LunarSite {
        LunarSite::new(lat, lon, 0.0).unwrap()
    }

    #[test]
    fn site_position_on_reference_sphere() {
        let p = site(0.0, 90.0).position_km();
        assert_abs_diff_eq!(p.x, 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(p.y, MOON_RADIUS_KM, epsilon = 1e-9);
        let pole = site(-90.0, 0.0).position_km();
        assert_abs_diff_eq!(pole.z, -MOON_RADIUS_KM, epsilon = 1e-9);
    }

    #[test]
    fn zenith_and_nadir() {
        let s = site(0.0, 0.0);
        let up = az_el(&s, Vector3::new(FAR, 0.0, 0.0));
        assert_abs_diff_eq!(up.elevation_deg, 90.0, epsilon = 1e-6);
        let down = az_el(&s, Vector3::new(-FAR, 0.0, 0.0));
        assert_abs_diff_eq!(down.elevation_deg, -90.0, epsilon = 1e-6);
    }

    #[test]
    fn azimuth_north_east_south_west() {
        let s = site(0.0, 0.0);
        let north = az_el(&s, Vector3::new(0.0, 0.0, FAR));
        assert_abs_diff_eq!(north.elevation_deg, 0.0, epsilon = 1e-3);
        assert!(north.azimuth_deg < 1e-3 || north.azimuth_deg > 360.0 - 1e-3);
        assert_abs_diff_eq!(az_el(&s, Vector3::new(0.0, FAR, 0.0)).azimuth_deg, 90.0, epsilon = 1e-3);
        assert_abs_diff_eq!(az_el(&s, Vector3::new(0.0, 0.0, -FAR)).azimuth_deg, 180.0, epsilon = 1e-3);
        assert_abs_diff_eq!(az_el(&s, Vector3::new(0.0, -FAR, 0.0)).azimuth_deg, 270.0, epsilon = 1e-3);
    }

    #[test]
    fn range_from_site_not_from_centre() {
        let s = site(0.0, 0.0);
        let r = az_el(&s, Vector3::new(MOON_RADIUS_KM + 100.0, 0.0, 0.0));
        assert_abs_diff_eq!(r.range_km, 100.0, epsilon = 1e-9);
    }

    #[test]
    fn pole_has_finite_azimuth() {
        let s = site(-90.0, 0.0);
        let a = az_el(&s, Vector3::new(FAR, 0.0, 0.0));
        assert!(a.azimuth_deg.is_finite());
        assert_abs_diff_eq!(a.elevation_deg, 0.0, epsilon = 1e-3);
        // Направление меридиана 0° на полюсе принято за начало отсчёта азимута.
        assert!(a.azimuth_deg < 1e-3 || a.azimuth_deg > 360.0 - 1e-3);
    }

    #[test]
    fn rejects_invalid_site() {
        assert!(LunarSite::new(91.0, 0.0, 0.0).is_err());
        assert!(LunarSite::new(0.0, 400.0, 0.0).is_err());
        assert!(LunarSite::new(0.0, 0.0, f64::NAN).is_err());
    }

    #[test]
    fn cli_epoch_format_parses() {
        use std::str::FromStr;
        let e = Epoch::from_str("2026-11-14T12:00:00 UTC").unwrap();
        assert_eq!(e, Epoch::from_gregorian_utc(2026, 11, 14, 12, 0, 0, 0));
    }

    #[test]
    fn series_rejects_bad_arguments() {
        use anise::time::Unit;
        let almanac = Almanac::default();
        let s = site(0.0, 0.0);
        let t0 = Epoch::from_gregorian_utc_at_midnight(2026, 11, 14);
        assert!(sky_series(&almanac, &s, t0, t0 + 1 * Unit::Hour, Duration::ZERO).is_err());
        assert!(sky_series(&almanac, &s, t0, t0 - 1 * Unit::Hour, 1 * Unit::Hour).is_err());
        // Без ядер расчёт должен вернуть ошибку, а не выдуманные числа.
        assert!(sky_at(&almanac, &s, t0).is_err());
    }
}
