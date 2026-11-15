//! Прямой перелёт (Direct Transfer): задача Ламберта в поле Земли.
//!
//! Для каждого прохода космодрома через плоскость перебирается время на опорной орбите
//! (с шагом `coast_step`); из точки разгона решается задача Ламберта до положения Луны в момент
//! прилёта (DE440) при любом угле перелёта. Время полёта — в диапазоне
//! `tof_min_days..tof_max_days`; из вариантов выбирается лучший по критерию.

use serde::{Deserialize, Serialize};

use super::common::{CrossingGeometry, Leg, LegCounters, LegModel, Transfer, TransferCommonConfig};
use crate::error::{Error, Result};
use crate::lambert;

/// Прямой перелёт.
pub type DirectTransfer = Transfer<DirectTransferConfig>;

/// Свои параметры прямого перелёта (секция `[direct_transfer]`).
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct DirectTransferConfig {
    /// Допустимое время перелёта от разгона до Луны, сутки.
    pub tof_min_days: f64,
    pub tof_max_days: f64,
}

impl Default for DirectTransferConfig {
    fn default() -> Self {
        Self { tof_min_days: 2.5, tof_max_days: 5.5 }
    }
}

impl DirectTransferConfig {
    pub fn validate(&self, _common: &TransferCommonConfig) -> Result<()> {
        if !(self.tof_min_days >= 1.0 && self.tof_max_days > self.tof_min_days && self.tof_max_days <= 10.0) {
            return Err(Error::Config("direct_transfer: need 1 <= tof_min_days < tof_max_days <= 10".into()));
        }
        Ok(())
    }
}

impl LegModel for DirectTransferConfig {
    const NAME: &'static str = "direct_transfer";
    const SECTION: &'static str = "direct_transfer";
    const DESCRIPTION: &'static str =
        "direct transfer (Lambert arc at any transfer angle), patched conic (Earth gravity only on the way, Moon ephemeris DE440)";

    fn validate(&self, common: &TransferCommonConfig) -> Result<()> {
        DirectTransferConfig::validate(self, common)
    }

    fn time_of_flight_bounds_s(&self, _common: &TransferCommonConfig) -> (f64, f64) {
        (self.tof_min_days * 86_400.0, self.tof_max_days * 86_400.0)
    }

    fn legs(&self, g: &CrossingGeometry, counters: &mut LegCounters, out: &mut dyn FnMut(Leg)) {
        let step_s = g.coast_step_s;
        let k_min = (g.coast_min_s / step_s).ceil() as i64;
        let k_max = (g.coast_max_s / step_s).floor() as i64;
        // Внутренний цикл — в секундах (f64): арифметика дат здесь была дороже задачи Ламберта.
        let launch_to_arrival_s = g.launch_to_arrival_s();
        let (tof_min_s, tof_max_s) = (self.tof_min_days * 86_400.0, self.tof_max_days * 86_400.0);
        for k in k_min..=k_max {
            let coast_s = step_s * k as f64;
            if coast_s < g.ascent_s + 60.0 {
                continue;
            }
            let tof_s = launch_to_arrival_s - coast_s;
            if tof_s < tof_min_s || tof_s > tof_max_s {
                continue;
            }
            let ang = g.ascent_phase_rad + g.n_park_rad_s * coast_s;
            let r1_hat = g.parking_direction(ang);
            let r1 = r1_hat * g.parking_radius_km;
            counters.solves += 1;
            let Ok(sol) = lambert::solve(&r1, &g.moon_position_km, tof_s, g.mu_earth_km3_s2, &g.plane_normal) else {
                counters.failures += 1;
                continue;
            };
            out(Leg {
                coast_s,
                coast_angle_rad: ang,
                r1,
                v_park: g.plane_normal.cross(&r1_hat) * g.circular_speed_km_s,
                v1: sol.v1,
                v2: sol.v2,
                transfer_angle_rad: sol.transfer_angle_rad,
                hohmann: None,
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::*;
    use super::*;
    use crate::ephemeris::{AlmanacSource, load_almanac_blocking};
    use crate::spaceport::SpaceportCatalog;
    use anise::prelude::Almanac;
    use anise::time::{Epoch, Unit};
    use approx::assert_abs_diff_eq;
    use std::path::Path;

    #[test]
    fn config_validation() {
        let common = TransferCommonConfig::default();
        assert!(DirectTransferConfig::default().validate(&common).is_ok());
        assert!(DirectTransferConfig { tof_min_days: 6.0, tof_max_days: 5.5 }.validate(&common).is_err());
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
        assert!(fast.departure.delta_v_magnitude_km_s <= TransferCommonConfig::default().max_tli_delta_v_km_s);

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
