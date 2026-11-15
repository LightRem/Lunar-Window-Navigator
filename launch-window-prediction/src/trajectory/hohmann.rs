//! Перелёт Гомана (Hohmann transfer) к реальной Луне.
//!
//! Классический перелёт Гомана — полуэллипс между двумя круговыми орбитами: разгон по касательной
//! в перигее (на опорной орбите), прилёт в апогее, угол перелёта 180°, время полёта — половина
//! периода эллипса. Здесь вторая «орбита» — положение Луны по эфемеридам DE440:
//!
//! 1. Для момента прилёта `t_a` точка разгона — на опорной орбите строго напротив Луны:
//!    `r1 = −r̂_Луны(t_a) · r_p` (плоскость опорной орбиты и так проходит через Луну).
//! 2. Разгон по касательной (угол наклона траектории 0°), апогей — в направлении Луны.
//! 3. Время полёта задано геометрией: `T = π·√(a³/μ)`, `a = (r_p + |r_Луны|)/2` (~5 суток).
//!    Момент разгона — это момент, когда аппарат после старта проходит точку напротив Луны
//!    (каждый виток опорной орбиты), поэтому время полёта совпадает с `T` только с точностью
//!    до допуска `timing_tolerance_min`. Расхождение компенсируется большой полуосью: эллипс
//!    подбирается так, чтобы за фактическое время полёта прийти в апогей ровно в `t_a`;
//!    апогей при этом отличается от расстояния до Луны на `apogee_miss_km`
//!    (при допуске 2 мин — до ~70 км, после доводки — обычно единицы километров).
//! 4. Лучший вариант доводится: момент прилёта сдвигается, пока расхождение времени не станет
//!    меньше секунды (или перестанет уменьшаться).
//!
//! Скорость у Луны — скорость в апогее по касательной; торможение на окололунную орбиту, посадка,
//! бюджет аппарата и численное уточнение — общие для всех траекторий ([`super::common`]).
//! Численное уточнение ([`crate::refine`]) подбирает скорость после разгона с учётом притяжения
//! Луны и Солнца, поэтому уточнённый разгон немного отличается от чисто касательного.

use std::f64::consts::{PI, TAU};

use serde::{Deserialize, Serialize};

use super::common::{CrossingGeometry, Leg, LegCounters, LegModel, MU_EARTH_DEFAULT, Transfer, TransferCommonConfig};
use crate::error::{Error, Result};

/// Перелёт Гомана.
pub type HohmannTransfer = Transfer<HohmannConfig>;

/// Наименьшее и наибольшее расстояние до Луны, км (с запасом): для окон поиска.
const MOON_DISTANCE_MIN_KM: f64 = 352_000.0;
const MOON_DISTANCE_MAX_KM: f64 = 410_000.0;
/// Доводка момента прилёта продолжается, пока расхождение времени больше этого, с.
const POLISH_TOLERANCE_S: f64 = 1.0;

/// Свои параметры перелёта Гомана (секция `[hohmann]`).
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct HohmannConfig {
    /// Допустимое расхождение времени полёта с временем полёта по эллипсу Гомана
    /// (апогей на расстоянии Луны), минуты. По умолчанию 2.
    ///
    /// В грубом переборе моментов прилёта допуск автоматически не меньше половины шага перебора.
    pub timing_tolerance_min: f64,
}

impl Default for HohmannConfig {
    fn default() -> Self {
        Self { timing_tolerance_min: 2.0 }
    }
}

impl HohmannConfig {
    pub fn validate(&self, _common: &TransferCommonConfig) -> Result<()> {
        if !(self.timing_tolerance_min >= 0.1 && self.timing_tolerance_min <= 60.0) {
            return Err(Error::Config("hohmann: timing_tolerance_min must be in [0.1, 60]".into()));
        }
        Ok(())
    }

    /// Допуск по времени для шага перебора моментов прилёта `arrival_step_s`, с.
    fn tolerance_s(&self, arrival_step_s: f64) -> f64 {
        (self.timing_tolerance_min * 60.0).max(arrival_step_s / 2.0 + 60.0)
    }
}

/// Время полёта по эллипсу Гомана (половина периода) между радиусами `r_p` и `r_a`, с.
pub fn hohmann_time_of_flight_s(r_p: f64, r_a: f64, mu: f64) -> f64 {
    let a = 0.5 * (r_p + r_a);
    PI * (a.powi(3) / mu).sqrt()
}

/// Данные перелёта Гомана для найденного варианта.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct HohmannInfo {
    /// Время полёта по эллипсу Гомана с апогеем на расстоянии Луны в момент прилёта, с.
    pub nominal_time_of_flight_s: f64,
    /// Фактическое время полёта минус номинальное, с (в пределах допуска).
    pub timing_error_s: f64,
    /// Допуск по времени, с которым вариант принят, с.
    pub timing_tolerance_s: f64,
    /// Большая полуось эллипса перелёта, км.
    pub semi_major_axis_km: f64,
    /// Радиус перигея (опорная орбита) и апогея эллипса перелёта, км.
    pub perigee_radius_km: f64,
    pub apogee_radius_km: f64,
    /// Расстояние до центра Луны в момент прилёта, км.
    pub moon_distance_km: f64,
    /// Апогей минус расстояние до Луны, км (из-за расхождения времени).
    pub apogee_miss_km: f64,
    /// Скорость в перигее после разгона и в апогее (относительно Земли), км/с.
    pub perigee_speed_km_s: f64,
    pub apogee_speed_km_s: f64,
}

impl LegModel for HohmannConfig {
    const NAME: &'static str = "hohmann";
    const SECTION: &'static str = "hohmann";
    const DESCRIPTION: &'static str = "Hohmann transfer (tangential burn at perigee, apogee at the Moon, 180 deg), patched conic (Earth gravity only on the way, Moon ephemeris DE440)";

    fn validate(&self, common: &TransferCommonConfig) -> Result<()> {
        HohmannConfig::validate(self, common)
    }

    fn time_of_flight_bounds_s(&self, common: &TransferCommonConfig) -> (f64, f64) {
        let r_p = common.parking_radius_km();
        let tol = self.tolerance_s(common.search_step_min * 60.0);
        (
            hohmann_time_of_flight_s(r_p, MOON_DISTANCE_MIN_KM, MU_EARTH_DEFAULT) - tol,
            hohmann_time_of_flight_s(r_p, MOON_DISTANCE_MAX_KM, MU_EARTH_DEFAULT) + tol,
        )
    }

    fn legs(&self, g: &CrossingGeometry, counters: &mut LegCounters, out: &mut dyn FnMut(Leg)) {
        let mu = g.mu_earth_km3_s2;
        let r_p = g.parking_radius_km;
        let moon_distance = g.moon_position_km.norm();
        // Точка разгона — напротив Луны; она лежит в плоскости опорной орбиты.
        let u = -g.moon_position_km / moon_distance;
        let theta = u.dot(&g.along).atan2(u.dot(&g.l_hat));
        let nominal = hohmann_time_of_flight_s(r_p, moon_distance, mu);
        let tol = self.tolerance_s(g.arrival_step_s);
        let period = TAU / g.n_park_rad_s;
        let launch_to_arrival_s = g.launch_to_arrival_s();
        let earliest = g.coast_min_s.max(g.ascent_s + 60.0);
        // Первый проход точки напротив Луны после старта, затем каждый виток.
        let mut coast_s = (theta - g.ascent_phase_rad).rem_euclid(TAU) / g.n_park_rad_s;
        while coast_s <= g.coast_max_s {
            if coast_s >= earliest {
                counters.solves += 1;
                let tof = launch_to_arrival_s - coast_s;
                let error = tof - nominal;
                if tof <= 0.0 || error.abs() > tol {
                    counters.rejected_timing += 1;
                } else {
                    // Эллипс, который за фактическое время полёта приходит в апогей.
                    let a = (mu * (tof / PI).powi(2)).cbrt();
                    let r_a = 2.0 * a - r_p;
                    if r_a <= r_p {
                        counters.failures += 1;
                    } else {
                        let v_p = (mu * (2.0 / r_p - 1.0 / a)).sqrt();
                        let v_a = (mu * (2.0 / r_a - 1.0 / a)).sqrt();
                        let t_hat = g.plane_normal.cross(&u);
                        out(Leg {
                            coast_s,
                            coast_angle_rad: g.ascent_phase_rad + g.n_park_rad_s * coast_s,
                            r1: u * r_p,
                            v_park: t_hat * g.circular_speed_km_s,
                            v1: t_hat * v_p,
                            // В апогее (точка −u) скорость направлена против t_hat.
                            v2: -t_hat * v_a,
                            transfer_angle_rad: PI,
                            hohmann: Some(HohmannInfo {
                                nominal_time_of_flight_s: nominal,
                                timing_error_s: error,
                                timing_tolerance_s: tol,
                                semi_major_axis_km: a,
                                perigee_radius_km: r_p,
                                apogee_radius_km: r_a,
                                moon_distance_km: moon_distance,
                                apogee_miss_km: r_a - moon_distance,
                                perigee_speed_km_s: v_p,
                                apogee_speed_km_s: v_a,
                            }),
                        });
                    }
                }
            }
            coast_s += period;
        }
    }

    fn timing_reason(&self) -> &'static str {
        "no launch opportunity reaches the Hohmann departure point (opposite the Moon) at the right time: the half-ellipse flight time does not match this arrival within timing_tolerance_min"
    }

    fn arrival_correction_s(&self, leg: &Leg) -> Option<f64> {
        let h = leg.hohmann?;
        // Разгон почти не сдвигается с моментом прилёта, поэтому прилёт сдвигается на расхождение.
        (h.timing_error_s.abs() > POLISH_TOLERANCE_S).then_some(-h.timing_error_s)
    }
}

#[cfg(test)]
mod tests {
    use super::super::*;
    use super::*;
    use crate::ephemeris::{AlmanacSource, load_almanac_blocking};
    use crate::spaceport::SpaceportCatalog;
    use anise::math::Vector3;
    use anise::prelude::Almanac;
    use anise::time::{Epoch, Unit};
    use approx::assert_abs_diff_eq;
    use std::path::Path;

    /// Синтетический проход: плоскость XY, космодром на оси X, Луна под углом 30° от −X.
    fn geometry(arrival_offset_s: f64) -> (CrossingGeometry, f64) {
        let mu = MU_EARTH_DEFAULT;
        let r_p = TransferCommonConfig::default().parking_radius_km();
        let n = (mu / r_p.powi(3)).sqrt();
        let moon_dir = Vector3::new(-(30f64.to_radians().cos()), -(30f64.to_radians().sin()), 0.0);
        let moon = moon_dir * 384_400.0;
        // Точка напротив Луны — под углом 30° от X; выход на орбиту — сразу (фаза 0).
        let coast = 30f64.to_radians() / n + TAU / n; // второй виток
        let nominal = hohmann_time_of_flight_s(r_p, 384_400.0, mu);
        let launch = Epoch::from_gregorian_utc_at_midnight(2027, 1, 1);
        let g = CrossingGeometry {
            arrival: launch + (coast + nominal + arrival_offset_s) * Unit::Second,
            launch,
            moon_position_km: moon,
            moon_velocity_km_s: Vector3::new(0.5, -0.87, 0.0),
            plane_normal: Vector3::new(0.0, 0.0, 1.0),
            l_hat: Vector3::new(1.0, 0.0, 0.0),
            along: Vector3::new(0.0, 1.0, 0.0),
            ascent_phase_rad: 0.0,
            ascent_s: 600.0,
            n_park_rad_s: n,
            parking_radius_km: r_p,
            circular_speed_km_s: (mu / r_p).sqrt(),
            mu_earth_km3_s2: mu,
            coast_min_s: 20.0 * 60.0,
            coast_max_s: 180.0 * 60.0,
            coast_step_s: 30.0,
            arrival_step_s: 60.0,
        };
        (g, coast)
    }

    fn legs(g: &CrossingGeometry) -> (Vec<Leg>, LegCounters) {
        let mut out = Vec::new();
        let mut c = LegCounters::default();
        HohmannConfig::default().legs(g, &mut c, &mut |l| out.push(l));
        (out, c)
    }

    #[test]
    fn leg_is_a_tangential_half_ellipse_to_the_moon() {
        let (g, coast) = geometry(30.0);
        let (found, c) = legs(&g);
        // Первый проход точки — раньше coast_min, третий — позже coast_max: подходит только второй.
        assert_eq!(found.len(), 1, "one revolution matches the timing");
        assert_eq!((c.solves, c.rejected_timing), (1, 0));
        let leg = found[0];
        assert_abs_diff_eq!(leg.coast_s, coast, epsilon = 1e-6);
        let h = leg.hohmann.unwrap();
        assert_abs_diff_eq!(h.timing_error_s, 30.0, epsilon = 1e-6);
        // Точка разгона напротив Луны, разгон по касательной, скорость по касательной и в апогее.
        let moon_hat = g.moon_position_km.normalize();
        assert_abs_diff_eq!(leg.r1.normalize().dot(&moon_hat), -1.0, epsilon = 1e-12);
        assert_abs_diff_eq!(leg.r1.dot(&leg.v1), 0.0, epsilon = 1e-6);
        assert_abs_diff_eq!(leg.v1.normalize().dot(&leg.v_park.normalize()), 1.0, epsilon = 1e-12);
        assert_abs_diff_eq!(leg.transfer_angle_rad, PI, epsilon = 1e-15);
        // Энергия и момент импульса: апогей эллипса — в направлении Луны на радиусе r_a.
        let mu = g.mu_earth_km3_s2;
        let r1 = leg.r1.norm();
        let energy = leg.v1.norm_squared() / 2.0 - mu / r1;
        let a = -mu / (2.0 * energy);
        assert_abs_diff_eq!(a, h.semi_major_axis_km, epsilon = 1e-6);
        let hv = leg.r1.cross(&leg.v1);
        let e_vec = leg.v1.cross(&hv) / mu - leg.r1 / r1;
        let apogee_dir = -e_vec.normalize();
        assert_abs_diff_eq!(apogee_dir.dot(&moon_hat), 1.0, epsilon = 1e-9);
        assert_abs_diff_eq!(a * (1.0 + e_vec.norm()), h.apogee_radius_km, epsilon = 1e-6);
        // Расхождение 30 с даёт промах апогея ~20 км; время полёта — ровно половина периода.
        assert!(h.apogee_miss_km > 0.0 && h.apogee_miss_km < 40.0, "{}", h.apogee_miss_km);
        assert_abs_diff_eq!(PI * (a.powi(3) / mu).sqrt(), h.nominal_time_of_flight_s + 30.0, epsilon = 1e-6);
        // Скорость в апогее — по касательной, по направлению движения (h = z).
        assert_abs_diff_eq!(leg.v2.dot(&moon_hat), 0.0, epsilon = 1e-9);
        assert!(moon_hat.cross(&leg.v2).z > 0.0);
        // Δv разгона с 185 км к Луне — около 3.1 км/с.
        let dv = (leg.v1 - leg.v_park).norm();
        assert!(dv > 3.05 && dv < 3.2, "{dv}");
        // Доводка: сдвинуть прилёт на −30 с.
        assert_abs_diff_eq!(HohmannConfig::default().arrival_correction_s(&leg).unwrap(), -30.0, epsilon = 1e-9);
    }

    #[test]
    fn timing_outside_tolerance_is_rejected() {
        let (g, _) = geometry(10.0 * 60.0);
        let (found, c) = legs(&g);
        assert!(found.is_empty());
        assert!(c.rejected_timing >= 1);
        // В грубом переборе (шаг 30 мин) допуск — не меньше половины шага.
        let coarse = CrossingGeometry { arrival_step_s: 30.0 * 60.0, ..g };
        assert_eq!(legs(&coarse).0.len(), 1);
    }

    #[test]
    fn bounds_cover_the_moon_distance_range() {
        let cfg = HohmannConfig::default();
        let common = TransferCommonConfig::default();
        let (lo, hi) = cfg.time_of_flight_bounds_s(&common);
        let r_p = common.parking_radius_km();
        for d in [356_400.0, 384_400.0, 406_700.0] {
            let t = hohmann_time_of_flight_s(r_p, d, MU_EARTH_DEFAULT);
            assert!(t > lo && t < hi);
        }
        assert!(lo > 4.0 * 86_400.0 && hi < 5.7 * 86_400.0, "{lo} {hi}");
        assert!(cfg.validate(&common).is_ok());
        assert!(HohmannConfig { timing_tolerance_min: 0.0 }.validate(&common).is_err());
    }

    // ---- С настоящими ядрами из data/kernels: `cargo test -- --ignored` ----

    fn real() -> (Almanac, SpaceportCatalog) {
        let almanac = load_almanac_blocking(AlmanacSource::LocalDir("data/kernels".into())).unwrap();
        let text = std::fs::read_to_string("data/spaceports.toml").unwrap();
        let catalog = SpaceportCatalog::parse(&text, Path::new("data/spaceports.toml")).unwrap();
        (almanac, catalog)
    }

    fn check_hohmann(s: &TransferSolution) {
        assert_eq!(s.trajectory, "hohmann");
        let h = s.hohmann.expect("hohmann info");
        assert!(h.timing_error_s.abs() <= h.timing_tolerance_s);
        assert!(h.apogee_miss_km.abs() < 100.0, "apogee miss {}", h.apogee_miss_km);
        assert_abs_diff_eq!(s.departure.transfer_angle_deg, 180.0, epsilon = 1e-9);
        assert_abs_diff_eq!(s.departure.flight_path_angle_deg, 0.0, epsilon = 1e-6);
        assert_abs_diff_eq!(s.departure.transfer_perigee_km, s.parking_orbit.radius_km, epsilon = 1e-3);
        assert_abs_diff_eq!(s.departure.transfer_apogee_km.unwrap(), h.apogee_radius_km, epsilon = 1e-3);
        let tof = s.departure.time_of_flight_s.to_seconds();
        assert_abs_diff_eq!(tof, h.nominal_time_of_flight_s + h.timing_error_s, epsilon = 1e-3);
        assert!(tof > 4.3 * 86_400.0 && tof < 5.5 * 86_400.0, "tof {tof}");
        assert!(s.departure.delta_v_magnitude_km_s > 3.0 && s.departure.delta_v_magnitude_km_s < 3.2);
        assert!(s.arrival.delta_v_loi_km_s > 0.6 && s.arrival.delta_v_loi_km_s < 1.2);
        // Точка разгона напротив Луны в момент прилёта.
        let cos = s.departure.position_km.normalize().dot(&s.arrival.moon_position_km.normalize());
        assert_abs_diff_eq!(cos, -1.0, epsilon = 1e-9);
        assert_eq!(s.launch.epoch + s.departure.coast_s, s.departure.epoch);
        assert_eq!(s.departure.epoch + s.departure.time_of_flight_s, s.arrival.epoch);
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_hohmann_transfer_both_ways() {
        let (almanac, catalog) = real();
        let ksc = catalog.get("ksc").unwrap();
        let tr = HohmannTransfer::default();
        let day = Epoch::from_gregorian_utc_at_midnight(2026, 12, 1);
        let fwd = tr.solve(&almanac, ksc, DateConstraint::LaunchDay(day), None, None).unwrap();
        assert_eq!(fwd.trajectory, "hohmann");
        assert!(fwd.config.get("hohmann").is_some() && fwd.config.get("transfer").is_some());
        let s = fwd.solution.clone().expect("hohmann solution");
        check_hohmann(&s);
        // После доводки время совпадает с точностью до секунд.
        assert!(s.hohmann.unwrap().timing_error_s.abs() < 5.0, "{:?}", s.hohmann);
        assert!(s.launch.epoch >= day && s.launch.epoch < day + 1 * Unit::Day);
        assert!(ksc.azimuth_allowed(s.launch.azimuth_deg));

        // Гоман медленнее и не экономичнее по разгону, чем лучший прямой перелёт в эти сутки
        // (прямой перелёт с тем же временем полёта — частный случай поиска Ламберта).
        let direct = DirectTransfer::default().solve(&almanac, ksc, DateConstraint::LaunchDay(day), None, None).unwrap();
        let d = direct.solution.unwrap();
        assert!(d.delta_v_total_km_s <= s.delta_v_total_km_s + 0.02, "direct {} vs hohmann {}", d.delta_v_total_km_s, s.delta_v_total_km_s);

        // Обратная задача: по суткам прилёта находится перелёт Гомана с прилётом в эти сутки.
        let (y, m, dd, ..) = s.arrival.epoch.to_gregorian_utc();
        let arr_day = Epoch::from_gregorian_utc_at_midnight(y, m, dd);
        let inv = tr.solve(&almanac, ksc, DateConstraint::ArrivalDay(arr_day), None, None).unwrap();
        let si = inv.solution.expect("inverse hohmann solution");
        check_hohmann(&si);
        assert!(si.arrival.epoch >= arr_day && si.arrival.epoch < arr_day + 1 * Unit::Day);

        // solve_at по моменту прилёта из перебора даёт вариант Гомана.
        let p = fwd.coarse_scan.iter().find(|p| p.feasible).unwrap();
        let alt = tr.solve_at(&almanac, ksc, DateConstraint::LaunchDay(day), None, None, p.arrival).unwrap().unwrap();
        check_hohmann(&alt);
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_hohmann_with_landing_site_and_budget() {
        use crate::lunar_landing::{LandingTarget, LandingWindowConfig};
        use crate::visibility::{SiteVisibility, VisibilityConfig};
        let (almanac, catalog) = real();
        let ksc = catalog.get("ksc").unwrap();
        let tr = HohmannTransfer::default();
        let target = LandingTarget {
            site: SiteVisibility::new(None, 5.0, -10.0, 0.0, VisibilityConfig::default()).unwrap(),
            window: LandingWindowConfig::default(),
            descent: crate::lunar_landing::DescentModel::default(),
        };
        let day = Epoch::from_gregorian_utc_at_midnight(2027, 1, 1);
        let rep = tr.solve(&almanac, ksc, DateConstraint::ArrivalDay(day), None, Some(&target)).unwrap();
        let s = rep.solution.expect("landing on New Year's day with a Hohmann transfer");
        check_hohmann(&s);
        let l = s.landing.expect("landing info");
        assert!(l.epoch >= day && l.epoch < day + 1 * Unit::Day, "landing {}", l.epoch);

        // Слишком маленький бюджет: перелёта нет, причина — бюджет.
        let tiny = DeltaVBudget::new("tiny", 1500.0, &MissionConfig::default()).unwrap();
        let none = tr.solve(&almanac, ksc, DateConstraint::LaunchDay(day), Some(&tiny), None).unwrap();
        assert!(none.solution.is_none());
        assert_eq!(none.no_solution_reason, Some(REASON_BUDGET));
    }

    #[test]
    #[ignore = "needs kernels in data/kernels"]
    fn real_hohmann_refines_numerically() {
        let (almanac, catalog) = real();
        let ksc = catalog.get("ksc").unwrap();
        let day = Epoch::from_gregorian_utc_at_midnight(2026, 12, 1);
        let mut s = HohmannTransfer::default()
            .solve(&almanac, ksc, DateConstraint::LaunchDay(day), None, None)
            .unwrap()
            .solution
            .unwrap();
        let cfg = crate::refine::NumericalConfig::default();
        let refined = crate::refine::refine(&almanac, &cfg, &s, None, None, None).unwrap().expect("refined");
        // Уточнённый разгон близок к касательному по формулам.
        assert!((refined.tli_delta_v_km_s - s.departure.delta_v_magnitude_km_s).abs() < 0.1);
        s.refined = Some(refined);
    }
}
