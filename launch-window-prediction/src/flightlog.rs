//! Разделы-сводки журнала информации о полёте ([`crate::runlog`]): место старта,
//! место посадки и короткий итог. Всё здесь собирается из уже посчитанных результатов;
//! полные данные остаются в журнале событий.

use anise::constants::frames::{EARTH_ITRF93, EARTH_J2000, IAU_EARTH_FRAME, SUN_J2000};
use anise::math::Vector3;
use anise::prelude::Almanac;
use anise::time::{Epoch, Unit};
use serde_json::{Value, json};

use crate::landing::{LandingAssessment, LandingConfig};
use crate::runlog::{human_duration, iso_short};
use crate::spaceport::Spaceport;
use crate::trajectory::{TransferReport, TransferSolution, geodetic_to_ecef};
use crate::visibility::{BodyVisibility, SiteVisibility, VisibilitySummary};

/// Угловая скорость вращения Земли, рад/с.
const EARTH_ROTATION_RAD_S: f64 = 7.292_115_146_7e-5;

/// Направление по азимуту: N, NNE, NE, ... (16 румбов).
pub fn compass(azimuth_deg: f64) -> &'static str {
    const NAMES: [&str; 16] =
        ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"];
    let i = ((azimuth_deg.rem_euclid(360.0) / 22.5).round() as usize) % 16;
    NAMES[i]
}

/// `5.0000° N, 10.0000° W`.
pub fn lat_lon_text(lat_deg: f64, lon_deg: f64) -> String {
    format!(
        "{:.4}° {}, {:.4}° {}",
        lat_deg.abs(),
        if lat_deg >= 0.0 { "N" } else { "S" },
        lon_deg.abs(),
        if lon_deg >= 0.0 { "E" } else { "W" }
    )
}

/// Местное среднее солнечное время по долготе: `HH:MM`.
fn local_mean_solar_time(t: Epoch, lon_deg: f64) -> String {
    let (_, _, _, h, m, s, ns) = t.to_gregorian_utc();
    let utc_h = h as f64 + m as f64 / 60.0 + (s as f64 + ns as f64 * 1e-9) / 3600.0;
    let local = (utc_h + lon_deg / 15.0).rem_euclid(24.0);
    let minutes = (local * 60.0).round() as u32 % (24 * 60);
    format!("{:02}:{:02}", minutes / 60, minutes % 60)
}

/// Время от старта для человека: `T+1 h 17 min 30 s`.
fn t_plus(launch: Epoch, t: Epoch) -> String {
    format!("T+{}", human_duration((t - launch).to_seconds()))
}

/// Солнце над космодромом в момент `t`: азимут и высота над горизонтом (геодезическим), градусы.
pub fn sun_at_earth_site(almanac: &Almanac, p: &Spaceport, t: Epoch) -> Option<(f64, f64)> {
    let sun = almanac.translate(SUN_J2000, EARTH_J2000, t, None).ok()?.radius_km;
    let dcm = almanac
        .rotate(EARTH_ITRF93, EARTH_J2000, t)
        .or_else(|_| almanac.rotate(IAU_EARTH_FRAME, EARTH_J2000, t))
        .ok()?
        .rot_mat;
    let sun_fixed = dcm.transpose() * sun;
    let site = geodetic_to_ecef(p.lat_deg, p.lon_deg, p.altitude_m.unwrap_or(0.0) / 1000.0);
    let d = sun_fixed - site;
    let (la, lo) = (p.lat_deg.to_radians(), p.lon_deg.to_radians());
    let east = Vector3::new(-lo.sin(), lo.cos(), 0.0);
    let north = Vector3::new(-la.sin() * lo.cos(), -la.sin() * lo.sin(), la.cos());
    let up = Vector3::new(la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin());
    let el = (d.dot(&up) / d.norm()).asin().to_degrees();
    let az = d.dot(&east).atan2(d.dot(&north)).to_degrees().rem_euclid(360.0);
    Some((az, el))
}

fn daylight(sun_elevation_deg: f64) -> &'static str {
    match sun_elevation_deg {
        e if e > 0.0 => "day (Sun above the horizon)",
        e if e > -6.0 => "civil twilight",
        e if e > -18.0 => "twilight",
        _ => "night",
    }
}

/// Сведения о космодроме, не зависящие от варианта перелёта.
pub fn spaceport_info(p: &Spaceport) -> Value {
    let alt_km = p.altitude_m.unwrap_or(0.0) / 1000.0;
    let ecef = geodetic_to_ecef(p.lat_deg, p.lon_deg, alt_km);
    let rho = ecef.x.hypot(ecef.y);
    json!({
        "id": p.id,
        "name": p.name,
        "country": p.country,
        "pads": p.pads,
        "location": lat_lon_text(p.lat_deg, p.lon_deg),
        "latitude_deg": p.lat_deg,
        "longitude_deg": p.lon_deg,
        "altitude_m": p.altitude_m,
        "altitude_note": if p.altitude_m.is_none() { Some("altitude not given in spaceports.toml; 0 m (sea level) is used") } else { None },
        "position": {
            "ellipsoid": "WGS84",
            "earth_fixed_km": [ecef.x, ecef.y, ecef.z],
            "earth_fixed_frame": "ITRF93",
            "distance_from_earth_centre_km": ecef.norm(),
            "geocentric_latitude_deg": (ecef.z / ecef.norm()).asin().to_degrees(),
        },
        "earth_rotation_speed_km_s": EARTH_ROTATION_RAD_S * rho,
        "lowest_orbit_inclination_without_plane_change_deg": p.lat_deg.abs(),
        "range_limits": {
            "azimuth_min_deg": p.azimuth_min_deg,
            "azimuth_max_deg": p.azimuth_max_deg,
            "inclination_min_deg": p.inclination_min_deg,
            "inclination_max_deg": p.inclination_max_deg,
            "allowed_inclinations_deg": p.allowed_inclinations_deg,
            "source": p.limits_source,
        },
    })
}

/// Раздел `launch_site`: космодром и всё о старте с него в найденном варианте.
pub fn launch_site(r: &TransferReport, almanac: Option<&Almanac>) -> Value {
    let p = &r.spaceport;
    let mut v = spaceport_info(p);
    let obj = v.as_object_mut().expect("object");
    obj.insert(
        "launch_window_searched".into(),
        json!({ "from_utc": iso_short(&r.launch_window_start), "to_utc": iso_short(&r.launch_window_end) }),
    );
    let Some(s) = &r.solution else {
        obj.insert(
            "launch".into(),
            json!({ "possible": false, "reason": r.no_solution_reason.unwrap_or("no transfer found") }),
        );
        return v;
    };
    let l = &s.launch;
    let a = &s.ascent;
    let po = &s.parking_orbit;
    let sun = almanac.and_then(|al| sun_at_earth_site(al, p, l.epoch));
    let (start, _) = crate::trajlog::endpoints(r);
    let sim = a.simulation.as_ref();
    let launch = json!({
        "possible": true,
        "time_utc": iso_short(&l.epoch),
        "local_mean_solar_time": local_mean_solar_time(l.epoch, p.lon_deg),
        "sun_at_site": sun.map(|(az, el)| json!({
            "azimuth_deg": az,
            "elevation_deg": el,
            "lighting": daylight(el),
        })),
        "azimuth_deg": l.azimuth_deg,
        "direction": compass(l.azimuth_deg),
        "pass": if l.ascending_pass { "ascending (the site crosses the orbit plane moving north)" } else { "descending (the site crosses the orbit plane moving south)" },
        "earth_rotation_gain_km_s": l.earth_rotation_gain_km_s,
        "rocket": a.rocket,
        "ascent": {
            "model": a.model,
            "model_note": a.model_note,
            "duration_s": a.duration_s,
            "duration": human_duration(a.duration_s),
            "downrange_km": a.downrange_km,
            "orbit_insertion_utc": iso_short(&a.insertion_epoch),
            "liftoff_mass_kg": sim.map(|x| x.liftoff_mass_kg),
            "liftoff_thrust_to_weight": sim.map(|x| x.liftoff_thrust_to_weight),
            "max_q_pa": sim.map(|x| x.max_q_pa),
            "max_acceleration_g": sim.map(|x| x.max_acceleration_g),
            "mass_in_orbit_kg": sim.map(|x| x.mass_in_orbit_kg),
        },
        "parking_orbit": {
            "altitude_km": po.altitude_km,
            "inclination_deg": po.inclination_deg,
            "raan_deg": po.raan_deg,
            "period_min": po.period_s / 60.0,
            "speed_km_s": po.circular_speed_km_s,
        },
        "liftoff_point": start,
    });
    obj.insert("launch".into(), launch);
    v
}

/// Начало раздела `landing_site`: запрошенная точка.
pub fn landing_request(lat_deg: f64, lon_deg: f64) -> Value {
    json!({
        "location": lat_lon_text(lat_deg, lon_deg),
        "latitude_deg": lat_deg,
        "longitude_deg": lon_deg,
        "side_of_moon": if lon_deg.abs() <= 90.0 { "near side (faces the Earth)" } else { "far side (faces away from the Earth)" },
        "coordinates": "Moon mean-Earth/polar-axis frame (MOON_ME), longitude positive east",
    })
}

/// Проверка площадки по рельефу для раздела `landing_site`.
pub fn landing_terrain(a: &LandingAssessment, cfg: &LandingConfig) -> Value {
    let centre = a.samples.iter().find(|p| p.east_m == 0.0 && p.north_m == 0.0).or(a.samples.first());
    json!({
        "terrain_check": {
            "verdict": a.verdict,
            "explanation": a.explain(cfg),
            "map": a.map_name,
            "map_resolution_m": a.resolution_m,
            "surface_height_m": centre.map(|p| p.height_m),
            "surface_height_reference": "above the 1737.4 km lunar sphere",
            "slope_baseline_m": a.baseline_m,
            "max_slope_deg": a.max_slope_deg,
            "mean_slope_deg": a.mean_slope_deg,
            "steep_fraction": a.steep_fraction,
            "steepest_at": a.steepest_at.map(|(la, lo)| json!({ "latitude_deg": la, "longitude_deg": lo })),
            "relief_m": a.relief_m,
            "points_checked": a.points,
            "rules": cfg,
        },
    })
}

/// Положение, высота и горизонт точки для раздела `landing_site`.
pub fn landing_surface(site: &SiteVisibility) -> Value {
    let pos = site.site.position_km();
    let horizon = site.horizon_mask_deg().filter(|m| !m.is_empty()).map(|m| {
        let step = 360.0 / m.len() as f64;
        let (imax, max) = m.iter().copied().enumerate().fold((0, f64::MIN), |b, (i, x)| if x > b.1 { (i, x) } else { b });
        let (imin, min) = m.iter().copied().enumerate().fold((0, f64::MAX), |b, (i, x)| if x < b.1 { (i, x) } else { b });
        json!({
            "highest_deg": max,
            "highest_at_azimuth_deg": imax as f64 * step,
            "highest_direction": compass(imax as f64 * step),
            "lowest_deg": min,
            "lowest_at_azimuth_deg": imin as f64 * step,
            "mean_deg": m.iter().sum::<f64>() / m.len() as f64,
            "azimuth_step_deg": step,
            "full_mask": "event visibility_site in the log",
        })
    });
    json!({
        "observer": {
            "height_km": site.site.height_km,
            "height_note": "above the 1737.4 km sphere, including the observer height above the surface",
            "observer_height_above_surface_m": site.config.horizon.observer_height_m,
            "moon_fixed_km": [pos.x, pos.y, pos.z],
            "distance_from_moon_centre_km": pos.norm(),
        },
        "surface": site.terrain,
        "horizon": horizon,
        "visibility_rules": site.config,
    })
}

fn body(b: &BodyVisibility) -> Value {
    json!({
        "azimuth_deg": b.azimuth_deg,
        "direction": compass(b.azimuth_deg),
        "elevation_deg": b.elevation_deg,
        "terrain_horizon_deg": b.horizon_deg,
        "visible_fraction": b.visible_fraction,
        "visible": b.visible,
    })
}

/// Посадка в найденном варианте для раздела `landing_site`.
pub fn landing_result(r: &TransferReport) -> Value {
    let Some(s) = &r.solution else {
        return json!({
            "landing": { "possible": false, "reason": r.no_solution_reason.unwrap_or("no transfer found") },
            "landing_window_rules": r.landing_window,
        });
    };
    let Some(f) = &s.landing else {
        return json!({ "landing": { "possible": false, "reason": "landing was not calculated" } });
    };
    let launch = s.launch.epoch;
    let refined = s.refined.as_ref();
    let rl = refined.and_then(|n| n.landing.as_ref());
    let (_, end) = crate::trajlog::endpoints(r);
    let arrival = refined.map_or(s.arrival.epoch, |n| n.arrival_epoch);
    let doi = rl.map_or(f.doi_epoch, |l| l.doi_epoch);
    let pdi = rl.map_or(f.pdi_epoch, |l| l.pdi_epoch);
    let touchdown = rl.map_or(f.epoch, |l| l.touchdown_epoch);
    let wait_s = rl.map_or(f.wait_after_arrival_s, |l| l.wait_in_orbit_s);
    let (sun, earth) = rl.map_or((&f.sun, &f.earth), |l| (&l.sun, &l.earth));
    let sun_v = {
        let mut b = body(sun);
        // По азимуту Солнца: на востоке — лунное утро, на западе — вечер.
        let part = if sun.azimuth_deg.rem_euclid(360.0) < 180.0 { "lunar morning (Sun in the east)" } else { "lunar evening (Sun in the west)" };
        b.as_object_mut().expect("object").insert("time_of_lunar_day".into(), json!(part));
        b
    };
    json!({
        "landing": {
            "possible": true,
            "values_from": if rl.is_some() { "numerical refinement" } else { "formulas (patched conics)" },
            "touchdown_utc": iso_short(&touchdown),
            "touchdown_t_plus": t_plus(launch, touchdown),
            "lunar_orbit_insertion_utc": iso_short(&arrival),
            "deorbit_burn_utc": iso_short(&doi),
            "powered_descent_start_utc": iso_short(&pdi),
            "wait_in_lunar_orbit_s": wait_s,
            "wait_in_lunar_orbit": human_duration(wait_s),
            "orbit_pass": f.orbit_pass + rl.map_or(0, |l| l.extra_orbits),
            "lunar_orbit": {
                "altitude_km": s.arrival.lunar_orbit_altitude_km,
                "inclination_to_lunar_equator_deg": rl.map_or(f.lunar_orbit_inclination_deg, |l| l.orbit_inclination_deg),
                "period_min": f.lunar_orbit_period_s / 60.0,
                "min_altitude_km": rl.map(|l| l.orbit_min_altitude_km),
                "max_altitude_km": rl.map(|l| l.orbit_max_altitude_km),
            },
            "descent": {
                "deorbit_delta_v_m_s": rl.map_or(f.descent.doi_delta_v_m_s, |l| l.doi_delta_v_m_s),
                "powered_descent_start_altitude_km": rl.map(|l| l.pdi_altitude_km),
                "powered_descent_start_height_above_site_km": rl.map(|l| l.pdi_height_above_site_km),
                "powered_descent_start_speed_km_s": rl.map_or(f.descent.pdi_speed_km_s, |l| l.pdi_speed_km_s),
                "braking_duration_s": rl.map_or(f.descent.braking_duration_s, |l| l.braking.duration_s),
                "braking_range_km": rl.map_or(f.descent.powered_descent_range_km, |l| l.braking.range_km),
                "terminal_phase_s": rl.map_or(f.descent.powered_descent_duration_s - f.descent.braking_duration_s, |l| l.terminal_duration_s),
                "descent_delta_v_m_s": rl.map_or(f.descent.total_delta_v_m_s, |l| l.descent_delta_v_m_s),
                "mass_after_braking_kg": rl.map(|l| l.braking.final_mass_kg),
            },
            "miss_distance_m": rl.map(|l| l.miss_km * 1000.0),
            "touchdown_point": end,
        },
        "conditions_at_touchdown": {
            "sun": sun_v,
            "earth": body(earth),
            "all_conditions_met": rl.map_or(true, |l| l.conditions_met),
            "unmet_reason": rl.and_then(|l| l.unmet_reason),
        },
        "landing_window_rules": r.landing_window,
    })
}

/// Видимость Солнца и Земли после посадки для раздела `landing_site`.
pub fn landing_visibility(start: Epoch, end: Epoch, step_s: f64, s: &VisibilitySummary) -> Value {
    json!({
        "visibility_after_landing": {
            "from_utc": iso_short(&start),
            "to_utc": iso_short(&end),
            "step_min": step_s / 60.0,
            "samples": s.samples,
            "sun_visible_percent": s.sun_pct,
            "earth_visible_percent": s.earth_pct,
            "both_visible_percent": s.both_pct,
            "longest_without_sun": human_duration(s.longest_sun_gap.to_seconds()),
            "longest_without_earth": human_duration(s.longest_earth_gap.to_seconds()),
            "longest_without_both": human_duration(s.longest_both_gap.to_seconds()),
            "every_sample": "events visibility_sample in the log",
        },
    })
}

fn timeline(s: &TransferSolution) -> Vec<Value> {
    let launch = s.launch.epoch;
    let refined = s.refined.as_ref();
    let rl = refined.and_then(|n| n.landing.as_ref());
    let mut rows = vec![
        ("launch", "liftoff from the pad", launch),
        ("parking_orbit_insertion", "end of ascent, in the parking orbit", s.ascent.insertion_epoch),
        ("tli", "trans-lunar injection burn", s.departure.epoch),
    ];
    // Разгон конечной длительности: начало и конец работы двигателя вместо момента импульса.
    if let Some(b) = refined.and_then(|n| n.tli_burn.as_ref()) {
        rows.pop();
        rows.push(("tli_start", "trans-lunar injection burn starts", b.start_epoch));
        rows.push(("tli_end", "trans-lunar injection burn ends", b.end_epoch));
    }
    rows.push(("loi", "lunar orbit insertion burn", refined.map_or(s.arrival.epoch, |n| n.arrival_epoch)));
    if let Some(f) = &s.landing {
        rows.push(("doi", "deorbit burn", rl.map_or(f.doi_epoch, |l| l.doi_epoch)));
        rows.push(("pdi", "powered descent starts", rl.map_or(f.pdi_epoch, |l| l.pdi_epoch)));
        if let Some(l) = rl {
            rows.push((
                "braking_end",
                "braking ends above the site",
                l.pdi_epoch + l.braking.duration_s * Unit::Second,
            ));
        }
        rows.push(("touchdown", "touchdown at the site", rl.map_or(f.epoch, |l| l.touchdown_epoch)));
    }
    rows.into_iter()
        .map(|(event, what, t)| json!({ "event": event, "what": what, "time_utc": iso_short(&t), "t_plus": t_plus(launch, t) }))
        .collect()
}

/// Раздел `summary`: короткий итог. `reports` отсортированы (лучший первым).
pub fn summary(reports: &[TransferReport], spacecraft: Option<&str>, rocket: Option<&str>) -> Value {
    let best = reports.iter().find(|r| r.solution.is_some()).or(reports.first());
    let Some(best) = best else { return json!({ "status": "no spaceports were searched" }) };
    let alternatives: Vec<Value> = reports
        .iter()
        .map(|r| match &r.solution {
            Some(s) => json!({
                "spaceport": r.spaceport.id,
                "feasible": true,
                "launch_utc": iso_short(&s.launch.epoch),
                "landing_utc": s
                    .refined
                    .as_ref()
                    .and_then(|n| n.landing.as_ref().map(|l| l.touchdown_epoch))
                    .or(s.landing.map(|l| l.epoch))
                    .map(|t| iso_short(&t)),
                "delta_v_total_km_s": s.delta_v_total_km_s,
                "refined": s.refined.is_some(),
            }),
            None => json!({ "spaceport": r.spaceport.id, "feasible": false, "reason": r.no_solution_reason }),
        })
        .collect();
    let Some(s) = &best.solution else {
        return json!({
            "status": "no flight found",
            "reason": best.no_solution_reason,
            "spaceports_searched": alternatives,
        });
    };
    let refined = s.refined.as_ref();
    let rl = refined.and_then(|n| n.landing.as_ref());
    let launch = s.launch.epoch;
    let arrival = refined.map_or(s.arrival.epoch, |n| n.arrival_epoch);
    let touchdown = rl.map(|l| l.touchdown_epoch).or(s.landing.map(|l| l.epoch));
    let tof_s = refined.map_or(s.departure.time_of_flight_s.to_seconds(), |n| n.time_of_flight_s);
    let tli = refined.map_or(s.departure.delta_v_magnitude_km_s, |n| n.tli_delta_v_km_s);
    let loi = refined.map_or(s.arrival.delta_v_loi_km_s, |n| n.loi_delta_v_km_s);
    let descent_m_s = rl.map(|l| l.descent_delta_v_m_s).or(s.landing.map(|l| l.descent.total_delta_v_m_s));
    let budget = refined.and_then(|n| n.budget.as_ref()).or(s.budget.as_ref());
    let site = best.landing_site.map(|l| lat_lon_text(l.lat_deg, l.lon_deg));
    json!({
        "status": if touchdown.is_some() { "flight and landing found" } else { "flight found" },
        "values_from": if refined.is_some() { "numerical refinement" } else { "formulas (patched conics)" },
        "not_refined_because": s.refine_failure.as_deref().or(s.refine_skipped.as_deref()),
        "spaceport": format!("{} ({})", best.spaceport.name, best.spaceport.id),
        "rocket": rocket.map(str::to_string).or_else(|| s.ascent.rocket.clone()),
        "spacecraft": spacecraft,
        "landing_site": site,
        "date_constraint": best.constraint,
        "launch_utc": iso_short(&launch),
        "lunar_orbit_insertion_utc": iso_short(&arrival),
        "touchdown_utc": touchdown.map(|t| iso_short(&t)),
        "timeline": timeline(s),
        "durations": {
            "earth_parking_orbit": human_duration(s.departure.coast_s.to_seconds()),
            "earth_to_moon": human_duration(tof_s),
            "launch_to_lunar_orbit": human_duration((arrival - launch).to_seconds()),
            "launch_to_touchdown": touchdown.map(|t| human_duration((t - launch).to_seconds())),
            "earth_to_moon_h": tof_s / 3600.0,
        },
        "delta_v": {
            "tli_m_s": tli * 1000.0,
            "tli_burn": refined.and_then(|n| n.tli_burn.as_ref()).map(|b| json!({
                "by": b.engine.name,
                "start_utc": iso_short(&b.start_epoch),
                "end_utc": iso_short(&b.end_epoch),
                "duration": human_duration(b.duration_s),
                "engine_delta_v_m_s": b.delta_v_km_s * 1000.0,
                "impulsive_equivalent_m_s": b.impulsive_delta_v_km_s * 1000.0,
                "gravity_loss_m_s": b.gravity_loss_m_s,
                "propellant_used_kg": b.propellant_used_kg,
                "propellant_available_kg": b.engine.propellant_kg,
            })),
            "lunar_orbit_insertion_m_s": loi * 1000.0,
            "descent_m_s": descent_m_s,
            "tli_plus_loi_m_s": (tli + loi) * 1000.0,
            "spacecraft_budget": budget.map(|b| json!({
                "maneuvers": b.maneuvers,
                "required_with_reserve_m_s": b.required_with_reserve_m_s,
                "available_m_s": b.available_m_s,
                "margin_m_s": b.margin_m_s,
                "enough": b.margin_m_s >= 0.0,
            })),
        },
        "moon_distance_at_arrival_km": s.arrival.moon_distance_km,
        "variants_checked_by_refinement": (!best.refine_attempts.is_empty()).then_some(&best.refine_attempts),
        "spaceports_compared": alternatives,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compass_and_text() {
        assert_eq!(compass(0.0), "N");
        assert_eq!(compass(96.3), "E");
        assert_eq!(compass(268.7), "W");
        assert_eq!(compass(359.0), "N");
        assert_eq!(lat_lon_text(5.0, -10.0), "5.0000° N, 10.0000° W");
    }

    #[test]
    fn local_time_from_longitude() {
        let t = Epoch::from_gregorian_utc_hms(2026, 12, 29, 5, 24, 12);
        assert_eq!(local_mean_solar_time(t, -80.604333), "00:02");
        assert_eq!(local_mean_solar_time(t, 110.95), "12:48");
    }
}
