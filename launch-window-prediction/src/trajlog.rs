//! Журнал траектории для построения полёта в 3D.
//!
//! Отдельно от журнала информации о полёте (`flight-<время>.json`) каждый запуск пишет
//! `trajectory-<время>.json` — один JSON-документ только с тем, что нужно для 3D:
//!
//! - система координат, единицы, радиусы и гравитационные параметры тел;
//! - космодром и точка посадки;
//! - для каждого найденного перелёта: ключевые события (старт, выход на опорную орбиту,
//!   разгон, выход на окололунную орбиту, сход с орбиты, начало торможения, остановка, касание)
//!   и точки траектории по этапам с положением Луны и ориентацией Земли и Луны в каждой точке.
//!
//! Точки берутся из численного интегрирования ([`crate::refine`]) — того же расчёта, что
//! уточняет перелёт. Если уточнение выключено или не удалось, точек нет, и причина указана.
//!
//! Все положения и скорости — в системе J2000 (ICRF) с началом в центре Земли, км и км/с.
//! Для этапов у Луны дополнительно дано положение относительно центра Луны (`r_moon`).
//! Ориентация тел — кватернион `[w, x, y, z]` поворота из связанной с телом системы в J2000:
//! вектор в J2000 = q · вектор в системе тела · q⁻¹.

use anise::constants::frames::{EARTH_ITRF93, EARTH_J2000, IAU_EARTH_FRAME, MOON_J2000, MOON_ME_DE440_ME421_FRAME, SUN_J2000};
use anise::math::{Matrix3, Vector3};
use anise::prelude::Almanac;
use anise::time::{Epoch, Unit};
use serde::Serialize;

use crate::error::{Error, Result};
use crate::propagate::EARTH_J2_RADIUS_KM;
use crate::runlog::iso;
use crate::sky::MOON_RADIUS_KM;
use crate::trajectory::TransferReport;

/// Версия формата файла.
pub const FORMAT: &str = "lwp-trajectory/1";

/// Этап полёта, к которому относится точка.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Phase {
    /// Подъём от стартового стола до выхода на опорную орбиту (кинематическая модель).
    Ascent,
    /// Опорная орбита Земли от выведения до разгона.
    ParkingOrbit,
    /// Разгон к Луне с работающим двигателем (если разгон моделируется с конечной длительностью).
    TliBurn,
    /// Перелёт от разгона до периселения (выхода на окололунную орбиту).
    Transfer,
    /// Окололунная орбита от выхода на неё до схода.
    LunarOrbit,
    /// Пассивный полёт по эллипсу спуска от схода с орбиты до начала торможения.
    DescentCoast,
    /// Торможение двигателем до остановки над точкой посадки.
    Braking,
    /// Заход и вертикальный участок: в расчёте это запас времени у точки, не моделирование.
    Terminal,
}

/// Точка траектории.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Sample {
    /// Момент, UTC.
    pub t: String,
    /// Секунды от старта.
    pub t_s: f64,
    pub phase: Phase,
    /// Положение относительно центра Земли, J2000, км.
    pub r: [f64; 3],
    /// Скорость относительно центра Земли, J2000, км/с.
    pub v: [f64; 3],
    /// Положение центра Луны относительно центра Земли, J2000, км.
    pub moon_r: [f64; 3],
    /// Скорость центра Луны относительно центра Земли, J2000, км/с.
    pub moon_v: [f64; 3],
    /// Положение относительно центра Луны, J2000, км (только у этапов возле Луны).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub r_moon: Option<[f64; 3]>,
    /// Скорость относительно центра Луны, J2000, км/с (только у этапов возле Луны).
    /// Именно её показывать как «скорость» у Луны: `v` относительно Земли качается
    /// на ±1 км/с с каждым витком из-за движения самой Луны.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub v_moon: Option<[f64; 3]>,
    /// Ориентация Земли (только на подъёме и опорной орбите).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub earth_q: Option<[f64; 4]>,
    /// Ориентация Луны, система MOON_ME (только у этапов возле Луны).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub moon_q: Option<[f64; 4]>,
}

/// Ключевое событие полёта.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Event {
    /// launch, parking_insertion, tli, loi, doi, pdi, braking_end, touchdown.
    pub name: &'static str,
    pub description: String,
    pub t: String,
    pub t_s: f64,
    pub r: [f64; 3],
    /// Скорость после события (после импульса, если он есть), км/с.
    pub v: [f64; 3],
    /// Импульс скорости события, J2000, км/с (если это манёвр).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub delta_v: Option<[f64; 3]>,
    pub moon_r: [f64; 3],
    /// Скорость центра Луны относительно центра Земли, J2000, км/с.
    pub moon_v: [f64; 3],
    #[serde(skip_serializing_if = "Option::is_none")]
    pub r_moon: Option<[f64; 3]>,
    /// Скорость после события относительно центра Луны, J2000, км/с (у событий возле Луны).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub v_moon: Option<[f64; 3]>,
    pub earth_q: [f64; 4],
    pub moon_q: [f64; 4],
    /// Направление на Солнце из центра Земли, единичный вектор J2000.
    pub sun_dir: [f64; 3],
}

/// Точки и события одного перелёта (до объединения в журнал).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TrajectoryPath {
    pub events: Vec<Event>,
    pub samples: Vec<Sample>,
}

fn r3(v: &Vector3, scale: f64) -> [f64; 3] {
    [(v.x * scale).round() / scale, (v.y * scale).round() / scale, (v.z * scale).round() / scale]
}

/// Положение с точностью до метра.
pub fn km(v: &Vector3) -> [f64; 3] {
    r3(v, 1e3)
}

/// Скорость с точностью до мм/с.
pub fn km_s(v: &Vector3) -> [f64; 3] {
    r3(v, 1e6)
}

/// Кватернион `[w, x, y, z]` из матрицы поворота.
pub fn quaternion(m: &Matrix3) -> [f64; 4] {
    let tr = m[(0, 0)] + m[(1, 1)] + m[(2, 2)];
    let (w, x, y, z) = if tr > 0.0 {
        let s = (tr + 1.0).sqrt() * 2.0;
        (0.25 * s, (m[(2, 1)] - m[(1, 2)]) / s, (m[(0, 2)] - m[(2, 0)]) / s, (m[(1, 0)] - m[(0, 1)]) / s)
    } else if m[(0, 0)] > m[(1, 1)] && m[(0, 0)] > m[(2, 2)] {
        let s = (1.0 + m[(0, 0)] - m[(1, 1)] - m[(2, 2)]).sqrt() * 2.0;
        ((m[(2, 1)] - m[(1, 2)]) / s, 0.25 * s, (m[(0, 1)] + m[(1, 0)]) / s, (m[(0, 2)] + m[(2, 0)]) / s)
    } else if m[(1, 1)] > m[(2, 2)] {
        let s = (1.0 + m[(1, 1)] - m[(0, 0)] - m[(2, 2)]).sqrt() * 2.0;
        ((m[(0, 2)] - m[(2, 0)]) / s, (m[(0, 1)] + m[(1, 0)]) / s, 0.25 * s, (m[(1, 2)] + m[(2, 1)]) / s)
    } else {
        let s = (1.0 + m[(2, 2)] - m[(0, 0)] - m[(1, 1)]).sqrt() * 2.0;
        ((m[(1, 0)] - m[(0, 1)]) / s, (m[(0, 2)] + m[(2, 0)]) / s, (m[(1, 2)] + m[(2, 1)]) / s, 0.25 * s)
    };
    let sign = if w < 0.0 { -1.0 } else { 1.0 };
    let q = [w * sign, x * sign, y * sign, z * sign];
    q.map(|c| (c * 1e12).round() / 1e12)
}

/// Положение и ориентация тел на момент — для точек и событий.
pub struct Bodies<'a> {
    pub almanac: &'a Almanac,
    pub t0: Epoch,
}

impl Bodies<'_> {
    pub fn epoch(&self, t_s: f64) -> Epoch {
        self.t0 + t_s * Unit::Second
    }

    pub fn moon(&self, t_s: f64) -> Result<(Vector3, Vector3)> {
        let e = self.epoch(t_s);
        let s = self
            .almanac
            .translate(MOON_J2000, EARTH_J2000, e, None)
            .map_err(|err| Error::Almanac(format!("Moon from Earth at {e}: {err}")))?;
        Ok((s.radius_km, s.velocity_km_s))
    }

    pub fn sun_dir(&self, t_s: f64) -> Result<Vector3> {
        let e = self.epoch(t_s);
        let s = self
            .almanac
            .translate(SUN_J2000, EARTH_J2000, e, None)
            .map_err(|err| Error::Almanac(format!("Sun from Earth at {e}: {err}")))?;
        Ok(s.radius_km.normalize())
    }

    /// Матрица перехода из связанной с Землёй системы (ITRF93, а без точных данных IAU_EARTH) в J2000.
    pub fn earth_dcm(&self, t_s: f64) -> Result<Matrix3> {
        let e = self.epoch(t_s);
        let dcm = self
            .almanac
            .rotate(EARTH_ITRF93, EARTH_J2000, e)
            .or_else(|_| self.almanac.rotate(IAU_EARTH_FRAME, EARTH_J2000, e))
            .map_err(|err| Error::Almanac(format!("Earth orientation at {e}: {err}")))?;
        Ok(dcm.rot_mat)
    }

    /// Матрица перехода из MOON_ME в J2000.
    pub fn moon_dcm(&self, t_s: f64) -> Result<Matrix3> {
        let e = self.epoch(t_s);
        let dcm = self
            .almanac
            .rotate(MOON_ME_DE440_ME421_FRAME, MOON_J2000, e)
            .map_err(|err| Error::Almanac(format!("Moon orientation at {e}: {err}")))?;
        Ok(dcm.rot_mat)
    }

    /// Точка по состоянию относительно Земли (`near_moon = false`) или Луны (`true`).
    pub fn sample(&self, phase: Phase, t_s: f64, r: &Vector3, v: &Vector3, near_moon: bool) -> Result<Sample> {
        let (mr, mv) = self.moon(t_s)?;
        let (re, ve, rm) = if near_moon { (r + mr, v + mv, Some(km(r))) } else { (*r, *v, None) };
        Ok(Sample {
            t: iso(&self.epoch(t_s)),
            t_s: (t_s * 1e3).round() / 1e3,
            phase,
            r: km(&re),
            v: km_s(&ve),
            moon_r: km(&mr),
            moon_v: km_s(&mv),
            r_moon: rm,
            v_moon: near_moon.then(|| km_s(v)),
            earth_q: matches!(phase, Phase::Ascent | Phase::ParkingOrbit).then(|| self.earth_dcm(t_s).map(|m| quaternion(&m))).transpose()?,
            moon_q: near_moon.then(|| self.moon_dcm(t_s).map(|m| quaternion(&m))).transpose()?,
        })
    }

    /// Событие по состоянию относительно Земли или Луны; `dv` — импульс (км/с).
    #[allow(clippy::too_many_arguments)]
    pub fn event(
        &self,
        name: &'static str,
        description: impl Into<String>,
        t_s: f64,
        r: &Vector3,
        v: &Vector3,
        dv: Option<Vector3>,
        near_moon: bool,
    ) -> Result<Event> {
        let (mr, mv) = self.moon(t_s)?;
        let (re, ve, rm) = if near_moon { (r + mr, v + mv, Some(km(r))) } else { (*r, *v, None) };
        Ok(Event {
            name,
            description: description.into(),
            t: iso(&self.epoch(t_s)),
            t_s: (t_s * 1e3).round() / 1e3,
            r: km(&re),
            v: km_s(&ve),
            delta_v: dv.map(|d| km_s(&d)),
            moon_r: km(&mr),
            moon_v: km_s(&mv),
            r_moon: rm,
            v_moon: near_moon.then(|| km_s(v)),
            earth_q: quaternion(&self.earth_dcm(t_s)?),
            moon_q: quaternion(&self.moon_dcm(t_s)?),
            sun_dir: r3(&self.sun_dir(t_s)?, 1e9),
        })
    }
}

#[derive(Serialize)]
struct Body {
    radius_km: f64,
    mu_km3_s2: f64,
    orientation_frame: &'static str,
}

#[derive(Serialize)]
struct Site {
    lat_deg: f64,
    lon_deg: f64,
    height_km: f64,
    /// Положение в системе, связанной с телом (ITRF93/IAU_EARTH для Земли, MOON_ME для Луны), км.
    body_fixed_km: [f64; 3],
}

/// Начальная или конечная точка полёта.
#[derive(Debug, Clone, Serialize)]
pub struct Endpoint {
    /// launch — стол на космодроме; touchdown — посадка; иначе имя последнего события.
    event: &'static str,
    /// Тело, на котором (или у которого) точка: "earth" или "moon".
    body: &'static str,
    t: String,
    t_s: f64,
    /// Широта и долгота на теле, градусы; высота над опорной сферой (Луна) или эллипсоидом WGS84 (Земля), км.
    lat_deg: f64,
    lon_deg: f64,
    height_km: f64,
    /// Положение в системе, связанной с телом (ITRF93/IAU_EARTH или MOON_ME), км.
    body_fixed_km: [f64; 3],
    /// Положение относительно центра Земли, J2000, км.
    r: [f64; 3],
    /// Положение относительно центра Луны, J2000, км (для точки у Луны).
    #[serde(skip_serializing_if = "Option::is_none")]
    r_moon: Option<[f64; 3]>,
    /// Скорость относительно центра Земли, J2000, км/с.
    v: [f64; 3],
}

/// Поворот вектора кватернионом, обратным `q` (из J2000 в систему тела).
fn rotate_inverse(q: &[f64; 4], v: &Vector3) -> Vector3 {
    let [w, x, y, z] = *q;
    let qv = Vector3::new(-x, -y, -z);
    let t = qv.cross(v) * 2.0;
    v + t * w + qv.cross(&t)
}

fn start_point(ev: &Event, p: &crate::spaceport::Spaceport) -> Endpoint {
    let r = Vector3::new(ev.r[0], ev.r[1], ev.r[2]);
    Endpoint {
        event: ev.name,
        body: "earth",
        t: ev.t.clone(),
        t_s: ev.t_s,
        lat_deg: p.lat_deg,
        lon_deg: p.lon_deg,
        height_km: p.altitude_m.unwrap_or(0.0) / 1000.0,
        body_fixed_km: km(&rotate_inverse(&ev.earth_q, &r)),
        r: ev.r,
        r_moon: None,
        v: ev.v,
    }
}

fn end_point(ev: &Event) -> Endpoint {
    let rm = ev.r_moon.map(|m| Vector3::new(m[0], m[1], m[2])).unwrap_or_else(|| {
        Vector3::new(ev.r[0] - ev.moon_r[0], ev.r[1] - ev.moon_r[1], ev.r[2] - ev.moon_r[2])
    });
    let b = rotate_inverse(&ev.moon_q, &rm);
    let n = b.norm();
    Endpoint {
        event: ev.name,
        body: "moon",
        t: ev.t.clone(),
        t_s: ev.t_s,
        lat_deg: (b.z / n).asin().to_degrees(),
        lon_deg: b.y.atan2(b.x).to_degrees(),
        height_km: n - MOON_RADIUS_KM,
        body_fixed_km: km(&b),
        r: ev.r,
        r_moon: Some(km(&rm)),
        v: ev.v,
    }
}

/// Начальная (стартовый стол) и конечная (посадка или последнее событие у Луны) точки
/// уточнённой траектории лучшего варианта отчёта, если точки записаны.
pub fn endpoints(r: &TransferReport) -> (Option<Endpoint>, Option<Endpoint>) {
    let events = r
        .solution
        .as_ref()
        .and_then(|s| s.refined.as_ref())
        .and_then(|n| n.path.as_ref())
        .map_or(&[][..], |p| &p.events[..]);
    (
        events.first().map(|e| start_point(e, &r.spaceport)),
        events.last().filter(|e| e.name != "launch").map(end_point),
    )
}

#[derive(Serialize)]
struct Entry<'a> {
    /// Тип траектории перелёта: "direct_transfer" или "hohmann".
    trajectory: &'static str,
    spaceport_id: &'a str,
    spaceport_name: &'a str,
    spaceport: Site,
    /// Начальная точка (стартовый стол) и конечная (посадка или последнее событие у Луны).
    #[serde(skip_serializing_if = "Option::is_none")]
    start: Option<Endpoint>,
    #[serde(skip_serializing_if = "Option::is_none")]
    end: Option<Endpoint>,
    /// "numerical" — точки из численного интегрирования; "none" — точек нет.
    source: &'static str,
    /// Модель подъёма: "simulation" или "kinematic".
    ascent_model: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    rocket: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    unavailable_reason: Option<String>,
    /// Периоды этапов: [фаза, начало, конец] в секундах от старта.
    phases: Vec<(Phase, f64, f64)>,
    events: &'a [Event],
    samples: &'a [Sample],
}

#[derive(Serialize)]
struct Document<'a> {
    format: &'static str,
    created_utc: String,
    /// Журнал информации о полёте, записанный тем же запуском (имя файла в той же папке).
    flight_log: String,
    frame: &'static str,
    units: &'static str,
    time: &'static str,
    quaternions: &'static str,
    earth: Body,
    moon: Body,
    #[serde(skip_serializing_if = "Option::is_none")]
    landing_site: Option<Site>,
    notes: Vec<&'static str>,
    trajectories: Vec<Entry<'a>>,
}

/// Широта, долгота, высота → декартовы координаты сферы радиуса `radius_km`.
fn body_fixed(lat_deg: f64, lon_deg: f64, radius_km: f64) -> [f64; 3] {
    let (la, lo) = (lat_deg.to_radians(), lon_deg.to_radians());
    km(&(Vector3::new(la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin()) * radius_km))
}

/// Пишет журнал траектории в `path`.
pub fn write(
    path: &std::path::Path,
    flight_log: &std::path::Path,
    reports: &[TransferReport],
    landing_site: Option<(f64, f64, f64)>,
) -> Result<()> {
    let mut mu_e = 398_600.435_436;
    let mut mu_m = crate::gravity::GravityField::moon_grail(0).map(|f| f.gm).unwrap_or(4902.8);
    let mut entries = Vec::new();
    for r in reports {
        let Some(s) = &r.solution else { continue };
        mu_e = s.mu_earth_km3_s2;
        let refined = s.refined.as_ref();
        let path = refined.and_then(|n| n.path.as_ref());
        if let Some(n) = refined {
            mu_m = n.moon_gm_km3_s2;
        }
        let reason = match (refined, path) {
            (_, Some(_)) => None,
            (Some(_), None) => Some("trajectory points were not recorded".to_string()),
            (None, _) => Some(
                s.refine_failure
                    .clone()
                    .map(|f| format!("numerical refinement failed: {f}"))
                    .or_else(|| s.refine_skipped.clone())
                    .unwrap_or_else(|| "numerical refinement is disabled ([numerical] enabled = false)".into()),
            ),
        };
        let (events, samples) = path.map_or((&[][..], &[][..]), |p| (&p.events[..], &p.samples[..]));
        let mut phases: Vec<(Phase, f64, f64)> = Vec::new();
        for smp in samples {
            match phases.last_mut() {
                Some((ph, _, end)) if *ph == smp.phase => *end = smp.t_s,
                _ => phases.push((smp.phase, smp.t_s, smp.t_s)),
            }
        }
        let p = &r.spaceport;
        let alt_km = p.altitude_m.unwrap_or(0.0) / 1000.0;
        entries.push(Entry {
            trajectory: s.trajectory,
            spaceport_id: &p.id,
            spaceport_name: &p.name,
            spaceport: Site {
                lat_deg: p.lat_deg,
                lon_deg: p.lon_deg,
                height_km: alt_km,
                // Эллипсоид WGS84 — так же, как в моделировании подъёма.
                body_fixed_km: km(&crate::trajectory::geodetic_to_ecef(p.lat_deg, p.lon_deg, alt_km)),
            },
            start: events.first().map(|e| start_point(e, p)),
            end: events.last().filter(|e| e.name != "launch").map(end_point),
            source: if path.is_some() { "numerical" } else { "none" },
            ascent_model: s.ascent.model,
            rocket: s.ascent.rocket.clone(),
            unavailable_reason: reason,
            phases,
            events,
            samples,
        });
    }
    let doc = Document {
        format: FORMAT,
        created_utc: iso(&Epoch::now().unwrap_or_else(|_| Epoch::from_gregorian_utc_at_midnight(1970, 1, 1))),
        // Только имя файла: журнал полёта лежит в той же папке.
        flight_log: flight_log.file_name().map_or_else(|| flight_log.display().to_string(), |n| n.to_string_lossy().into_owned()),
        frame: "J2000 (ICRF), origin at the Earth's centre; r_moon fields are relative to the Moon's centre",
        units: "km, km/s, seconds",
        time: "t: UTC ISO 8601; t_s: seconds since launch",
        quaternions: "[w, x, y, z], rotation from the body-fixed frame to J2000 (v_j2000 = q * v_body * q^-1)",
        earth: Body {
            radius_km: EARTH_J2_RADIUS_KM,
            mu_km3_s2: mu_e,
            orientation_frame: "ITRF93 (IAU_EARTH when high-precision Earth data is not loaded)",
        },
        moon: Body { radius_km: MOON_RADIUS_KM, mu_km3_s2: mu_m, orientation_frame: "MOON_ME (DE440)" },
        landing_site: landing_site.map(|(lat, lon, h)| Site {
            lat_deg: lat,
            lon_deg: lon,
            height_km: h,
            body_fixed_km: body_fixed(lat, lon, MOON_RADIUS_KM + h),
        }),
        notes: vec![
            "Ascent: when the rocket has stage data, a numerical simulation (stages, thrust vs ambient pressure, US Standard Atmosphere 1976 drag, Earth J2 and rotation, gravity turn then closed-loop guidance); otherwise a kinematic profile (altitude h = H*(6t^2 - 8t^3 + 3t^4), downrange angle cubic, duration and downrange from transfer.toml). See ascent_model in each trajectory.",
            "Points come from numerical integration: Earth with J2, Moon, Sun during the transfer; the Moon's gravity field with mascons (GRAIL), Earth and Sun in lunar orbit.",
            "TLI: when the engine is known (rocket upper stage or the spacecraft), the burn is simulated with constant thrust and decreasing mass between events tli and tli_end (phase tli_burn); otherwise it is impulsive at event tli. LOI and DOI are impulsive: the velocity jumps at the event.",
            "Speed near the Moon: use v_moon (relative to the Moon's centre). v is relative to the Earth's centre and swings between ~0.6 and ~2.6 km/s every lunar orbit, because the Moon itself moves at ~1 km/s around the Earth. moon_v is the Moon's own velocity (v = moon_v + v_moon).",
            "On the surface (terminal phase) v_moon is not zero: the point turns with the Moon, ~4.6 m/s at the lunar equator in J2000. Speed over the ground is v_moon minus that rotation.",
            "Braking is simulated with a point-mass Moon and the Moon's rotation; the terminal phase is a time allowance at the site, not a simulation.",
            "Sampling: 1 s during a simulated ascent (10 s for the kinematic profile), 30 s in the parking orbit, 10 s during the TLI burn, 300 s in the transfer (30 s within 70000 km of the Moon), 60 s in lunar orbit, 10 s on the descent ellipse, 5 s while braking.",
        ],
        trajectories: entries,
    };
    let text = serde_json::to_string(&doc).map_err(|e| Error::Config(format!("trajectory log: {e}")))?;
    std::fs::write(path, text).map_err(|source| Error::Io { path: path.to_path_buf(), source })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quaternion_round_trip() {
        // Поворот на 90° вокруг z.
        let m = Matrix3::new(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let q = quaternion(&m);
        let h = std::f64::consts::FRAC_1_SQRT_2;
        assert!((q[0] - h).abs() < 1e-12 && q[1].abs() < 1e-12 && q[2].abs() < 1e-12 && (q[3] - h).abs() < 1e-12);
        // Произвольный поворот: восстановленная матрица совпадает.
        let a = anise::math::rotation::r1(0.3) * anise::math::rotation::r3(-1.2) * anise::math::rotation::r2(2.5);
        let [w, x, y, z] = quaternion(&a);
        let back = Matrix3::new(
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - w * z),
            2.0 * (x * z + w * y),
            2.0 * (x * y + w * z),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - w * x),
            2.0 * (x * z - w * y),
            2.0 * (y * z + w * x),
            1.0 - 2.0 * (x * x + y * y),
        );
        assert!((back - a).norm() < 1e-9, "{}", (back - a).norm());
    }

    #[test]
    fn rounding() {
        assert_eq!(km(&Vector3::new(1.23456789, -2.0004, 0.0)), [1.235, -2.0, 0.0]);
        assert_eq!(body_fixed(0.0, 90.0, 10.0), [0.0, 10.0, 0.0]);
    }
}
