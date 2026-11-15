//! Численное интегрирование движения: метод Дормана — Принса 5(4) с автоматическим шагом
//! и остановкой по событию, и модели сил для перелёта и окололунной орбиты.
//!
//! Состояние — `[x, y, z, vx, vy, vz]` в км и км/с; время — секунды от опорной эпохи.

use anise::constants::frames::{EARTH_J2000, IAU_EARTH_FRAME, MOON_J2000, MOON_PA_DE440_FRAME, SUN_J2000};
use anise::math::{Matrix3, Vector3};
use anise::prelude::{Almanac, Frame};
use anise::time::{Epoch, Unit};

use crate::error::{Error, Result};
use crate::gravity::GravityField;

pub type State = [f64; 6];

/// Сжатие Земли J2 и её экваториальный радиус (EGM2008), км.
pub const EARTH_J2: f64 = 1.082_626_68e-3;
pub const EARTH_J2_RADIUS_KM: f64 = 6378.1363;
const MU_SUN_DEFAULT: f64 = 1.327_124_400_41e11;

pub fn pos(y: &State) -> Vector3 {
    Vector3::new(y[0], y[1], y[2])
}

pub fn vel(y: &State) -> Vector3 {
    Vector3::new(y[3], y[4], y[5])
}

pub fn state(r: &Vector3, v: &Vector3) -> State {
    [r.x, r.y, r.z, v.x, v.y, v.z]
}

/// Правая часть: ускорение (км/с²) по времени и состоянию.
pub trait Dynamics: Sync {
    fn acceleration(&self, t: f64, y: &State) -> Result<Vector3>;
}

fn deriv(d: &dyn Dynamics, t: f64, y: &State) -> Result<State> {
    let a = d.acceleration(t, y)?;
    Ok([y[3], y[4], y[5], a.x, a.y, a.z])
}

/// Один шаг Дормана — Принса: новое состояние и оценка относительной ошибки.
fn dopri_step(d: &dyn Dynamics, t: f64, y: &State, h: f64) -> Result<(State, f64)> {
    const A: [[f64; 6]; 6] = [
        [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
        [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
        [19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0],
        [9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0],
        [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0],
    ];
    const C: [f64; 6] = [1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];
    const B4: [f64; 7] =
        [5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0];
    let mut k = [[0.0f64; 6]; 7];
    k[0] = deriv(d, t, y)?;
    let mut y5 = *y;
    for s in 0..6 {
        let mut ys = *y;
        for (j, kj) in k.iter().enumerate().take(s + 1) {
            for i in 0..6 {
                ys[i] += h * A[s][j] * kj[i];
            }
        }
        k[s + 1] = deriv(d, t + C[s] * h, &ys)?;
        if s == 5 {
            y5 = ys;
        }
    }
    let mut err = [0.0f64; 6];
    for i in 0..6 {
        let y4: f64 = (0..7).map(|j| B4[j] * k[j][i]).sum::<f64>() * h + y[i];
        err[i] = y5[i] - y4;
    }
    let scale_r = pos(&y5).norm().max(1.0);
    let scale_v = vel(&y5).norm().max(1e-3);
    let e = (Vector3::new(err[0], err[1], err[2]).norm() / scale_r)
        .max(Vector3::new(err[3], err[4], err[5]).norm() / scale_v);
    Ok((y5, e))
}

/// Чем закончилось интегрирование.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Outcome {
    pub t: f64,
    pub y: State,
    /// `true` — остановка по событию, `false` — дошли до `t_end`.
    pub event: bool,
    pub steps: u32,
}

/// Интегрирует от `t0` до `t_end` (вперёд) с относительной точностью `tol`.
/// `stop(t, y)` — функция события: интегрирование останавливается, когда она переходит
/// от положительных значений к неположительным. `observe` вызывается после каждого шага.
pub fn integrate(
    d: &dyn Dynamics,
    t0: f64,
    y0: &State,
    t_end: f64,
    tol: f64,
    mut stop: Option<&mut dyn FnMut(f64, &State) -> Result<f64>>,
    mut observe: Option<&mut dyn FnMut(f64, &State)>,
) -> Result<Outcome> {
    let (mut t, mut y) = (t0, *y0);
    let mut g_prev = match stop.as_mut() {
        Some(g) => g(t, &y)?,
        None => 1.0,
    };
    // Начальный шаг — доля характерного времени r/v.
    let mut h = (0.01 * pos(&y).norm() / vel(&y).norm().max(1e-6)).clamp(0.1, 600.0);
    let mut steps = 0u32;
    while t < t_end {
        if steps > 2_000_000 {
            return Err(Error::Config("numerical integration: too many steps".into()));
        }
        let h_try = h.min(t_end - t);
        let (y_new, err) = dopri_step(d, t, &y, h_try)?;
        if !err.is_finite() {
            return Err(Error::Config("numerical integration diverged".into()));
        }
        if err > tol {
            h = h_try * (0.9 * (tol / err).powf(0.2)).clamp(0.1, 0.9);
            continue;
        }
        steps += 1;
        let t_new = t + h_try;
        if let Some(g) = stop.as_mut() {
            let g_new = g(t_new, &y_new)?;
            if g_prev > 0.0 && g_new <= 0.0 {
                // Уточняем момент события делением шага (метод ложного положения с подстраховкой).
                let (mut lo, mut hi, mut g_lo, mut g_hi) = (0.0, h_try, g_prev, g_new);
                let mut best = (t_new, y_new);
                for _ in 0..60 {
                    let mid = if g_lo - g_hi != 0.0 { lo + (hi - lo) * g_lo / (g_lo - g_hi) } else { 0.5 * (lo + hi) };
                    let mid = mid.clamp(lo + 0.05 * (hi - lo), hi - 0.05 * (hi - lo));
                    let (ym, _) = dopri_step(d, t, &y, mid)?;
                    let gm = g(t + mid, &ym)?;
                    if gm > 0.0 {
                        (lo, g_lo) = (mid, gm);
                    } else {
                        (hi, g_hi) = (mid, gm);
                        best = (t + mid, ym);
                    }
                    if hi - lo < 1e-4 {
                        break;
                    }
                }
                if let Some(o) = observe.as_mut() {
                    o(best.0, &best.1);
                }
                return Ok(Outcome { t: best.0, y: best.1, event: true, steps });
            }
            g_prev = g_new;
        }
        (t, y) = (t_new, y_new);
        if let Some(o) = observe.as_mut() {
            o(t, &y);
        }
        h = h_try * (0.9 * (tol / err.max(1e-300)).powf(0.2)).clamp(0.2, 5.0);
    }
    Ok(Outcome { t, y, event: false, steps })
}

fn mu_of(almanac: &Almanac, frame: Frame, default: f64) -> f64 {
    almanac.frame_info(frame).ok().and_then(|f| f.mu_km3_s2().ok()).unwrap_or(default)
}

fn body(almanac: &Almanac, target: Frame, observer: Frame, epoch: Epoch) -> Result<(Vector3, Vector3)> {
    let s = almanac
        .translate(target, observer, epoch, None)
        .map_err(|e| Error::Almanac(format!("{target} from {observer} at {epoch}: {e}")))?;
    Ok((s.radius_km, s.velocity_km_s))
}

/// Шаг узлов таблиц положений и поворотов, с.
const TRACK_STEP_S: f64 = 600.0;

/// Положение тела относительно другого: узлы из ANISE с шагом 10 мин (положение и скорость),
/// между узлами — кубическая интерполяция Эрмита. Ошибка для Луны около 1e-8 км (порядка
/// микрометров), для Солнца ещё меньше; запрос в ~30 раз быстрее, чем к ANISE.
/// Узлы считаются по мере надобности и хранятся.
pub struct BodyTrack<'a> {
    almanac: &'a Almanac,
    t0: Epoch,
    target: Frame,
    observer: Frame,
    nodes: std::sync::Mutex<std::collections::HashMap<i64, (Vector3, Vector3)>>,
}

impl<'a> BodyTrack<'a> {
    pub fn new(almanac: &'a Almanac, t0: Epoch, target: Frame, observer: Frame) -> Self {
        Self { almanac, t0, target, observer, nodes: Default::default() }
    }

    fn node(&self, k: i64) -> Result<(Vector3, Vector3)> {
        if let Some(n) = self.nodes.lock().unwrap_or_else(|e| e.into_inner()).get(&k) {
            return Ok(*n);
        }
        let n = body(self.almanac, self.target, self.observer, self.t0 + (k as f64 * TRACK_STEP_S) * Unit::Second)?;
        self.nodes.lock().unwrap_or_else(|e| e.into_inner()).insert(k, n);
        Ok(n)
    }

    /// Положение (км) и скорость (км/с) в момент `t` секунд от `t0`.
    pub fn at(&self, t: f64) -> Result<(Vector3, Vector3)> {
        let h = TRACK_STEP_S;
        let x = t / h;
        let k = x.floor();
        let s = x - k;
        let ((p0, v0), (p1, v1)) = (self.node(k as i64)?, self.node(k as i64 + 1)?);
        let (s2, s3) = (s * s, s * s * s);
        let p = p0 * (2.0 * s3 - 3.0 * s2 + 1.0)
            + v0 * (h * (s3 - 2.0 * s2 + s))
            + p1 * (-2.0 * s3 + 3.0 * s2)
            + v1 * (h * (s3 - s2));
        let v = (p0 * (6.0 * s2 - 6.0 * s) + p1 * (6.0 * s - 6.0 * s2)) / h
            + v0 * (3.0 * s2 - 4.0 * s + 1.0)
            + v1 * (3.0 * s2 - 2.0 * s);
        Ok((p, v))
    }
}

/// Поворот главных осей Луны в J2000: узлы из ANISE с шагом 10 мин, между ними — кубическая
/// интерполяция по четырём узлам (ошибка элементов матрицы ~1e-13).
struct RotationTrack<'a> {
    almanac: &'a Almanac,
    t0: Epoch,
    nodes: std::sync::Mutex<std::collections::HashMap<i64, Matrix3>>,
}

impl RotationTrack<'_> {
    fn node(&self, k: i64) -> Result<Matrix3> {
        if let Some(m) = self.nodes.lock().unwrap_or_else(|e| e.into_inner()).get(&k) {
            return Ok(*m);
        }
        let e = self.t0 + (k as f64 * TRACK_STEP_S) * Unit::Second;
        let m = self
            .almanac
            .rotate(MOON_PA_DE440_FRAME, MOON_J2000, e)
            .map_err(|err| Error::Almanac(format!("Moon principal axes at {e}: {err}")))?
            .rot_mat;
        self.nodes.lock().unwrap_or_else(|e| e.into_inner()).insert(k, m);
        Ok(m)
    }

    fn at(&self, t: f64) -> Result<Matrix3> {
        let x = t / TRACK_STEP_S;
        let k = x.floor() as i64;
        let s = x - k as f64;
        // Лагранж по узлам k-1, k, k+1, k+2 (s отсчитывается от узла k).
        let w = [
            -s * (s - 1.0) * (s - 2.0) / 6.0,
            (s + 1.0) * (s - 1.0) * (s - 2.0) / 2.0,
            -(s + 1.0) * s * (s - 2.0) / 2.0,
            (s + 1.0) * s * (s - 1.0) / 6.0,
        ];
        let mut m = Matrix3::zeros();
        for (i, wi) in w.iter().enumerate() {
            m += self.node(k - 1 + i as i64)? * *wi;
        }
        Ok(m)
    }
}

/// Притяжение третьего тела в системе с началом в центральном теле.
fn third_body(mu: f64, r: &Vector3, body_pos: &Vector3) -> Vector3 {
    let d = body_pos - r;
    d * (mu / d.norm().powi(3)) - body_pos * (mu / body_pos.norm().powi(3))
}

/// Какие силы учитывать.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Forces {
    pub earth_j2: bool,
    pub sun: bool,
}

/// Перелёт: система с началом в центре Земли (J2000). Земля (со сжатием J2), Луна, Солнце.
pub struct EarthCentered<'a> {
    pub almanac: &'a Almanac,
    pub t0: Epoch,
    pub mu_earth: f64,
    pub mu_moon: f64,
    pub mu_sun: f64,
    pub forces: Forces,
    /// Ось вращения Земли в J2000.
    pub pole: Vector3,
    /// Таблицы положений Луны и Солнца для правой части (интерполяция между узлами).
    moon_track: BodyTrack<'a>,
    sun_track: BodyTrack<'a>,
}

impl<'a> EarthCentered<'a> {
    pub fn new(almanac: &'a Almanac, t0: Epoch, mu_earth: f64, mu_moon: f64, forces: Forces) -> Result<Self> {
        let dcm = almanac
            .rotate(IAU_EARTH_FRAME, EARTH_J2000, t0)
            .map_err(|e| Error::Almanac(format!("Earth orientation at {t0}: {e}")))?;
        Ok(Self {
            almanac,
            t0,
            mu_earth,
            mu_moon,
            mu_sun: mu_of(almanac, SUN_J2000, MU_SUN_DEFAULT),
            forces,
            pole: dcm * Vector3::new(0.0, 0.0, 1.0),
            moon_track: BodyTrack::new(almanac, t0, MOON_J2000, EARTH_J2000),
            sun_track: BodyTrack::new(almanac, t0, SUN_J2000, EARTH_J2000),
        })
    }

    pub fn epoch(&self, t: f64) -> Epoch {
        self.t0 + t * Unit::Second
    }

    /// Положение и скорость Луны относительно Земли.
    pub fn moon(&self, t: f64) -> Result<(Vector3, Vector3)> {
        body(self.almanac, MOON_J2000, EARTH_J2000, self.epoch(t))
    }
}

impl Dynamics for EarthCentered<'_> {
    fn acceleration(&self, t: f64, y: &State) -> Result<Vector3> {
        let r = pos(y);
        let rn = r.norm();
        let mut a = -r * (self.mu_earth / rn.powi(3));
        if self.forces.earth_j2 {
            let z = r.dot(&self.pole);
            let k = 1.5 * EARTH_J2 * self.mu_earth * EARTH_J2_RADIUS_KM.powi(2) / rn.powi(5);
            a += (r * (5.0 * z * z / (rn * rn) - 1.0) - self.pole * (2.0 * z)) * k;
        }
        a += third_body(self.mu_moon, &r, &self.moon_track.at(t)?.0);
        if self.forces.sun {
            a += third_body(self.mu_sun, &r, &self.sun_track.at(t)?.0);
        }
        Ok(a)
    }
}

/// Окололунная орбита: система с началом в центре Луны (J2000). Поле Луны с масконами,
/// притяжение Земли и Солнца.
pub struct MoonCentered<'a> {
    pub almanac: &'a Almanac,
    pub t0: Epoch,
    pub field: &'a GravityField,
    pub mu_earth: f64,
    pub mu_sun: f64,
    pub forces: Forces,
    earth_track: BodyTrack<'a>,
    sun_track: BodyTrack<'a>,
    rotation: RotationTrack<'a>,
}

impl<'a> MoonCentered<'a> {
    pub fn new(almanac: &'a Almanac, t0: Epoch, field: &'a GravityField, mu_earth: f64, forces: Forces) -> Self {
        Self {
            almanac,
            t0,
            field,
            mu_earth,
            mu_sun: mu_of(almanac, SUN_J2000, MU_SUN_DEFAULT),
            forces,
            earth_track: BodyTrack::new(almanac, t0, EARTH_J2000, MOON_J2000),
            sun_track: BodyTrack::new(almanac, t0, SUN_J2000, MOON_J2000),
            rotation: RotationTrack { almanac, t0, nodes: Default::default() },
        }
    }

    pub fn epoch(&self, t: f64) -> Epoch {
        self.t0 + t * Unit::Second
    }

    /// Матрица перехода из главных осей Луны в J2000.
    pub fn pa_to_j2000(&self, t: f64) -> Result<Matrix3> {
        let e = self.epoch(t);
        Ok(self
            .almanac
            .rotate(MOON_PA_DE440_FRAME, MOON_J2000, e)
            .map_err(|err| Error::Almanac(format!("Moon principal axes at {e}: {err}")))?
            .rot_mat)
    }
}

impl Dynamics for MoonCentered<'_> {
    fn acceleration(&self, t: f64, y: &State) -> Result<Vector3> {
        let r = pos(y);
        let mut a = if self.field.degree == 0 {
            self.field.acceleration(&r)
        } else {
            let m = self.rotation.at(t)?;
            m * self.field.acceleration(&(m.transpose() * r))
        };
        a += third_body(self.mu_earth, &r, &self.earth_track.at(t)?.0);
        if self.forces.sun {
            a += third_body(self.mu_sun, &r, &self.sun_track.at(t)?.0);
        }
        Ok(a)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore = "needs data/kernels"]
    fn real_tracks_match_anise() {
        let almanac = crate::ephemeris::load_almanac_blocking(crate::ephemeris::AlmanacSource::LocalDir("data/kernels".into()))
            .unwrap();
        let t0 = Epoch::from_gregorian_utc_at_midnight(2027, 1, 1);
        let moon = BodyTrack::new(&almanac, t0, MOON_J2000, EARTH_J2000);
        let sun = BodyTrack::new(&almanac, t0, SUN_J2000, MOON_J2000);
        let rot = RotationTrack { almanac: &almanac, t0, nodes: Default::default() };
        for t in [-1234.5, 0.0, 17.0, 299.9, 86_400.0 * 2.37, 86_400.0 * 5.1] {
            let e = t0 + t * Unit::Second;
            let (p, v) = moon.at(t).unwrap();
            let (pe, ve) = body(&almanac, MOON_J2000, EARTH_J2000, e).unwrap();
            assert!((p - pe).norm() < 1e-6 && (v - ve).norm() < 1e-9, "moon {t}: {} km", (p - pe).norm());
            let (p, _) = sun.at(t).unwrap();
            let (pe, _) = body(&almanac, SUN_J2000, MOON_J2000, e).unwrap();
            assert!((p - pe).norm() < 1e-5, "sun {t}: {} km", (p - pe).norm());
            let m = rot.at(t).unwrap();
            let me = almanac.rotate(MOON_PA_DE440_FRAME, MOON_J2000, e).unwrap().rot_mat;
            assert!((m - me).norm() < 1e-11, "rotation {t}: {}", (m - me).norm());
        }
    }

    struct TwoBody(f64);
    impl Dynamics for TwoBody {
        fn acceleration(&self, _t: f64, y: &State) -> Result<Vector3> {
            let r = pos(y);
            Ok(-r * (self.0 / r.norm().powi(3)))
        }
    }

    #[test]
    fn kepler_orbit_closes_and_conserves_energy() {
        let mu = 398_600.4418_f64;
        // Эллипс 6600 x 40000 км: после одного периода возвращаемся в исходную точку.
        let (rp, ra) = (6600.0_f64, 40000.0_f64);
        let a = 0.5 * (rp + ra);
        let vp = (mu * (2.0 / rp - 1.0 / a)).sqrt();
        let y0 = state(&Vector3::new(rp, 0.0, 0.0), &Vector3::new(0.0, vp * 0.8, vp * 0.6));
        let period = std::f64::consts::TAU * (a.powi(3) / mu).sqrt();
        let out = integrate(&TwoBody(mu), 0.0, &y0, period, 1e-11, None, None).unwrap();
        assert!(!out.event);
        assert!((pos(&out.y) - pos(&y0)).norm() < 1e-3, "{}", (pos(&out.y) - pos(&y0)).norm());
        let energy = |y: &State| vel(y).norm_squared() / 2.0 - mu / pos(y).norm();
        assert!((energy(&out.y) - energy(&y0)).abs() < 1e-8);
    }

    #[test]
    fn stops_at_event() {
        let mu = 398_600.4418_f64;
        let (rp, ra) = (6600.0_f64, 40000.0_f64);
        let a = 0.5 * (rp + ra);
        let vp = (mu * (2.0 / rp - 1.0 / a)).sqrt();
        let y0 = state(&Vector3::new(rp, 0.0, 0.0), &Vector3::new(0.0, vp, 0.0));
        // Апоцентр: радиальная скорость меняет знак с плюса на минус.
        let mut apo = |_t: f64, y: &State| Ok(pos(y).dot(&vel(y)).max(-1.0).min(1.0) + if pos(y).norm() < 7000.0 { 1.0 } else { 0.0 });
        let out = integrate(&TwoBody(mu), 0.0, &y0, 1e6, 1e-11, Some(&mut apo), None).unwrap();
        assert!(out.event);
        let half = std::f64::consts::PI * (a.powi(3) / mu).sqrt();
        assert!((out.t - half).abs() < 0.01, "{} vs {half}", out.t);
        assert!((pos(&out.y).norm() - ra).abs() < 1e-3);
        // Наблюдатель видит каждый шаг.
        let mut n = 0u32;
        let mut obs = |_t: f64, _y: &State| n += 1;
        let out = integrate(&TwoBody(mu), 0.0, &y0, 1000.0, 1e-9, None, Some(&mut obs)).unwrap();
        assert_eq!(n, out.steps);
    }
}
