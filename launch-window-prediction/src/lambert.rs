//! Задача Ламберта (однооборотная): по двум положениям и времени перелёта найти скорости.
//!
//! Метод универсальных переменных (H. Curtis, «Orbital Mechanics for Engineering Students»,
//! алгоритм 5.2). Направление движения задаётся нормалью плоскости орбиты `normal`:
//! угол перелёта отсчитывается от `r1` к `r2` против часовой стрелки, если смотреть с конца
//! `normal`. Поэтому поддерживаются углы и меньше, и больше 180 градусов.

use anise::math::Vector3;

use crate::error::{Error, Result};

/// Решение задачи Ламберта.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LambertSolution {
    /// Скорость в начальной точке, км/с.
    pub v1: Vector3,
    /// Скорость в конечной точке, км/с.
    pub v2: Vector3,
    /// Угол перелёта, радианы (0, 2π).
    pub transfer_angle_rad: f64,
}

/// Функции Штумпфа C(z) и S(z).
fn stumpff(z: f64) -> (f64, f64) {
    if z.abs() < 1e-3 {
        // Ряды Тейлора: без потери точности около нуля.
        let c = 0.5 - z / 24.0 + z * z / 720.0 - z * z * z / 40_320.0;
        let s = 1.0 / 6.0 - z / 120.0 + z * z / 5_040.0 - z * z * z / 362_880.0;
        (c, s)
    } else if z > 0.0 {
        let sz = z.sqrt();
        ((1.0 - sz.cos()) / z, (sz - sz.sin()) / (sz * sz * sz))
    } else {
        let sz = (-z).sqrt();
        ((sz.cosh() - 1.0) / (-z), (sz.sinh() - sz) / (sz * sz * sz))
    }
}

/// Решает задачу Ламберта.
///
/// - `r1`, `r2` — положения, км; `tof_s` — время перелёта, с; `mu` — гравитационный параметр, км³/с²;
/// - `normal` — нормаль плоскости, задающая направление движения (не обязательно единичная).
pub fn solve(r1: &Vector3, r2: &Vector3, tof_s: f64, mu: f64, normal: &Vector3) -> Result<LambertSolution> {
    let (n1, n2) = (r1.norm(), r2.norm());
    if !(n1 > 0.0 && n2 > 0.0 && tof_s > 0.0 && mu > 0.0) {
        return Err(Error::Config("lambert: invalid input".into()));
    }
    let cos_dth = (r1.dot(r2) / (n1 * n2)).clamp(-1.0, 1.0);
    let sin_sign = r1.cross(r2).dot(normal);
    let mut dth = cos_dth.acos();
    if sin_sign < 0.0 {
        dth = std::f64::consts::TAU - dth;
    }
    // Вырожденные случаи: точки на одной прямой с центром — плоскость перелёта не определена.
    if dth < 1e-6 || (dth - std::f64::consts::PI).abs() < 1e-6 || dth > std::f64::consts::TAU - 1e-6 {
        return Err(Error::Config("lambert: transfer angle is degenerate (0 or 180 deg)".into()));
    }

    let a = dth.sin() * (n1 * n2 / (1.0 - cos_dth)).sqrt();
    let sqrt_mu = mu.sqrt();
    let y_of = |z: f64| {
        let (c, s) = stumpff(z);
        (n1 + n2 + a * (z * s - 1.0) / c.sqrt(), c, s)
    };
    // F(z) монотонно растёт с z; там, где y <= 0, перелёт заведомо «слишком быстрый».
    let f_of = |z: f64| -> f64 {
        let (y, c, s) = y_of(z);
        if y <= 0.0 {
            return f64::NEG_INFINITY;
        }
        (y / c).powf(1.5) * s + a * y.sqrt() - sqrt_mu * tof_s
    };

    // У 4π² функции Штумпфа вырождаются (1 - cos √z теряет точность), поэтому граница чуть ниже.
    let z_max = 4.0 * std::f64::consts::PI * std::f64::consts::PI * (1.0 - 1e-4);
    let mut lo = -100.0;
    let mut tries = 0;
    while f_of(lo) > 0.0 {
        lo *= 2.0;
        tries += 1;
        if tries > 40 {
            return Err(Error::Config("lambert: no solution (time of flight too short)".into()));
        }
    }
    let mut hi = z_max;
    if !(f_of(hi) >= 0.0) {
        return Err(Error::Config("lambert: no single-revolution solution (time of flight too long)".into()));
    }
    // F и её производная по z (Curtis, формула 5.43) за одно вычисление функций Штумпфа;
    // около z = 0 формула производной теряет точность — там производной нет.
    let f_df = |z: f64| -> (f64, Option<f64>) {
        let (y, c, s) = y_of(z);
        if y <= 0.0 {
            return (f64::NEG_INFINITY, None);
        }
        let yc = (y / c).powf(1.5);
        let f = yc * s + a * y.sqrt() - sqrt_mu * tof_s;
        if z.abs() < 1e-3 {
            return (f, None);
        }
        let d = yc * ((c - 1.5 * s / c) / (2.0 * z) + 0.75 * s * s / c)
            + a / 8.0 * (3.0 * s / c * y.sqrt() + a * (c / y).sqrt());
        (f, (d.is_finite() && d > 0.0).then_some(d))
    };
    // Метод Ньютона внутри интервала [lo, hi], где F меняет знак; шаг вне интервала
    // (или без производной) заменяется делением пополам. Сходится за 5-10 шагов вместо ~45.
    // Начало — чуть выше z = 0 (почти параболическая орбита, типичная для перелёта к Луне).
    let mut z = 0.5f64.clamp(lo, hi);
    for _ in 0..200 {
        let (fz, dfz) = f_df(z);
        if fz < 0.0 {
            lo = z;
        } else {
            hi = z;
        }
        let newton = dfz.map(|d| z - fz / d).filter(|zn| *zn > lo && *zn < hi);
        let next = newton.unwrap_or(0.5 * (lo + hi));
        let done = (next - z).abs() < 1e-12 * (1.0 + z.abs()) || hi - lo < 1e-12 * (1.0 + hi.abs());
        z = next;
        if done {
            break;
        }
    }
    let (y, _, _) = y_of(z);
    if !(y > 0.0) {
        return Err(Error::Config("lambert: did not converge".into()));
    }
    let f = 1.0 - y / n1;
    let g = a * (y / mu).sqrt();
    let gdot = 1.0 - y / n2;
    Ok(LambertSolution {
        v1: (r2 - r1 * f) / g,
        v2: (r2 * gdot - r1) / g,
        transfer_angle_rad: dth,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    const MU_EARTH: f64 = 398_600.0;
    const Z: Vector3 = Vector3::new(0.0, 0.0, 1.0);

    #[test]
    fn curtis_example_5_2() {
        let r1 = Vector3::new(5000.0, 10000.0, 2100.0);
        let r2 = Vector3::new(-14600.0, 2500.0, 7000.0);
        let s = solve(&r1, &r2, 3600.0, MU_EARTH, &r1.cross(&r2)).unwrap();
        assert_abs_diff_eq!(s.v1.x, -5.9925, epsilon = 1e-3);
        assert_abs_diff_eq!(s.v1.y, 1.9254, epsilon = 1e-3);
        assert_abs_diff_eq!(s.v1.z, 3.2456, epsilon = 1e-3);
        assert_abs_diff_eq!(s.v2.x, -3.3125, epsilon = 1e-3);
        assert_abs_diff_eq!(s.v2.y, -4.1966, epsilon = 1e-3);
        assert_abs_diff_eq!(s.v2.z, -0.38529, epsilon = 1e-3);
    }

    #[test]
    fn circular_orbit_short_and_long_way() {
        let r = 7000.0;
        let v_circ = (MU_EARTH / r).sqrt();
        let n = (MU_EARTH / r.powi(3)).sqrt();
        for theta_deg in [60.0_f64, 150.0, 210.0, 300.0] {
            let th = theta_deg.to_radians();
            let r1 = Vector3::new(r, 0.0, 0.0);
            let r2 = Vector3::new(r * th.cos(), r * th.sin(), 0.0);
            let s = solve(&r1, &r2, th / n, MU_EARTH, &Z).unwrap_or_else(|e| panic!("{theta_deg}: {e}"));
            assert_abs_diff_eq!(s.transfer_angle_rad, th, epsilon = 1e-9);
            assert_abs_diff_eq!(s.v1.x, 0.0, epsilon = 1e-6);
            assert_abs_diff_eq!(s.v1.y, v_circ, epsilon = 1e-6);
            assert_abs_diff_eq!(s.v2.norm(), v_circ, epsilon = 1e-6);
        }
    }

    #[test]
    fn hohmann_like_transfer_energy() {
        // Почти полвитка от 6563 км до 384400 км: скорость в перигее близка к гомановской (~10.92 км/с).
        let (rp, ra) = (6563.0_f64, 384_400.0_f64);
        let a = 0.5 * (rp + ra);
        let half_period = std::f64::consts::PI * (a.powi(3) / MU_EARTH).sqrt();
        let th = 179.0_f64.to_radians();
        let r1 = Vector3::new(rp, 0.0, 0.0);
        let r2 = Vector3::new(ra * th.cos(), ra * th.sin(), 0.0);
        let s = solve(&r1, &r2, half_period, MU_EARTH, &Z).unwrap();
        let v_hohmann = (MU_EARTH * (2.0 / rp - 1.0 / a)).sqrt();
        assert_abs_diff_eq!(s.v1.norm(), v_hohmann, epsilon = 0.01);
    }

    #[test]
    fn energy_and_momentum_conserved() {
        let r1 = Vector3::new(7000.0, 1000.0, 500.0);
        let r2 = Vector3::new(-20000.0, 30000.0, 4000.0);
        let s = solve(&r1, &r2, 5.0 * 3600.0, MU_EARTH, &r1.cross(&r2)).unwrap();
        let e1 = s.v1.norm_squared() / 2.0 - MU_EARTH / r1.norm();
        let e2 = s.v2.norm_squared() / 2.0 - MU_EARTH / r2.norm();
        assert_abs_diff_eq!(e1, e2, epsilon = 1e-9);
        let (h1, h2) = (r1.cross(&s.v1), r2.cross(&s.v2));
        assert_abs_diff_eq!((h1 - h2).norm(), 0.0, epsilon = 1e-6);
    }

    #[test]
    fn rejects_degenerate_geometry() {
        let r1 = Vector3::new(7000.0, 0.0, 0.0);
        let r2 = Vector3::new(-9000.0, 0.0, 0.0);
        assert!(solve(&r1, &r2, 3600.0, MU_EARTH, &Z).is_err());
        assert!(solve(&r1, &r2, -1.0, MU_EARTH, &Z).is_err());
    }
}
