//! Поле тяготения Луны с неоднородностями (масконами) по данным миссии GRAIL.
//!
//! Поле задаётся разложением по сферическим гармоникам. Коэффициенты — модель GL0660B
//! (JPL, файл `jggrx_0660b_sha.tab` из архива PDS Geosciences), обрезанная до степени 60:
//! `data/gravity/grail_gl0660b_degree60.tab` (встроена в программу). Оси — главные оси
//! инерции Луны (кадр `MOON_PA_DE440`).
//!
//! Ускорение считается рекурсией Каннингема (Montenbruck & Gill, «Satellite Orbits», 3.2.4-3.2.5):
//! она не имеет особенностей на полюсах, что важно для полярных орбит.

use anise::math::Vector3;

use crate::error::{Error, Result};

/// Встроенные коэффициенты GRAIL GL0660B до степени 60.
const GL0660B_DEGREE_60: &str = include_str!("../data/gravity/grail_gl0660b_degree60.tab");
/// Наибольшая поддерживаемая степень (выше ненормированная рекурсия теряет точность).
pub const MAX_DEGREE: usize = 100;

#[derive(Debug, Clone, PartialEq)]
pub struct GravityField {
    /// Гравитационный параметр, км³/с².
    pub gm: f64,
    /// Опорный радиус разложения, км.
    pub radius_km: f64,
    pub degree: usize,
    /// Ненормированные коэффициенты, индекс `n(n+1)/2 + m`.
    c: Vec<f64>,
    s: Vec<f64>,
}

fn idx(n: usize, m: usize) -> usize {
    n * (n + 1) / 2 + m
}

impl GravityField {
    /// Точечная масса (без неоднородностей).
    pub fn point_mass(gm: f64, radius_km: f64) -> Self {
        Self { gm, radius_km, degree: 0, c: vec![1.0], s: vec![0.0] }
    }

    /// Поле Луны GRAIL GL0660B до степени `degree` (0 — точечная масса с тем же GM).
    pub fn moon_grail(degree: usize) -> Result<Self> {
        Self::parse(GL0660B_DEGREE_60, degree)
    }

    /// Разбирает таблицу формата PDS SHADR: первая строка — опорный радиус (км) и GM,
    /// далее строки `n, m, C, S` с полностью нормированными коэффициентами.
    pub fn parse(text: &str, degree: usize) -> Result<Self> {
        let bad = |m: String| Error::Config(format!("gravity field: {m}"));
        if degree > MAX_DEGREE {
            return Err(bad(format!("degree {degree} is above the supported maximum {MAX_DEGREE}")));
        }
        let mut lines = text.lines().filter(|l| !l.trim().is_empty());
        let head: Vec<f64> = lines
            .next()
            .ok_or_else(|| bad("empty table".into()))?
            .split(',')
            .take(2)
            .map(|x| x.trim().parse::<f64>())
            .collect::<std::result::Result<_, _>>()
            .map_err(|e| bad(format!("header: {e}")))?;
        let [radius_km, gm] = head[..] else { return Err(bad("header needs radius and GM".into())) };
        let size = idx(degree, degree) + 1;
        let (mut c, mut s) = (vec![0.0; size], vec![0.0; size]);
        c[0] = 1.0;
        let mut max_seen = 0usize;
        for line in lines {
            let p: Vec<&str> = line.split(',').map(str::trim).collect();
            if p.len() < 4 {
                return Err(bad(format!("bad line: {line}")));
            }
            let num = |x: &str| x.parse::<f64>().map_err(|e| bad(format!("{x}: {e}")));
            let (n, m) = (num(p[0])? as usize, num(p[1])? as usize);
            max_seen = max_seen.max(n);
            if n > degree || m > n {
                continue;
            }
            // Нормировка: C = C̄·sqrt((2-δ0m)(2n+1)(n-m)!/(n+m)!).
            let mut ratio = 1.0f64;
            for k in (n - m + 1)..=(n + m) {
                ratio /= k as f64;
            }
            let norm = ((if m == 0 { 1.0 } else { 2.0 }) * (2 * n + 1) as f64 * ratio).sqrt();
            c[idx(n, m)] = num(p[2])? * norm;
            s[idx(n, m)] = num(p[3])? * norm;
        }
        if degree > max_seen {
            return Err(bad(format!("degree {degree} requested, the table has only {max_seen}")));
        }
        Ok(Self { gm, radius_km, degree, c, s })
    }

    /// Ненормированный J2 = −C20.
    pub fn j2(&self) -> f64 {
        if self.degree >= 2 { -self.c[idx(2, 0)] } else { 0.0 }
    }

    /// V и W Каннингема до степени `nmax` в точке `r` (в осях поля).
    fn vw(&self, r: &Vector3, nmax: usize) -> (Vec<f64>, Vec<f64>) {
        let size = idx(nmax, nmax) + 1;
        let (mut v, mut w) = (vec![0.0; size], vec![0.0; size]);
        let r2 = r.norm_squared();
        let rr = self.radius_km;
        let (x, y, z) = (r.x * rr / r2, r.y * rr / r2, r.z * rr / r2);
        let q = rr * rr / r2;
        v[0] = rr / r2.sqrt();
        for m in 0..=nmax {
            if m > 0 {
                let f = (2 * m - 1) as f64;
                let (vp, wp) = (v[idx(m - 1, m - 1)], w[idx(m - 1, m - 1)]);
                v[idx(m, m)] = f * (x * vp - y * wp);
                w[idx(m, m)] = f * (x * wp + y * vp);
            }
            for n in (m + 1)..=nmax {
                let a = (2 * n - 1) as f64 / (n - m) as f64;
                let b = (n + m - 1) as f64 / (n - m) as f64;
                let (v1, w1) = (v[idx(n - 1, m)], w[idx(n - 1, m)]);
                let (v2, w2) = if n >= m + 2 { (v[idx(n - 2, m)], w[idx(n - 2, m)]) } else { (0.0, 0.0) };
                v[idx(n, m)] = a * z * v1 - b * q * v2;
                w[idx(n, m)] = a * z * w1 - b * q * w2;
            }
        }
        (v, w)
    }

    /// Степень, выше которой вклад членов на расстоянии `r_km` ниже точности f64.
    fn effective_degree(&self, r_km: f64) -> usize {
        let q = self.radius_km / r_km;
        if q >= 0.5 {
            return self.degree;
        }
        // (R/r)^n < 1e-16  =>  n > ln(1e-16) / ln(R/r).
        let n = ((1e-16f64).ln() / q.ln()).ceil() as usize;
        n.clamp(2.min(self.degree), self.degree)
    }

    /// Потенциал, км²/с² (положительный).
    pub fn potential(&self, r: &Vector3) -> f64 {
        let (v, w) = self.vw(r, self.degree);
        let mut u = 0.0;
        for i in 0..self.c.len() {
            u += self.c[i] * v[i] + self.s[i] * w[i];
        }
        u * self.gm / self.radius_km
    }

    /// Ускорение, км/с², в точке `r` (км) — всё в осях поля.
    pub fn acceleration(&self, r: &Vector3) -> Vector3 {
        if self.degree == 0 {
            return -r * (self.gm / r.norm().powi(3));
        }
        // Вдали от Луны члены степени n убывают как (R/r)^n: те, что меньше 1e-16 от основного
        // члена (ниже точности f64), не считаются. Результат тот же, но в 5-30 раз быстрее
        // на подлёте; у поверхности считается полная степень.
        let nmax = self.effective_degree(r.norm());
        let (v, w) = self.vw(r, nmax + 1);
        let (mut ax, mut ay, mut az) = (0.0, 0.0, 0.0);
        for n in 0..=nmax {
            for m in 0..=n {
                let (c, s) = (self.c[idx(n, m)], self.s[idx(n, m)]);
                if c == 0.0 && s == 0.0 {
                    continue;
                }
                let up = idx(n + 1, m + 1);
                let mid = idx(n + 1, m);
                if m == 0 {
                    ax -= c * v[up];
                    ay -= c * w[up];
                } else {
                    let dn = idx(n + 1, m - 1);
                    let f = ((n - m + 1) * (n - m + 2)) as f64;
                    ax += 0.5 * (-c * v[up] - s * w[up] + f * (c * v[dn] + s * w[dn]));
                    ay += 0.5 * (-c * w[up] + s * v[up] + f * (-c * w[dn] + s * v[dn]));
                }
                az += (n - m + 1) as f64 * (-c * v[mid] - s * w[mid]);
            }
        }
        Vector3::new(ax, ay, az) * (self.gm / (self.radius_km * self.radius_km))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_mass_matches_newton() {
        let f = GravityField::moon_grail(0).unwrap();
        let r = Vector3::new(1200.0, -900.0, 800.0);
        let a = f.acceleration(&r);
        let expect = -r * (f.gm / r.norm().powi(3));
        assert!((a - expect).norm() < 1e-15);
        assert!((f.gm - 4902.8003).abs() < 1e-3 && (f.radius_km - 1738.0).abs() < 1e-9);
    }

    #[test]
    fn j2_matches_closed_form() {
        let f = GravityField::moon_grail(2).unwrap();
        // Оставляем только C20.
        let mut only = f.clone();
        for i in 0..only.c.len() {
            if i != 0 && i != idx(2, 0) {
                only.c[i] = 0.0;
                only.s[i] = 0.0;
            }
        }
        let j2 = only.j2();
        assert!((j2 - 2.0321e-4).abs() < 2e-7, "{j2}");
        let r = Vector3::new(1500.0, 700.0, -900.0);
        let rn = r.norm();
        let k = 1.5 * j2 * only.gm * only.radius_km.powi(2) / rn.powi(5);
        let zz = 5.0 * r.z * r.z / (rn * rn);
        let expect = -r * (only.gm / rn.powi(3))
            + Vector3::new(k * r.x * (zz - 1.0), k * r.y * (zz - 1.0), k * r.z * (zz - 3.0));
        assert!((only.acceleration(&r) - expect).norm() < 1e-14, "{:?}", only.acceleration(&r) - expect);
    }

    #[test]
    fn acceleration_is_gradient_of_potential() {
        // Полное поле степени 60, включая точки над полюсом и ниже опорного радиуса.
        let f = GravityField::moon_grail(60).unwrap();
        for r in [
            Vector3::new(1100.0, 900.0, 1100.0),
            Vector3::new(0.0, 0.0, 1837.4),
            Vector3::new(1e-9, 0.0, -1752.4),
            Vector3::new(-1736.0, 30.0, 12.0),
        ] {
            let a = f.acceleration(&r);
            let d = 1e-3;
            let mut g = Vector3::zeros();
            for i in 0..3 {
                let mut e = Vector3::zeros();
                e[i] = d;
                g[i] = (f.potential(&(r + e)) - f.potential(&(r - e))) / (2.0 * d);
            }
            assert!((a - g).norm() < 2e-8 * a.norm(), "{r:?}: {:?}", a - g);
        }
        // Масконы заметны: на высоте 15 км отличие от точечной массы — доли процента.
        let r = Vector3::new(1752.4, 0.0, 0.0);
        let rel = (f.acceleration(&r) + r * (f.gm / r.norm().powi(3))).norm() / (f.gm / r.norm_squared());
        assert!(rel > 1e-4 && rel < 1e-2, "{rel}");
    }

    #[test]
    fn far_field_truncation_is_below_double_precision() {
        // Сравнение с полной суммой всех 60 степеней на разных расстояниях.
        let f = GravityField::moon_grail(60).unwrap();
        for k in [2.1, 3.0, 5.0, 10.0, 40.0] {
            let r = Vector3::new(0.6, -0.5, 0.62).normalize() * (f.radius_km * k);
            assert!(f.effective_degree(r.norm()) < 60, "{k}");
            let (v, w) = f.vw(&r, 61);
            let (mut ax, mut ay, mut az) = (0.0, 0.0, 0.0);
            for n in 0..=60 {
                for m in 0..=n {
                    let (c, s) = (f.c[idx(n, m)], f.s[idx(n, m)]);
                    let (up, mid) = (idx(n + 1, m + 1), idx(n + 1, m));
                    if m == 0 {
                        ax -= c * v[up];
                        ay -= c * w[up];
                    } else {
                        let dn = idx(n + 1, m - 1);
                        let fct = ((n - m + 1) * (n - m + 2)) as f64;
                        ax += 0.5 * (-c * v[up] - s * w[up] + fct * (c * v[dn] + s * w[dn]));
                        ay += 0.5 * (-c * w[up] + s * v[up] + fct * (-c * w[dn] + s * v[dn]));
                    }
                    az += (n - m + 1) as f64 * (-c * v[mid] - s * w[mid]);
                }
            }
            let full = Vector3::new(ax, ay, az) * (f.gm / (f.radius_km * f.radius_km));
            let a = f.acceleration(&r);
            assert!((a - full).norm() <= 1e-15 * full.norm(), "{k}: {:e}", (a - full).norm() / full.norm());
        }
    }

    #[test]
    fn rejects_bad_input() {
        assert!(GravityField::moon_grail(61).is_err());
        assert!(GravityField::parse("", 2).is_err());
        assert!(GravityField::parse("1738, 4902.8\n2,0,abc,0", 2).is_err());
    }
}
