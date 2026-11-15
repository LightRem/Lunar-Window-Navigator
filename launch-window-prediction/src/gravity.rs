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
    /// Множители рекурсии Каннингема `(2n-1)/(n-m)` и `(n+m-1)/(n-m)` (индекс `n(n+1)/2 + m`,
    /// до степени `degree + 1`): считаются один раз, а не при каждом вызове.
    rec_a: Vec<f64>,
    rec_b: Vec<f64>,
}

/// Множители рекурсии до степени `nmax` включительно.
fn recursion_tables(nmax: usize) -> (Vec<f64>, Vec<f64>) {
    let size = idx(nmax, nmax) + 1;
    let (mut a, mut b) = (vec![0.0; size], vec![0.0; size]);
    for m in 0..=nmax {
        for n in (m + 1)..=nmax {
            a[idx(n, m)] = (2 * n - 1) as f64 / (n - m) as f64;
            b[idx(n, m)] = (n + m - 1) as f64 / (n - m) as f64;
        }
    }
    (a, b)
}

thread_local! {
    /// Рабочие массивы V и W (без выделения памяти на каждый вызов).
    static SCRATCH: std::cell::RefCell<(Vec<f64>, Vec<f64>)> = const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
}

fn idx(n: usize, m: usize) -> usize {
    n * (n + 1) / 2 + m
}

impl GravityField {
    /// Точечная масса (без неоднородностей).
    pub fn point_mass(gm: f64, radius_km: f64) -> Self {
        let (rec_a, rec_b) = recursion_tables(1);
        Self { gm, radius_km, degree: 0, c: vec![1.0], s: vec![0.0], rec_a, rec_b }
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
        let (rec_a, rec_b) = recursion_tables(degree + 1);
        Ok(Self { gm, radius_km, degree, c, s, rec_a, rec_b })
    }

    /// Ненормированный J2 = −C20.
    pub fn j2(&self) -> f64 {
        if self.degree >= 2 { -self.c[idx(2, 0)] } else { 0.0 }
    }

    /// Вызывает `f` с V и W Каннингема до степени `nmax` (не больше `degree + 1`) в точке `r`
    /// (в осях поля). Массивы — рабочие массивы потока: каждый элемент до `nmax` перезаписывается.
    fn with_vw<T>(&self, r: &Vector3, nmax: usize, f: impl FnOnce(&[f64], &[f64]) -> T) -> T {
        SCRATCH.with(|cell| {
            let mut guard = cell.borrow_mut();
            let (v, w) = &mut *guard;
            let size = idx(nmax, nmax) + 1;
            if v.len() < size {
                v.resize(size, 0.0);
                w.resize(size, 0.0);
            }
            self.vw_into(r, nmax, &mut v[..size], &mut w[..size]);
            f(&v[..size], &w[..size])
        })
    }

    /// Рекурсия по строкам (степень n — подряд в памяти): каждый V(n,m), W(n,m) считается по той
    /// же формуле из тех же значений, что и в рекурсии по столбцам (Montenbruck & Gill 3.2.4).
    fn vw_into(&self, r: &Vector3, nmax: usize, v: &mut [f64], w: &mut [f64]) {
        let (ta, tb) = (&self.rec_a, &self.rec_b);
        let r2 = r.norm_squared();
        let rr = self.radius_km;
        let (x, y, z) = (r.x * rr / r2, r.y * rr / r2, r.z * rr / r2);
        let q = rr * rr / r2;
        v[0] = rr / r2.sqrt();
        w[0] = 0.0;
        for n in 1..=nmax {
            let (p0, p1) = (idx(n, 0), idx(n - 1, 0));
            let (vh, vc) = v.split_at_mut(p0);
            let (wh, wc) = w.split_at_mut(p0);
            let (vc, wc) = (&mut vc[..=n], &mut wc[..=n]);
            let (v1, w1) = (&vh[p1..p0], &wh[p1..p0]);
            let (ar, br) = (&ta[p0..p0 + n], &tb[p0..p0 + n]);
            if n >= 2 {
                let p2 = idx(n - 2, 0);
                let (v2, w2) = (&vh[p2..p1], &wh[p2..p1]);
                for m in 0..n - 1 {
                    let (a, b) = (ar[m], br[m]);
                    vc[m] = a * z * v1[m] - b * q * v2[m];
                    wc[m] = a * z * w1[m] - b * q * w2[m];
                }
            }
            let m = n - 1;
            let (a, b) = (ar[m], br[m]);
            let (v2, w2) = (0.0, 0.0);
            vc[m] = a * z * v1[m] - b * q * v2;
            wc[m] = a * z * w1[m] - b * q * w2;
            // Диагональ: из (n-1, n-1).
            let f = (2 * n - 1) as f64;
            let (vp, wp) = (v1[n - 1], w1[n - 1]);
            vc[n] = f * (x * vp - y * wp);
            wc[n] = f * (x * wp + y * vp);
        }
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
        let u = self.with_vw(r, self.degree, |v, w| {
            let mut u = 0.0;
            for i in 0..self.c.len() {
                u += self.c[i] * v[i] + self.s[i] * w[i];
            }
            u
        });
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
        let (ax, ay, az) = self.with_vw(r, nmax + 1, |v, w| self.sum_acceleration(v, w, nmax));
        Vector3::new(ax, ay, az) * (self.gm / (self.radius_km * self.radius_km))
    }

    fn sum_acceleration(&self, v: &[f64], w: &[f64], nmax: usize) -> (f64, f64, f64) {
        let (mut ax, mut ay, mut az) = (0.0, 0.0, 0.0);
        for n in 0..=nmax {
            let (pc, pv) = (idx(n, 0), idx(n + 1, 0));
            let (cr, sr) = (&self.c[pc..=pc + n], &self.s[pc..=pc + n]);
            let (vr, wr) = (&v[pv..=pv + n + 1], &w[pv..=pv + n + 1]);
            for m in 0..=n {
                let (c, s) = (cr[m], sr[m]);
                if c == 0.0 && s == 0.0 {
                    continue;
                }
                let (vu, wu, vm, wm) = (vr[m + 1], wr[m + 1], vr[m], wr[m]);
                if m == 0 {
                    ax -= c * vu;
                    ay -= c * wu;
                } else {
                    let (vd, wd) = (vr[m - 1], wr[m - 1]);
                    let f = ((n - m + 1) * (n - m + 2)) as f64;
                    ax += 0.5 * (-c * vu - s * wu + f * (c * vd + s * wd));
                    ay += 0.5 * (-c * wu + s * vu + f * (-c * wd + s * vd));
                }
                az += (n - m + 1) as f64 * (-c * vm - s * wm);
            }
        }
        (ax, ay, az)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Прежняя рекурсия по столбцам (эталон для рекурсии по строкам).
    fn vw_by_columns(f: &GravityField, r: &Vector3, nmax: usize) -> (Vec<f64>, Vec<f64>) {
        let size = idx(nmax, nmax) + 1;
        let (mut v, mut w) = (vec![0.0; size], vec![0.0; size]);
        let r2 = r.norm_squared();
        let rr = f.radius_km;
        let (x, y, z) = (r.x * rr / r2, r.y * rr / r2, r.z * rr / r2);
        let q = rr * rr / r2;
        v[0] = rr / r2.sqrt();
        for m in 0..=nmax {
            if m > 0 {
                let k = (2 * m - 1) as f64;
                let (vp, wp) = (v[idx(m - 1, m - 1)], w[idx(m - 1, m - 1)]);
                v[idx(m, m)] = k * (x * vp - y * wp);
                w[idx(m, m)] = k * (x * wp + y * vp);
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

    #[test]
    fn row_recursion_matches_column_recursion_exactly() {
        let f = GravityField::moon_grail(60).unwrap();
        for (k, r) in [1.02, 1.1, 1.6, 4.0].iter().zip([
            Vector3::new(0.3, -0.8, 0.52),
            Vector3::new(-0.1, 0.05, -0.99),
            Vector3::new(0.7, 0.7, 0.1),
            Vector3::new(0.0, 0.0, 1.0),
        ]) {
            let r = r.normalize() * (f.radius_km * k);
            let (v0, w0) = vw_by_columns(&f, &r, 61);
            let (v1, w1) = f.with_vw(&r, 61, |v, w| (v.to_vec(), w.to_vec()));
            assert_eq!(v0, v1);
            assert_eq!(w0, w1);
            // Повторный вызов (рабочие массивы уже заполнены) даёт то же самое и для меньшей степени.
            let (v2, _) = f.with_vw(&r, 10, |v, w| (v.to_vec(), w.to_vec()));
            assert_eq!(&v0[..v2.len()], &v2[..]);
        }
    }

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
            let (v, w) = f.with_vw(&r, 61, |v, w| (v.to_vec(), w.to_vec()));
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
