//! Рельеф Луны: карты высот LOLA GDR (PDS3, пара `.LBL` + `.IMG`) и маска горизонта.
//!
//! Можно подключать сразу несколько карт (папка `data/terrain`): глобальную и полярные.
//! Для каждой точки берётся самая подробная карта, которая её покрывает. Поддерживаются
//! проекции `POLAR STEREOGRAPHIC` (например, `LDEM_80S_80M`, `LDEM_80N_80M`) и
//! `SIMPLE CYLINDRICAL` (например, глобальная `LDEM_64` или её тайлы).
//!
//! Все карты LOLA заданы в системе Mean Earth / Polar Axis (DE421), той же, в которой
//! [`crate::sky`] считает положение Солнца и Земли. Высоты приводятся к общей опорной
//! сфере [`REFERENCE_RADIUS_M`].

use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};

use anise::math::Vector3;
use memmap2::Mmap;
use rayon::prelude::*;
use serde::Serialize;

use crate::error::{Error, Result};
use crate::sky::{LunarSite, local_frame};

/// Общая опорная сфера для всех высот, м (та же, что у LOLA и в [`crate::sky::MOON_RADIUS_KM`]).
pub const REFERENCE_RADIUS_M: f64 = 1_737_400.0;

/// Проекция карты.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub enum Projection {
    /// Полярная стереографическая (сферическая), масштаб верен на полюсе.
    PolarStereographic {
        /// true = южный полюс, false = северный.
        south: bool,
        /// Размер пикселя на полюсе, м.
        map_scale_m: f64,
    },
    /// Равнопромежуточная цилиндрическая (широта/долгота на регулярной сетке).
    SimpleCylindrical {
        /// Пикселей на градус.
        pixels_per_deg: f64,
        center_lon_deg: f64,
        min_lat_deg: f64,
        max_lat_deg: f64,
        west_lon_deg: f64,
        east_lon_deg: f64,
    },
}

/// Параметры карты из PDS3-метки.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct GridMeta {
    pub lines: usize,
    pub samples: usize,
    /// Высота, м = DN * scaling_factor (относительно `reference_radius_m`).
    pub scaling_factor: f64,
    /// Радиус опорной сферы этой карты, м.
    pub reference_radius_m: f64,
    pub line_offset: f64,
    pub sample_offset: f64,
    pub projection: Projection,
}

impl GridMeta {
    /// Разбирает текст PDS3-метки (`KEY = value <unit>`).
    pub fn parse_label(text: &str) -> Result<Self> {
        let mut kv: HashMap<&str, &str> = HashMap::new();
        for line in text.lines() {
            if let Some((k, v)) = line.split_once('=') {
                kv.entry(k.trim()).or_insert(v.trim());
            }
        }
        let raw = |key: &str| -> Result<&str> {
            kv.get(key)
                .copied()
                .ok_or_else(|| Error::Config(format!("label: missing key {key}")))
        };
        let num = |key: &str| -> Result<f64> {
            let v = raw(key)?;
            let token = v.split_whitespace().next().unwrap_or("").trim_end_matches('.');
            token
                .parse::<f64>()
                .map_err(|_| Error::Config(format!("label: {key} = '{v}' is not a number")))
        };

        // Только карты высот: в той же папке на сервере лежат и карты числа измерений (LDAC/LDEC).
        let name = raw("NAME")?.trim_matches('"');
        if name != "HEIGHT" {
            return Err(Error::Config(format!("label: image '{name}' is not a height map")));
        }
        let sample_type = raw("SAMPLE_TYPE")?;
        let bits = num("SAMPLE_BITS")?;
        if sample_type != "LSB_INTEGER" || bits != 16.0 {
            return Err(Error::Config(format!(
                "label: unsupported sample format {sample_type}/{bits} (need LSB_INTEGER/16)"
            )));
        }

        let projection = match raw("MAP_PROJECTION_TYPE")?.trim_matches('"') {
            "POLAR STEREOGRAPHIC" => {
                let center_lat = num("CENTER_LATITUDE")?;
                let south = match center_lat {
                    -90.0 => true,
                    90.0 => false,
                    _ => {
                        return Err(Error::Config(format!(
                            "label: CENTER_LATITUDE {center_lat} is not a pole"
                        )));
                    }
                };
                if num("CENTER_LONGITUDE")? != 0.0 {
                    return Err(Error::Config("label: CENTER_LONGITUDE must be 0".into()));
                }
                Projection::PolarStereographic { south, map_scale_m: num("MAP_SCALE")? }
            }
            "SIMPLE CYLINDRICAL" => {
                if num("CENTER_LATITUDE")? != 0.0 {
                    return Err(Error::Config("label: CENTER_LATITUDE must be 0".into()));
                }
                Projection::SimpleCylindrical {
                    pixels_per_deg: num("MAP_RESOLUTION")?,
                    center_lon_deg: num("CENTER_LONGITUDE")?,
                    min_lat_deg: num("MINIMUM_LATITUDE")?,
                    max_lat_deg: num("MAXIMUM_LATITUDE")?,
                    west_lon_deg: num("WESTERNMOST_LONGITUDE")?,
                    east_lon_deg: num("EASTERNMOST_LONGITUDE")?,
                }
            }
            other => {
                return Err(Error::Config(format!("label: unsupported projection '{other}'")));
            }
        };

        Ok(Self {
            lines: num("LINES")? as usize,
            samples: num("LINE_SAMPLES")? as usize,
            scaling_factor: num("SCALING_FACTOR")?,
            reference_radius_m: num("OFFSET")?,
            line_offset: num("LINE_PROJECTION_OFFSET")?,
            sample_offset: num("SAMPLE_PROJECTION_OFFSET")?,
            projection,
        })
    }

    /// Карта покрывает все долготы (тогда по долготе индексы заворачиваются).
    fn wraps_longitude(&self) -> bool {
        matches!(self.projection, Projection::SimpleCylindrical { west_lon_deg, east_lon_deg, .. }
            if (east_lon_deg - west_lon_deg - 360.0).abs() < 1e-9)
    }

    /// (широта, долгота) в градусах -> дробные индексы (sample, line), считая от 0.
    /// Целые значения соответствуют центрам пикселей. `None`, если точка вне диапазона карты.
    pub fn lat_lon_to_pixel(&self, lat_deg: f64, lon_deg: f64) -> Option<(f64, f64)> {
        match self.projection {
            Projection::PolarStereographic { south, map_scale_m } => {
                let colat = if south { 90.0 + lat_deg } else { 90.0 - lat_deg };
                let r = 2.0 * self.reference_radius_m * (colat.to_radians() / 2.0).tan();
                let lon = lon_deg.to_radians();
                let x = r * lon.sin();
                let y = if south { r * lon.cos() } else { -r * lon.cos() };
                Some((x / map_scale_m + self.sample_offset, self.line_offset - y / map_scale_m))
            }
            Projection::SimpleCylindrical {
                pixels_per_deg,
                center_lon_deg,
                min_lat_deg,
                max_lat_deg,
                west_lon_deg,
                east_lon_deg,
            } => {
                if lat_deg < min_lat_deg || lat_deg > max_lat_deg {
                    return None;
                }
                let lon = west_lon_deg + (lon_deg - west_lon_deg).rem_euclid(360.0);
                if lon > east_lon_deg {
                    return None;
                }
                Some((
                    self.sample_offset + (lon - center_lon_deg) * pixels_per_deg,
                    self.line_offset - lat_deg * pixels_per_deg,
                ))
            }
        }
    }

    /// Обратное преобразование: индексы (sample, line) -> (широта, долгота) в градусах.
    pub fn pixel_to_lat_lon(&self, sample: f64, line: f64) -> (f64, f64) {
        match self.projection {
            Projection::PolarStereographic { south, map_scale_m } => {
                let x = (sample - self.sample_offset) * map_scale_m;
                let y = (self.line_offset - line) * map_scale_m;
                let r = x.hypot(y);
                let colat = 2.0 * (r / (2.0 * self.reference_radius_m)).atan().to_degrees();
                if south {
                    (colat - 90.0, x.atan2(y).to_degrees())
                } else {
                    (90.0 - colat, x.atan2(-y).to_degrees())
                }
            }
            Projection::SimpleCylindrical { pixels_per_deg, center_lon_deg, .. } => (
                (self.line_offset - line) / pixels_per_deg,
                center_lon_deg + (sample - self.sample_offset) / pixels_per_deg,
            ),
        }
    }

    /// Размер пикселя на местности в точке, м (для цилиндрической — по меридиану).
    pub fn resolution_m(&self, lat_deg: f64) -> f64 {
        match self.projection {
            Projection::PolarStereographic { south, map_scale_m } => {
                let colat = if south { 90.0 + lat_deg } else { 90.0 - lat_deg };
                map_scale_m * (colat.to_radians() / 2.0).cos().powi(2)
            }
            Projection::SimpleCylindrical { pixels_per_deg, .. } => {
                std::f64::consts::PI * self.reference_radius_m / (180.0 * pixels_per_deg)
            }
        }
    }
}

/// Одна карта высот, отображённая в память.
pub struct ElevationMap {
    name: String,
    meta: GridMeta,
    data: Mmap,
}

impl std::fmt::Debug for ElevationMap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ElevationMap")
            .field("name", &self.name)
            .field("meta", &self.meta)
            .finish_non_exhaustive()
    }
}

impl ElevationMap {
    /// Открывает пару `.LBL` + `.IMG` (синхронно; из async вызывать через `spawn_blocking`).
    pub fn open(label: impl AsRef<Path>, image: impl AsRef<Path>) -> Result<Self> {
        let (label, image) = (label.as_ref(), image.as_ref());
        let io = |path: &Path| {
            let path: PathBuf = path.to_path_buf();
            move |source| Error::Io { path, source }
        };
        let text = std::fs::read_to_string(label).map_err(io(label))?;
        let meta = GridMeta::parse_label(&text)?;

        let file = File::open(image).map_err(io(image))?;
        // SAFETY: файл открыт только на чтение; предполагается, что его не меняют во время работы.
        let data = unsafe { Mmap::map(&file) }.map_err(io(image))?;

        let expected = meta.lines * meta.samples * 2;
        if data.len() != expected {
            return Err(Error::Config(format!(
                "{}: size {} bytes, label expects {expected}",
                image.display(),
                data.len()
            )));
        }
        let name = label
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| label.display().to_string());
        Ok(Self { name, meta, data })
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn meta(&self) -> &GridMeta {
        &self.meta
    }

    fn dn(&self, sample: usize, line: usize) -> f64 {
        let i = (line * self.meta.samples + sample) * 2;
        i16::from_le_bytes([self.data[i], self.data[i + 1]]) as f64
    }

    /// Размер пикселя на местности в точке, м, если карта её покрывает.
    pub fn resolution_at(&self, lat_deg: f64, lon_deg: f64) -> Option<f64> {
        self.height_m(lat_deg, lon_deg).map(|_| self.meta.resolution_m(lat_deg))
    }

    /// Высота над общей опорной сферой [`REFERENCE_RADIUS_M`], м (билинейная интерполяция).
    /// `None`, если точка вне карты.
    pub fn height_m(&self, lat_deg: f64, lon_deg: f64) -> Option<f64> {
        let (s, l) = self.meta.lat_lon_to_pixel(lat_deg, lon_deg)?;
        let (n_s, n_l) = (self.meta.samples, self.meta.lines);
        let wrap = self.meta.wraps_longitude();
        let global_lat = matches!(self.meta.projection, Projection::SimpleCylindrical { .. });

        // Индексы соседних центров пикселей.
        let s_floor = s.floor();
        let l_floor = l.floor();
        let fs = s - s_floor;
        let fl = l - l_floor;
        let idx_s = |i: f64| -> Option<usize> {
            if wrap {
                Some((i as i64).rem_euclid(n_s as i64) as usize)
            } else if i >= 0.0 && (i as usize) < n_s {
                Some(i as usize)
            } else {
                None
            }
        };
        // У цилиндрической карты крайние строки продолжаются до границы диапазона широт.
        let idx_l = |i: f64| -> Option<usize> {
            if global_lat {
                Some(i.clamp(0.0, (n_l - 1) as f64) as usize)
            } else if i >= 0.0 && (i as usize) < n_l {
                Some(i as usize)
            } else {
                None
            }
        };
        let (s0, s1) = (idx_s(s_floor)?, idx_s(s_floor + 1.0)?);
        let (l0, l1) = (idx_l(l_floor)?, idx_l(l_floor + 1.0)?);

        let top = self.dn(s0, l0) * (1.0 - fs) + self.dn(s1, l0) * fs;
        let bottom = self.dn(s0, l1) * (1.0 - fs) + self.dn(s1, l1) * fs;
        let local = (top * (1.0 - fl) + bottom * fl) * self.meta.scaling_factor;
        Some(local + self.meta.reference_radius_m - REFERENCE_RADIUS_M)
    }
}

/// Набор карт: для каждой точки используется самая подробная карта, которая её покрывает.
#[derive(Debug, Default)]
pub struct TerrainSet {
    maps: Vec<ElevationMap>,
}

/// Высота точки и откуда она взята.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TerrainSample<'a> {
    pub height_m: f64,
    /// Размер пикселя использованной карты в этой точке, м.
    pub resolution_m: f64,
    pub map_name: &'a str,
}

impl TerrainSet {
    pub fn new(maps: Vec<ElevationMap>) -> Self {
        Self { maps }
    }

    /// Открывает все карты высот из папки (каждый `*.LBL` с `.IMG` рядом).
    /// Файлы, которые не являются поддерживаемыми картами высот, пропускаются с предупреждением.
    /// Синхронно; из async вызывать через `spawn_blocking`.
    pub fn open_dir(dir: impl AsRef<Path>) -> Result<Self> {
        let dir = dir.as_ref();
        let entries =
            std::fs::read_dir(dir).map_err(|source| Error::Io { path: dir.to_path_buf(), source })?;
        let mut labels: Vec<PathBuf> = entries
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().is_some_and(|e| e.eq_ignore_ascii_case("lbl")))
            .collect();
        labels.sort();

        let mut maps = Vec::new();
        for label in labels {
            let image = ["IMG", "img"]
                .iter()
                .map(|ext| label.with_extension(ext))
                .find(|p| p.exists());
            let Some(image) = image else {
                tracing::warn!(label = %label.display(), "no .IMG next to label, skipped");
                continue;
            };
            match ElevationMap::open(&label, &image) {
                Ok(map) => {
                    tracing::info!(map = map.name(), "terrain map loaded");
                    maps.push(map);
                }
                Err(e) => tracing::warn!(label = %label.display(), error = %e, "terrain map skipped"),
            }
        }
        if maps.is_empty() {
            return Err(Error::Config(format!(
                "no supported LOLA height maps (.LBL + .IMG) found in {}",
                dir.display()
            )));
        }
        Ok(Self { maps })
    }

    pub fn maps(&self) -> &[ElevationMap] {
        &self.maps
    }

    /// Высота точки по самой подробной карте, которая её покрывает.
    pub fn sample(&self, lat_deg: f64, lon_deg: f64) -> Option<TerrainSample<'_>> {
        let mut best: Option<TerrainSample<'_>> = None;
        for map in &self.maps {
            let res = map.meta.resolution_m(lat_deg);
            if best.is_some_and(|b| b.resolution_m <= res) {
                continue;
            }
            if let Some(h) = map.height_m(lat_deg, lon_deg) {
                best = Some(TerrainSample { height_m: h, resolution_m: res, map_name: map.name() });
            }
        }
        best
    }

    pub fn height_m(&self, lat_deg: f64, lon_deg: f64) -> Option<f64> {
        self.sample(lat_deg, lon_deg).map(|s| s.height_m)
    }

    /// Уклон поверхности, градусы: центральные разности по локальным осям север/восток
    /// на базе одного пикселя самой подробной карты в этой точке.
    pub fn slope_deg(&self, lat_deg: f64, lon_deg: f64) -> Option<f64> {
        self.slope_over(lat_deg, lon_deg, 0.0).map(|s| s.slope_deg)
    }

    /// Уклон на заданной базе: высоты берутся в точках на `half_baseline_m` к северу/югу
    /// и востоку/западу (база = 2 * half_baseline_m). База не бывает мельче пикселя карты:
    /// если `half_baseline_m` меньше пикселя, используется пиксель.
    pub fn slope_over(&self, lat_deg: f64, lon_deg: f64, half_baseline_m: f64) -> Option<SlopeSample<'_>> {
        let here = self.sample(lat_deg, lon_deg)?;
        let map = self.maps.iter().find(|m| m.name() == here.map_name)?;
        let delta = half_baseline_m.max(here.resolution_m);
        let up = LunarSite::new(lat_deg, lon_deg, 0.0).ok()?.up();
        let (up, north, east) = local_frame(&up);
        let at = |dir: Vector3, sign: f64| {
            let (lat, lon) = offset_lat_lon(&up, &dir, sign * delta);
            map.height_m(lat, lon)
        };
        let dzdn = (at(north, 1.0)? - at(north, -1.0)?) / (2.0 * delta);
        let dzde = (at(east, 1.0)? - at(east, -1.0)?) / (2.0 * delta);
        Some(SlopeSample {
            slope_deg: dzdn.hypot(dzde).atan().to_degrees(),
            baseline_m: 2.0 * delta,
            height_m: here.height_m,
            resolution_m: here.resolution_m,
            map_name: here.map_name,
        })
    }
}

/// Уклон в точке и откуда он взят.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlopeSample<'a> {
    pub slope_deg: f64,
    /// База, на которой измерен уклон, м.
    pub baseline_m: f64,
    pub height_m: f64,
    pub resolution_m: f64,
    pub map_name: &'a str,
}

/// Точка на расстоянии `distance_m` по поверхности от точки с единичным вектором `up`
/// в направлении `dir` (единичный, касательный). Возвращает (широта, долгота) в градусах.
pub fn offset_lat_lon(up: &Vector3, dir: &Vector3, distance_m: f64) -> (f64, f64) {
    let theta = distance_m / REFERENCE_RADIUS_M;
    let u = up * theta.cos() + dir * theta.sin();
    (u.z.clamp(-1.0, 1.0).asin().to_degrees(), u.y.atan2(u.x).to_degrees())
}

/// Параметры построения маски горизонта.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct HorizonParams {
    /// Число направлений по азимуту (равномерно по кругу).
    pub azimuths: usize,
    /// Максимальная дальность луча, км.
    pub max_range_km: f64,
    /// Высота наблюдателя (антенны/панели) над поверхностью, м.
    pub observer_height_m: f64,
}

/// Маска горизонта: угол возвышения рельефа по азимутам вокруг точки.
#[derive(Debug, Clone)]
pub struct HorizonMask {
    /// Угол возвышения рельефа, градусы; индекс i соответствует азимуту i * 360 / N.
    pub elevation_deg: Vec<f64>,
    /// Высота поверхности в точке над опорной сферой, м.
    pub site_height_m: f64,
    /// Число лучей, которые вышли за край всех карт раньше `max_range_km`
    /// (в этих направлениях дальний рельеф не учтён).
    pub truncated_rays: usize,
}

impl HorizonMask {
    /// Строит маску для точки (`lat`, `lon`). Ошибка, если точка не покрыта ни одной картой.
    pub fn build(terrain: &TerrainSet, lat_deg: f64, lon_deg: f64, params: HorizonParams) -> Result<Self> {
        if params.azimuths < 4 || !(params.max_range_km > 0.0) || !(params.observer_height_m >= 0.0)
        {
            return Err(Error::Config("invalid horizon parameters".into()));
        }
        let here = terrain.sample(lat_deg, lon_deg).ok_or_else(|| {
            Error::Config(format!("site ({lat_deg}, {lon_deg}) is not covered by any terrain map"))
        })?;
        let site_height_m = here.height_m;

        let r_ref = REFERENCE_RADIUS_M;
        let site = LunarSite::new(lat_deg, lon_deg, 0.0)?;
        let (up, north, east) = local_frame(&site.up());
        let observer = up * (r_ref + site_height_m + params.observer_height_m);
        let max_range_m = params.max_range_km * 1000.0;

        let rays: Vec<(f64, bool)> = (0..params.azimuths)
            .into_par_iter()
            .map(|i| {
                let az = (i as f64) * std::f64::consts::TAU / (params.azimuths as f64);
                let dir = north * az.cos() + east * az.sin();
                let mut best = f64::NEG_INFINITY;
                let mut truncated = false;
                let mut d = here.resolution_m;
                while d <= max_range_m {
                    let theta = d / r_ref;
                    let u = up * theta.cos() + dir * theta.sin();
                    let lat = u.z.clamp(-1.0, 1.0).asin().to_degrees();
                    let lon = u.y.atan2(u.x).to_degrees();
                    let Some(s) = terrain.sample(lat, lon) else {
                        truncated = true;
                        break;
                    };
                    let los = u * (r_ref + s.height_m) - observer;
                    let el = (los.dot(&up) / los.norm()).clamp(-1.0, 1.0).asin();
                    best = best.max(el);
                    // Шаг: не мельче пикселя карты в этом месте и растёт с дальностью
                    // (угловой размер пикселя вдали меньше).
                    d += s.resolution_m.max(d * 0.005);
                }
                (best.to_degrees(), truncated)
            })
            .collect();

        Ok(Self {
            elevation_deg: rays.iter().map(|r| r.0).collect(),
            site_height_m,
            truncated_rays: rays.iter().filter(|r| r.1).count(),
        })
    }

    /// Угол возвышения рельефа в направлении `azimuth_deg` (линейная интерполяция).
    pub fn at(&self, azimuth_deg: f64) -> f64 {
        let n = self.elevation_deg.len();
        let pos = azimuth_deg.rem_euclid(360.0) / 360.0 * (n as f64);
        let i = (pos.floor() as usize) % n;
        let f = pos - pos.floor();
        self.elevation_deg[i] * (1.0 - f) + self.elevation_deg[(i + 1) % n] * f
    }
}

#[cfg(test)]
pub(crate) mod test_util {
    use super::*;
    use std::io::Write;

    pub fn polar_label(n: usize, south: bool) -> String {
        format!(
            "NAME = HEIGHT\nLINES = {n}\nLINE_SAMPLES = {n}\nSAMPLE_TYPE = LSB_INTEGER\nSAMPLE_BITS = 16\n\
             SCALING_FACTOR = 0.5\nOFFSET = 1737400.\nMAP_PROJECTION_TYPE = \"POLAR STEREOGRAPHIC\"\n\
             CENTER_LATITUDE = {lat} <deg>\nCENTER_LONGITUDE = 0 <deg>\nMAP_SCALE = 80 <m/pix>\n\
             LINE_PROJECTION_OFFSET = {o} <pix>\nSAMPLE_PROJECTION_OFFSET = {o} <pix>\n",
            lat = if south { -90 } else { 90 },
            o = (n as f64) / 2.0 - 0.5
        )
    }

    pub fn cyl_label(ppd: usize) -> String {
        let (lines, samples) = (180 * ppd, 360 * ppd);
        format!(
            "NAME = HEIGHT\nLINES = {lines}\nLINE_SAMPLES = {samples}\nSAMPLE_TYPE = LSB_INTEGER\n\
             SAMPLE_BITS = 16\nSCALING_FACTOR = 0.5\nOFFSET = 1737400.\n\
             MAP_PROJECTION_TYPE = \"SIMPLE CYLINDRICAL\"\nMAP_RESOLUTION = {ppd} <pix/deg>\n\
             CENTER_LATITUDE = 0 <deg>\nCENTER_LONGITUDE = 180 <deg>\nMAXIMUM_LATITUDE = 90 <deg>\n\
             MINIMUM_LATITUDE = -90 <deg>\nWESTERNMOST_LONGITUDE = 0 <deg>\n\
             EASTERNMOST_LONGITUDE = 360 <deg>\nLINE_PROJECTION_OFFSET = {lo} <pix>\n\
             SAMPLE_PROJECTION_OFFSET = {so} <pix>\n",
            lo = (lines as f64) / 2.0 - 0.5,
            so = (samples as f64) / 2.0 - 0.5,
        )
    }

    pub fn write_map(name: &str, label: &str, lines: usize, samples: usize, dn: impl Fn(usize, usize) -> i16) -> ElevationMap {
        let dir = std::env::temp_dir().join(format!("lwp-terrain-{name}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (lbl, img) = (dir.join(format!("{name}.lbl")), dir.join(format!("{name}.img")));
        std::fs::write(&lbl, label).unwrap();
        let mut buf = Vec::with_capacity(lines * samples * 2);
        for line in 0..lines {
            for sample in 0..samples {
                buf.write_all(&dn(sample, line).to_le_bytes()).unwrap();
            }
        }
        std::fs::write(&img, buf).unwrap();
        ElevationMap::open(&lbl, &img).unwrap()
    }

    pub fn polar(name: &str, n: usize, dn: impl Fn(usize, usize) -> i16) -> ElevationMap {
        write_map(name, &polar_label(n, true), n, n, dn)
    }

    /// Южная полярная карта n x n с заданным размером пикселя, м.
    pub fn polar_scaled(name: &str, n: usize, scale_m: f64, dn: impl Fn(usize, usize) -> i16) -> ElevationMap {
        let label = polar_label(n, true).replace("MAP_SCALE = 80 <m/pix>", &format!("MAP_SCALE = {scale_m} <m/pix>"));
        write_map(name, &label, n, n, dn)
    }
}

#[cfg(test)]
mod tests {
    use super::test_util::*;
    use super::*;
    use approx::assert_abs_diff_eq;

    const REAL_DIR: &str = "data/terrain";

    fn params() -> HorizonParams {
        HorizonParams { azimuths: 360, max_range_km: 10.0, observer_height_m: 2.0 }
    }

    #[test]
    fn polar_projection_axes_and_roundtrip() {
        let m = GridMeta::parse_label(&polar_label(7600, true)).unwrap();
        let (s, l) = m.lat_lon_to_pixel(-90.0, 0.0).unwrap();
        assert_abs_diff_eq!(s, 3799.5, epsilon = 1e-9);
        assert_abs_diff_eq!(l, 3799.5, epsilon = 1e-9);
        // Долгота 0 смотрит «вверх» карты (line убывает), долгота 90 — вправо (sample растёт).
        let (s0, l0) = m.lat_lon_to_pixel(-89.0, 0.0).unwrap();
        assert!(l0 < 3799.5 && (s0 - 3799.5).abs() < 1e-9);
        let (s90, l90) = m.lat_lon_to_pixel(-89.0, 90.0).unwrap();
        assert!(s90 > 3799.5 && (l90 - 3799.5).abs() < 1e-9);
        // 1 градус от полюса ~ 30.32 км ~ 379 пикселей по 80 м.
        assert_abs_diff_eq!(3799.5 - l0, 379.05, epsilon = 0.05);
        let (lat, lon) = m.pixel_to_lat_lon(s90 + 12.3, l90 - 45.6);
        let (s2, l2) = m.lat_lon_to_pixel(lat, lon).unwrap();
        assert_abs_diff_eq!(s2, s90 + 12.3, epsilon = 1e-6);
        assert_abs_diff_eq!(l2, l90 - 45.6, epsilon = 1e-6);
    }

    #[test]
    fn cylindrical_projection_pixel_centres_and_wrap() {
        let m = GridMeta::parse_label(&cyl_label(4)).unwrap();
        // Первый пиксель: центр на широте 90 - 0.125, долготе 0.125.
        let (s, l) = m.lat_lon_to_pixel(89.875, 0.125).unwrap();
        assert_abs_diff_eq!(s, 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(l, 0.0, epsilon = 1e-9);
        // Отрицательная долгота заворачивается: -0.125 == 359.875 -> последний столбец.
        let (s, _) = m.lat_lon_to_pixel(0.0, -0.125).unwrap();
        assert_abs_diff_eq!(s, 1439.0, epsilon = 1e-9);
        let (lat, lon) = m.pixel_to_lat_lon(100.0, 200.0);
        let (s2, l2) = m.lat_lon_to_pixel(lat, lon).unwrap();
        assert_abs_diff_eq!(s2, 100.0, epsilon = 1e-9);
        assert_abs_diff_eq!(l2, 200.0, epsilon = 1e-9);
    }

    #[test]
    fn label_rejects_unsupported_files() {
        let text = polar_label(10, true).replace("POLAR STEREOGRAPHIC", "ORTHOGRAPHIC");
        assert!(GridMeta::parse_label(&text).is_err());
        let counts = polar_label(10, true).replace("NAME = HEIGHT", "NAME = COUNTS");
        assert!(GridMeta::parse_label(&counts).is_err());
    }

    #[test]
    fn image_size_must_match_label() {
        let dir = std::env::temp_dir().join(format!("lwp-terrain-size-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("t.lbl"), polar_label(10, true)).unwrap();
        std::fs::write(dir.join("t.img"), [0u8; 50]).unwrap();
        assert!(ElevationMap::open(dir.join("t.lbl"), dir.join("t.img")).is_err());
    }

    #[test]
    fn cylindrical_map_covers_whole_moon_including_seam_and_poles() {
        // Высота = номер столбца (по долготе), чтобы проверить заворот у долготы 0/360.
        let map = write_map("cyl", &cyl_label(2), 360, 720, |s, _| s as i16);
        let set = TerrainSet::new(vec![map]);
        for (lat, lon) in [(0.0, 0.0), (90.0, 10.0), (-90.0, 200.0), (45.0, -170.0)] {
            assert!(set.height_m(lat, lon).is_some(), "({lat}, {lon})");
        }
        // Ровно на шве между последним (719) и первым (0) столбцом — среднее (719 + 0) / 2 * 0.5.
        assert_abs_diff_eq!(set.height_m(10.0, 0.0).unwrap(), 179.75, epsilon = 1e-9);
    }

    #[test]
    fn finest_map_wins_where_maps_overlap() {
        let global = write_map("glob", &cyl_label(2), 360, 720, |_, _| 1000);
        let pole = polar("pole", 100, |_, _| 200);
        let set = TerrainSet::new(vec![global, pole]);
        let at_pole = set.sample(-90.0, 0.0).unwrap();
        assert_eq!(at_pole.map_name, "pole");
        assert_abs_diff_eq!(at_pole.height_m, 100.0, epsilon = 1e-9);
        let elsewhere = set.sample(0.0, 0.0).unwrap();
        assert_eq!(elsewhere.map_name, "glob");
        assert_abs_diff_eq!(elsewhere.height_m, 500.0, epsilon = 1e-9);
    }

    #[test]
    fn heights_are_converted_to_common_reference() {
        let label = polar_label(100, true).replace("OFFSET = 1737400.", "OFFSET = 1737000.");
        let map = write_map("ref", &label, 100, 100, |_, _| 0);
        assert_abs_diff_eq!(map.height_m(-90.0, 0.0).unwrap(), -400.0, epsilon = 1e-9);
    }

    #[test]
    fn slope_of_flat_and_inclined_plane() {
        let flat = TerrainSet::new(vec![polar("flat", 100, |_, _| 200)]);
        assert_abs_diff_eq!(flat.slope_deg(-90.0, 0.0).unwrap(), 0.0, epsilon = 1e-9);
        // 40 м на пиксель 80 м вдоль sample: уклон atan(0.5) = 26.565 градуса.
        let ramp = TerrainSet::new(vec![polar("ramp", 100, |s, _| (s as i16) * 80)]);
        assert_abs_diff_eq!(ramp.slope_deg(-90.0, 0.0).unwrap(), 26.565, epsilon = 0.01);
    }

    #[test]
    fn slope_baseline_is_never_finer_than_a_pixel() {
        let ramp = TerrainSet::new(vec![polar("ramp2", 100, |s, _| (s as i16) * 80)]);
        let fine = ramp.slope_over(-90.0, 0.0, 10.0).unwrap();
        assert_abs_diff_eq!(fine.baseline_m, 160.0, epsilon = 1e-9);
        let wide = ramp.slope_over(-90.0, 0.0, 400.0).unwrap();
        assert_abs_diff_eq!(wide.baseline_m, 800.0, epsilon = 1e-9);
        // На плоскости уклон от базы не зависит.
        assert_abs_diff_eq!(fine.slope_deg, wide.slope_deg, epsilon = 1e-6);
    }

    #[test]
    fn flat_terrain_gives_slightly_negative_horizon() {
        let set = TerrainSet::new(vec![polar("flat2", 400, |_, _| 0)]);
        let mask = HorizonMask::build(&set, -90.0, 0.0, params()).unwrap();
        for el in &mask.elevation_deg {
            assert!(*el < 0.0 && *el > -1.5, "elevation {el}");
        }
    }

    #[test]
    fn wall_blocks_only_its_direction() {
        // Стена высотой 1000 м в верхней части карты (направление долготы 0), в 4 км от полюса.
        let set = TerrainSet::new(vec![polar("wall", 400, |_, line| if line < 150 { 2000 } else { 0 })]);
        let mask = HorizonMask::build(&set, -90.0, 0.0, params()).unwrap();
        // На полюсе азимут 0 = направление меридиана 0. atan(998 / ~4000) ~ 14 градусов.
        let towards = mask.at(0.0);
        assert!(towards > 12.0 && towards < 16.0, "towards {towards}");
        assert!(mask.at(180.0) < 0.0);
    }

    #[test]
    fn site_outside_all_maps_is_an_error() {
        let set = TerrainSet::new(vec![polar("out", 100, |_, _| 0)]);
        assert!(HorizonMask::build(&set, -80.0, 0.0, params()).is_err());
    }

    // ---- Проверки на настоящих картах из data/terrain: `cargo test -- --ignored` ----

    fn real() -> TerrainSet {
        TerrainSet::open_dir(REAL_DIR).unwrap()
    }

    #[test]
    #[ignore = "needs LOLA maps in data/terrain"]
    fn real_known_heights() {
        let t = real();
        // Дно кратера Шеклтон — на километры ниже опорной сферы, массив Малаперт — выше.
        assert!(t.height_m(-89.67, 129.78).unwrap() < -2500.0);
        assert!(t.height_m(-85.99, 2.9).unwrap() > 4500.0);
        // Самая высокая точка Луны (Smith et al., 2010): 5.4125 N, 201.3665 E, около 10.8 км.
        let summit = t.height_m(5.4125, 201.3665).unwrap();
        assert!(summit > 9500.0 && summit < 11500.0, "summit {summit}");
        // Море Спокойствия (место посадки Apollo 11) — около -1.9 км.
        let a11 = t.height_m(0.674, 23.473).unwrap();
        assert!(a11 > -2500.0 && a11 < -1500.0, "Apollo 11 {a11}");
    }

    #[test]
    #[ignore = "needs LOLA maps in data/terrain"]
    fn real_detailed_polar_map_matches_coarser_one() {
        // Карта 10 м и карта 80 м должны совпадать по положению (нет сдвига или поворота).
        let t = real();
        let fine = t.maps().iter().find(|m| m.name() == "LDEM_875S_10M");
        let coarse = t.maps().iter().find(|m| m.name() == "LDEM_80S_80M");
        let (Some(fine), Some(coarse)) = (fine, coarse) else { return };
        let mut diffs: Vec<f64> = (0..180)
            .filter_map(|i| {
                let lon = (i as f64) * 2.0 + 0.7;
                let lat = -88.0 - (i % 7) as f64 * 0.25;
                Some((fine.height_m(lat, lon)? - coarse.height_m(lat, lon)?).abs())
            })
            .collect();
        assert_eq!(diffs.len(), 180);
        diffs.sort_by(f64::total_cmp);
        let median = diffs[diffs.len() / 2];
        eprintln!("10 m vs 80 m: median |dh| = {median:.2} m, max = {:.1} m", diffs[diffs.len() - 1]);
        assert!(median < 15.0, "median |difference| {median} m");
    }

    #[test]
    #[ignore = "needs LOLA maps in data/terrain"]
    fn real_maps_agree_in_overlap() {
        // В зоне перекрытия (80–85 градусов) полярные и глобальная карта должны давать близкие высоты.
        // Это проверяет, что проекции и ориентация долготы у обеих карт поняты одинаково.
        let t = real();
        let global = t.maps().iter().find(|m| m.name() == "LDEM_64").expect("LDEM_64");
        for polar_name in ["LDEM_80S_80M", "LDEM_80N_80M"] {
            let Some(polar) = t.maps().iter().find(|m| m.name() == polar_name) else { continue };
            let lat = if polar_name.contains('S') { -83.0 } else { 83.0 };
            let mut diffs: Vec<f64> = (0..72)
                .filter_map(|i| {
                    let lon = (i as f64) * 5.0 + 1.3;
                    Some((polar.height_m(lat, lon)? - global.height_m(lat, lon)?).abs())
                })
                .collect();
            assert_eq!(diffs.len(), 72);
            diffs.sort_by(f64::total_cmp);
            let median = diffs[diffs.len() / 2];
            assert!(median < 60.0, "{polar_name}: median |difference| {median} m");
        }
    }
}
