//! Каталог космодромов (`data/spaceports.toml`).

use std::collections::HashSet;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

/// Допуск при сравнении наклонения с коридором из `allowed_inclinations_deg`, градусы.
/// Коридоры в источниках округлены (например, «51°» у Восточного при широте 51.9°).
pub const INCLINATION_CORRIDOR_TOLERANCE_DEG: f64 = 1.0;

/// Один космодром. Ограничения на пуск необязательны: если их нет, разрешён любой
/// пуск на восток. Источники ограничений указываются в `limits_source`.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Spaceport {
    pub id: String,
    pub name: String,
    #[serde(default)]
    pub pads: Vec<String>,
    pub country: String,
    /// Геодезическая широта, градусы (+ север).
    pub lat_deg: f64,
    /// Долгота, градусы (+ восток).
    pub lon_deg: f64,
    #[serde(default)]
    pub altitude_m: Option<f64>,
    /// Допустимые азимуты пуска, градусы (0 = север, 90 = восток). Минимум может быть
    /// отрицательным (например, -10.5 = 349.5): сектор идёт от min к max по часовой стрелке.
    #[serde(default)]
    pub azimuth_min_deg: Option<f64>,
    #[serde(default)]
    pub azimuth_max_deg: Option<f64>,
    /// Допустимый диапазон наклонений орбиты, градусы.
    #[serde(default)]
    pub inclination_min_deg: Option<f64>,
    #[serde(default)]
    pub inclination_max_deg: Option<f64>,
    /// Разрешённые коридоры наклонений (там, где пуск возможен только в отдельных
    /// направлениях из-за полей падения ступеней), градусы.
    #[serde(default)]
    pub allowed_inclinations_deg: Option<Vec<f64>>,
    /// Откуда взяты ограничения.
    #[serde(default)]
    pub limits_source: Option<String>,
}

impl Spaceport {
    /// Входит ли азимут в допустимый сектор (без ограничений — любой).
    pub fn azimuth_allowed(&self, azimuth_deg: f64) -> bool {
        match (self.azimuth_min_deg, self.azimuth_max_deg) {
            (Some(lo), Some(hi)) => (azimuth_deg - lo).rem_euclid(360.0) <= hi - lo + 1e-9,
            _ => true,
        }
    }

    /// Наклонения орбиты (градусы, по возрастанию), которые разрешены космодромом и не меньше
    /// `required_min_deg`. Для диапазона и отсутствия ограничений возвращается сетка с шагом
    /// `step_deg` до `max_deg` (обычно 1 градус до 89), для коридоров — сами коридоры.
    pub fn inclination_candidates(&self, required_min_deg: f64, step_deg: f64, max_deg: f64) -> Vec<f64> {
        if let Some(list) = &self.allowed_inclinations_deg {
            let tol = INCLINATION_CORRIDOR_TOLERANCE_DEG;
            let mut out: Vec<f64> = list
                .iter()
                .filter_map(|&c| {
                    let i = required_min_deg.max(c - tol);
                    (i <= c + tol && i <= max_deg).then_some(i)
                })
                .collect();
            out.sort_by(f64::total_cmp);
            out.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
            return out;
        }
        let lo = required_min_deg.max(self.inclination_min_deg.unwrap_or(f64::NEG_INFINITY));
        let hi = self.inclination_max_deg.unwrap_or(max_deg).min(max_deg);
        let mut out = Vec::new();
        let mut i = lo;
        while i <= hi + 1e-9 {
            out.push(i);
            i += step_deg;
        }
        out
    }
}

#[derive(Debug, Deserialize)]
struct CatalogFile {
    #[serde(default)]
    spaceport: Vec<Spaceport>,
}

/// Проверенный набор космодромов.
#[derive(Debug, Clone)]
pub struct SpaceportCatalog {
    ports: Vec<Spaceport>,
}

impl SpaceportCatalog {
    /// Разбирает TOML-текст. `origin` нужен только для сообщений об ошибках.
    pub fn parse(text: &str, origin: &Path) -> Result<Self> {
        let file: CatalogFile = toml::from_str(text).map_err(|source| Error::Toml {
            path: origin.to_path_buf(),
            source,
        })?;
        let catalog = Self { ports: file.spaceport };
        catalog.validate()?;
        Ok(catalog)
    }

    /// Асинхронно читает и проверяет файл.
    pub async fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = tokio::fs::read_to_string(path)
            .await
            .map_err(|source| Error::Io { path: path.to_path_buf(), source })?;
        Self::parse(&text, path)
    }

    pub fn get(&self, id: &str) -> Option<&Spaceport> {
        self.ports.iter().find(|p| p.id == id)
    }

    pub fn iter(&self) -> impl Iterator<Item = &Spaceport> {
        self.ports.iter()
    }

    pub fn len(&self) -> usize {
        self.ports.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ports.is_empty()
    }

    fn validate(&self) -> Result<()> {
        let mut seen = HashSet::new();
        for p in &self.ports {
            if !seen.insert(p.id.as_str()) {
                return Err(Error::Config(format!("duplicate spaceport id '{}'", p.id)));
            }
            if !(-90.0..=90.0).contains(&p.lat_deg) {
                return Err(Error::Config(format!(
                    "spaceport '{}': latitude {} out of range [-90, 90]",
                    p.id, p.lat_deg
                )));
            }
            if !(-180.0..=180.0).contains(&p.lon_deg) {
                return Err(Error::Config(format!(
                    "spaceport '{}': longitude {} out of range [-180, 180]",
                    p.id, p.lon_deg
                )));
            }
            let bad = |m: &str| Err(Error::Config(format!("spaceport '{}': {m}", p.id)));
            match (p.azimuth_min_deg, p.azimuth_max_deg) {
                (None, None) => {}
                (Some(lo), Some(hi)) if lo < hi && hi - lo <= 360.0 && lo >= -360.0 && hi <= 360.0 => {}
                _ => return bad("azimuth_min_deg and azimuth_max_deg must be set together, min < max"),
            }
            let incl_ok = |i: f64| i > 0.0 && i < 180.0;
            match (p.inclination_min_deg, p.inclination_max_deg) {
                (Some(lo), Some(hi)) if !(incl_ok(lo) && incl_ok(hi) && lo <= hi) => {
                    return bad("need 0 < inclination_min_deg <= inclination_max_deg < 180");
                }
                (Some(i), None) | (None, Some(i)) if !incl_ok(i) => {
                    return bad("inclination limits must be in (0, 180)");
                }
                _ => {}
            }
            if let Some(list) = &p.allowed_inclinations_deg {
                if list.is_empty() || !list.iter().all(|&i| incl_ok(i)) {
                    return bad("allowed_inclinations_deg must be non-empty, values in (0, 180)");
                }
                if p.inclination_min_deg.is_some() || p.inclination_max_deg.is_some() {
                    return bad("use either allowed_inclinations_deg or inclination_min/max, not both");
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const BUNDLED: &str = include_str!("../data/spaceports.toml");

    #[test]
    fn bundled_file_has_eleven_valid_spaceports() {
        let c = SpaceportCatalog::parse(BUNDLED, Path::new("data/spaceports.toml")).unwrap();
        assert_eq!(c.len(), 11);
        let kourou = c.get("kourou").unwrap();
        assert!((kourou.lat_deg - 5.26258).abs() < 1e-9);
        assert!(kourou.altitude_m.is_none());
    }

    #[test]
    fn bundled_limits_are_sourced() {
        let c = SpaceportCatalog::parse(BUNDLED, Path::new("data/spaceports.toml")).unwrap();
        for p in c.iter() {
            let has_limits = p.azimuth_min_deg.is_some()
                || p.inclination_min_deg.is_some()
                || p.allowed_inclinations_deg.is_some();
            assert_eq!(has_limits, p.limits_source.is_some(), "{}: limits without a source or vice versa", p.id);
        }
    }

    fn port(text: &str) -> Spaceport {
        let base = "[[spaceport]]\nid = \"x\"\nname = \"X\"\ncountry = \"X\"\nlat_deg = 28.5\nlon_deg = 0.0\n";
        let c = SpaceportCatalog::parse(&format!("{base}{text}"), Path::new("t.toml")).unwrap();
        c.iter().next().unwrap().clone()
    }

    #[test]
    fn azimuth_sector_with_negative_minimum() {
        let p = port("azimuth_min_deg = -10.5\nazimuth_max_deg = 91.5\n");
        assert!(p.azimuth_allowed(355.0));
        assert!(p.azimuth_allowed(0.0));
        assert!(p.azimuth_allowed(91.0));
        assert!(!p.azimuth_allowed(120.0));
        assert!(!p.azimuth_allowed(340.0));
        assert!(port("").azimuth_allowed(200.0));
    }

    #[test]
    fn inclination_candidates_by_kind() {
        // Без ограничений: от требуемого минимума с шагом.
        assert_eq!(port("").inclination_candidates(29.0, 1.0, 31.0), vec![29.0, 30.0, 31.0]);
        // Диапазон.
        let r = port("inclination_min_deg = 30.0\ninclination_max_deg = 31.0\n");
        assert_eq!(r.inclination_candidates(20.0, 1.0, 89.0), vec![30.0, 31.0]);
        assert!(r.inclination_candidates(35.0, 1.0, 89.0).is_empty());
        // Коридоры с допуском 1 градус.
        let c = port("allowed_inclinations_deg = [51.6, 64.8]\n");
        assert_eq!(c.inclination_candidates(50.0, 1.0, 89.0), vec![50.6, 63.8]);
        assert_eq!(c.inclination_candidates(52.0, 1.0, 89.0), vec![52.0, 63.8]);
        assert_eq!(c.inclination_candidates(60.0, 1.0, 89.0), vec![63.8]);
    }

    #[test]
    fn rejects_inconsistent_limits() {
        let base = "[[spaceport]]\nid = \"x\"\nname = \"X\"\ncountry = \"X\"\nlat_deg = 0.0\nlon_deg = 0.0\n";
        for bad in [
            "azimuth_min_deg = 40.0\n",
            "azimuth_min_deg = 100.0\nazimuth_max_deg = 40.0\n",
            "inclination_min_deg = 50.0\ninclination_max_deg = 40.0\n",
            "allowed_inclinations_deg = []\n",
            "allowed_inclinations_deg = [51.6]\ninclination_min_deg = 40.0\n",
            "unknown_field = 1\n",
        ] {
            assert!(SpaceportCatalog::parse(&format!("{base}{bad}"), Path::new("t.toml")).is_err(), "{bad}");
        }
    }

    #[test]
    fn rejects_duplicate_ids() {
        let text = r#"
[[spaceport]]
id = "a"
name = "A"
country = "X"
lat_deg = 0.0
lon_deg = 0.0

[[spaceport]]
id = "a"
name = "B"
country = "X"
lat_deg = 1.0
lon_deg = 1.0
"#;
        assert!(SpaceportCatalog::parse(text, Path::new("t.toml")).is_err());
    }

    #[test]
    fn rejects_out_of_range_latitude() {
        let text = r#"
[[spaceport]]
id = "a"
name = "A"
country = "X"
lat_deg = 91.0
lon_deg = 0.0
"#;
        assert!(SpaceportCatalog::parse(text, Path::new("t.toml")).is_err());
    }
}
