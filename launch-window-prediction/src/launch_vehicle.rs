//! Ракеты-носители: каталог и гибридный ввод (как у лунного аппарата, см. [`crate::spacecraft`]).
//!
//! - **Шаблоны** — встроенный каталог `data/launch_vehicles/presets.toml` (источник у каждой ракеты).
//! - **Переопределение** — любой параметр можно изменить (`--rk поле=значение`); профиль тогда
//!   получает статус `Custom (на базе <ракета>)`.
//! - **Свои профили** — `data/launch_vehicles/custom/<id>.toml` (`based_on` + изменённые поля).
//!
//! # Как ракета влияет на перелёт
//!
//! - Масса лунного аппарата не должна превышать грузоподъёмность ракеты к Луне (`payload_tli_kg`).
//! - Предельный Δv разгона к Луне для конкретной массы: паспортная грузоподъёмность к Луне
//!   соответствует разгону [`REFERENCE_TLI_DELTA_V_M_S`]; если груз легче, разгонная ступень
//!   может дать больше — прибавка считается по формуле Циолковского для этой ступени
//!   (с полной заправкой ступени, приближение).
//! - Старт: шаблон — только со своих космодромов; пользовательская ракета (свой профиль или
//!   шаблон с переопределениями) — с любого космодрома каталога.
//! - Время выведения и ожидания — если заданы.

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};
use crate::spacecraft::{G0, params_from_table, valid_id};

/// Δv разгона к Луне, которому соответствует паспортная грузоподъёмность к Луне, м/с.
/// Допущение: типичный разгон с низкой опорной орбиты к Луне (~3.1-3.2 км/с).
pub const REFERENCE_TLI_DELTA_V_M_S: f64 = 3150.0;

const PRESETS_TOML: &str = include_str!("../data/launch_vehicles/presets.toml");

/// Частичный набор параметров ракеты (шаблон, профиль или переопределения).
#[derive(Debug, Clone, Default, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct LaunchVehicleParams {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub spaceports: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub payload_tli_kg: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub payload_leo_kg: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub payload_gto_kg: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub upper_stage_name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub upper_stage_dry_mass_kg: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub upper_stage_propellant_kg: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub upper_stage_thrust_n: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub upper_stage_isp_s: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ascent_duration_s: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_coast_s: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parking_altitude_km: Option<f64>,
    /// Ступени для моделирования подъёма (заменяется целиком).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ascent: Option<crate::ascent::AscentVehicle>,
}

const NUMERIC_FIELDS: [&str; 10] = [
    "payload_tli_kg",
    "payload_leo_kg",
    "payload_gto_kg",
    "upper_stage_dry_mass_kg",
    "upper_stage_propellant_kg",
    "upper_stage_thrust_n",
    "upper_stage_isp_s",
    "ascent_duration_s",
    "max_coast_s",
    "parking_altitude_km",
];

impl LaunchVehicleParams {
    fn field_mut(&mut self, field: &str) -> Option<&mut Option<f64>> {
        Some(match field {
            "payload_tli_kg" => &mut self.payload_tli_kg,
            "payload_leo_kg" => &mut self.payload_leo_kg,
            "payload_gto_kg" => &mut self.payload_gto_kg,
            "upper_stage_dry_mass_kg" => &mut self.upper_stage_dry_mass_kg,
            "upper_stage_propellant_kg" => &mut self.upper_stage_propellant_kg,
            "upper_stage_thrust_n" => &mut self.upper_stage_thrust_n,
            "upper_stage_isp_s" => &mut self.upper_stage_isp_s,
            "ascent_duration_s" => &mut self.ascent_duration_s,
            "max_coast_s" => &mut self.max_coast_s,
            "parking_altitude_km" => &mut self.parking_altitude_km,
            _ => return None,
        })
    }

    fn field(&self, field: &str) -> Option<f64> {
        self.clone().field_mut(field).and_then(|v| *v)
    }

    /// Применяет `поле=значение`. Для `spaceports` — список через запятую.
    pub fn set_from_str(&mut self, assignment: &str) -> Result<()> {
        let (key, value) = assignment
            .split_once('=')
            .ok_or_else(|| Error::Config(format!("rocket override '{assignment}' must be field=value")))?;
        let (key, value) = (key.trim(), value.trim());
        match key {
            "name" => self.name = Some(value.to_string()),
            "upper_stage_name" => self.upper_stage_name = Some(value.to_string()),
            "spaceports" => {
                return Err(Error::Config(
                    "custom rockets launch from any spaceport; 'spaceports' cannot be overridden".into(),
                ));
            }
            _ => {
                let slot = self.field_mut(key).ok_or_else(|| {
                    Error::Config(format!(
                        "unknown rocket field '{key}'; known: name, upper_stage_name, {}",
                        NUMERIC_FIELDS.join(", ")
                    ))
                })?;
                *slot = Some(value.parse().map_err(|_| {
                    Error::Config(format!("rocket field '{key}': '{value}' is not a number"))
                })?);
            }
        }
        Ok(())
    }

    pub fn set_fields(&self) -> Vec<&'static str> {
        let mut out = Vec::new();
        if self.name.is_some() {
            out.push("name");
        }
        if self.spaceports.is_some() {
            out.push("spaceports");
        }
        if self.upper_stage_name.is_some() {
            out.push("upper_stage_name");
        }
        out.extend(NUMERIC_FIELDS.iter().copied().filter(|f| self.field(f).is_some()));
        if self.ascent.is_some() {
            out.push("ascent");
        }
        out
    }

    pub fn is_empty(&self) -> bool {
        self.set_fields().is_empty()
    }

    pub fn overlay(&self, over: &LaunchVehicleParams) -> LaunchVehicleParams {
        let mut out = self.clone();
        if over.name.is_some() {
            out.name = over.name.clone();
        }
        if over.spaceports.is_some() {
            out.spaceports = over.spaceports.clone();
        }
        if over.upper_stage_name.is_some() {
            out.upper_stage_name = over.upper_stage_name.clone();
        }
        if over.ascent.is_some() {
            out.ascent = over.ascent.clone();
        }
        for f in NUMERIC_FIELDS {
            if let Some(v) = over.field(f) {
                *out.field_mut(f).expect("known field") = Some(v);
            }
        }
        out
    }
}

/// Полная конфигурация ракеты.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LaunchVehicleConfig {
    pub name: String,
    /// Космодромы шаблона (для пользовательской ракеты — справочно, см. `any_spaceport`).
    pub spaceports: Vec<String>,
    /// Пользовательская ракета: старт с любого космодрома каталога.
    pub any_spaceport: bool,
    pub payload_tli_kg: f64,
    pub payload_leo_kg: Option<f64>,
    pub payload_gto_kg: Option<f64>,
    pub upper_stage_name: Option<String>,
    pub upper_stage_dry_mass_kg: f64,
    pub upper_stage_propellant_kg: f64,
    pub upper_stage_thrust_n: f64,
    pub upper_stage_isp_s: f64,
    pub ascent_duration_s: Option<f64>,
    pub max_coast_s: Option<f64>,
    pub parking_altitude_km: Option<f64>,
    /// Ступени для моделирования подъёма (если известны).
    pub ascent: Option<crate::ascent::AscentVehicle>,
}

impl LaunchVehicleConfig {
    /// Разгонная ступень в виде для моделирования подъёма.
    pub fn upper_stage(&self) -> crate::ascent::UpperStage {
        crate::ascent::UpperStage {
            dry_mass_kg: self.upper_stage_dry_mass_kg,
            propellant_kg: self.upper_stage_propellant_kg,
            thrust_n: self.upper_stage_thrust_n,
            isp_s: self.upper_stage_isp_s,
        }
    }
}

/// Разгон к Луне разгонной ступенью.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct TliBurn {
    pub duration_s: f64,
    pub propellant_used_kg: f64,
    /// По массам ступени топлива на такой разгон с круговой опорной орбиты не хватает.
    /// Так бывает, когда ступень начинает разгон с более высокой (вытянутой) орбиты, на которую
    /// её вывели нижние ступени (например, ICPS у SLS); паспортная грузоподъёмность это учитывает.
    pub exceeds_stage_propellant: bool,
}

/// Возможности ракеты для конкретной массы груза.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LaunchCapability {
    pub payload_kg: f64,
    pub payload_tli_kg: f64,
    /// Груз не тяжелее паспортной грузоподъёмности к Луне.
    pub within_rated_payload: bool,
    pub reference_tli_delta_v_m_s: f64,
    /// Предельный Δv разгона к Луне для этого груза, м/с.
    pub max_tli_delta_v_m_s: f64,
    /// Скорость истечения разгонной ступени, м/с.
    pub upper_stage_exhaust_velocity_m_s: f64,
}

impl LaunchVehicleConfig {
    fn from_params(p: &LaunchVehicleParams, fallback_name: &str, custom: bool) -> Result<Self> {
        let req = |f: &str| p.field(f);
        let missing: Vec<&str> = [
            "payload_tli_kg",
            "upper_stage_dry_mass_kg",
            "upper_stage_propellant_kg",
            "upper_stage_thrust_n",
            "upper_stage_isp_s",
        ]
        .into_iter()
        .filter(|f| req(f).is_none())
        .collect();
        let mut missing: Vec<String> = missing.into_iter().map(String::from).collect();
        if !custom && p.spaceports.as_ref().is_none_or(|s| s.is_empty()) {
            missing.push("spaceports".into());
        }
        if !missing.is_empty() {
            return Err(Error::Config(format!(
                "rocket '{fallback_name}': missing parameters: {}",
                missing.join(", ")
            )));
        }
        Ok(Self {
            name: p.name.clone().unwrap_or_else(|| fallback_name.to_string()),
            spaceports: p.spaceports.clone().unwrap_or_default(),
            any_spaceport: custom,
            payload_tli_kg: p.payload_tli_kg.unwrap_or_default(),
            payload_leo_kg: p.payload_leo_kg,
            payload_gto_kg: p.payload_gto_kg,
            upper_stage_name: p.upper_stage_name.clone(),
            upper_stage_dry_mass_kg: p.upper_stage_dry_mass_kg.unwrap_or_default(),
            upper_stage_propellant_kg: p.upper_stage_propellant_kg.unwrap_or_default(),
            upper_stage_thrust_n: p.upper_stage_thrust_n.unwrap_or_default(),
            upper_stage_isp_s: p.upper_stage_isp_s.unwrap_or_default(),
            ascent_duration_s: p.ascent_duration_s,
            max_coast_s: p.max_coast_s,
            parking_altitude_km: p.parking_altitude_km,
            ascent: p.ascent.clone(),
        })
    }

    /// Критические ошибки (пустой список — всё корректно).
    pub fn validate(&self) -> Vec<String> {
        let mut errs = Vec::new();
        if let Some(a) = &self.ascent {
            errs.extend(a.validate());
        }
        for (f, v) in [
            ("payload_tli_kg", self.payload_tli_kg),
            ("upper_stage_dry_mass_kg", self.upper_stage_dry_mass_kg),
            ("upper_stage_propellant_kg", self.upper_stage_propellant_kg),
            ("upper_stage_thrust_n", self.upper_stage_thrust_n),
            ("upper_stage_isp_s", self.upper_stage_isp_s),
        ] {
            if !(v > 0.0 && v.is_finite()) {
                errs.push(format!("{f} must be > 0 (got {v})"));
            }
        }
        for (f, v) in [
            ("payload_leo_kg", self.payload_leo_kg),
            ("payload_gto_kg", self.payload_gto_kg),
            ("ascent_duration_s", self.ascent_duration_s),
            ("max_coast_s", self.max_coast_s),
        ] {
            if let Some(v) = v
                && !(v > 0.0 && v.is_finite())
            {
                errs.push(format!("{f} must be > 0 (got {v})"));
            }
        }
        if let Some(h) = self.parking_altitude_km
            && !(120.0..=2000.0).contains(&h)
        {
            errs.push(format!("parking_altitude_km must be in [120, 2000] (got {h})"));
        }
        if !self.any_spaceport && self.spaceports.is_empty() {
            errs.push("spaceports must not be empty".into());
        }
        errs
    }

    /// Возможности для груза массой `payload_kg` (обычно стартовая масса лунного аппарата).
    pub fn capability(&self, payload_kg: f64) -> LaunchCapability {
        let c = self.upper_stage_isp_s * G0;
        let (dry, prop) = (self.upper_stage_dry_mass_kg, self.upper_stage_propellant_kg);
        let ratio = |p: f64| ((dry + prop + p) / (dry + p)).ln();
        let gain = c * (ratio(payload_kg) - ratio(self.payload_tli_kg));
        LaunchCapability {
            payload_kg,
            payload_tli_kg: self.payload_tli_kg,
            within_rated_payload: payload_kg <= self.payload_tli_kg,
            reference_tli_delta_v_m_s: REFERENCE_TLI_DELTA_V_M_S,
            max_tli_delta_v_m_s: REFERENCE_TLI_DELTA_V_M_S + gain,
            upper_stage_exhaust_velocity_m_s: c,
        }
    }

    /// Разгон к Луне для груза и Δv (ступень с полной заправкой в начале разгона, приближение).
    pub fn tli_burn(&self, payload_kg: f64, delta_v_m_s: f64) -> TliBurn {
        let c = self.upper_stage_isp_s * G0;
        let m0 = self.upper_stage_dry_mass_kg + self.upper_stage_propellant_kg + payload_kg;
        let used = m0 * (1.0 - (-delta_v_m_s / c).exp());
        TliBurn {
            duration_s: used * c / self.upper_stage_thrust_n,
            propellant_used_kg: used,
            exceeds_stage_propellant: used > self.upper_stage_propellant_kg,
        }
    }
}

// ---------------------------------------------------------------------------
// Каталог
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LaunchVehicleEntry {
    pub id: String,
    pub is_preset: bool,
    pub based_on: Option<String>,
    pub source: Option<String>,
    pub params: LaunchVehicleParams,
    pub path: Option<PathBuf>,
}

/// Готовая к расчёту ракета.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ResolvedLaunchVehicle {
    pub id: String,
    /// `Preset: …`, `Custom (на базе …)` или `Custom`.
    pub status_label: String,
    pub config: LaunchVehicleConfig,
    pub overridden: Vec<&'static str>,
    pub sources: Vec<String>,
    pub errors: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct LaunchVehicleLibrary {
    entries: Vec<LaunchVehicleEntry>,
    custom_dir: PathBuf,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct PresetFile {
    #[serde(default)]
    preset: Vec<toml::Table>,
}

fn take_string(t: &mut toml::Table, key: &str) -> Result<Option<String>> {
    match t.remove(key) {
        None => Ok(None),
        Some(toml::Value::String(s)) => Ok(Some(s)),
        Some(_) => Err(Error::Config(format!("'{key}' must be a string"))),
    }
}

impl LaunchVehicleLibrary {
    pub fn presets_only() -> Result<Self> {
        let file: PresetFile = toml::from_str(PRESETS_TOML)
            .map_err(|source| Error::Toml { path: PathBuf::from("data/launch_vehicles/presets.toml"), source })?;
        let mut entries = Vec::new();
        for mut row in file.preset {
            let id = take_string(&mut row, "id")?.ok_or_else(|| Error::Config("rocket preset without id".into()))?;
            let source = take_string(&mut row, "source")?;
            let params = params_from_table(row, &format!("rocket preset '{id}'"))?;
            entries.push(LaunchVehicleEntry { id, is_preset: true, based_on: None, source, params, path: None });
        }
        let lib = Self { entries, custom_dir: PathBuf::new() };
        lib.check_ids()?;
        Ok(lib)
    }

    /// Каталог + профили из `custom_dir` (папки может не быть). Синхронно.
    pub fn load(custom_dir: impl AsRef<Path>) -> Result<Self> {
        let custom_dir = custom_dir.as_ref().to_path_buf();
        let mut lib = Self::presets_only()?;
        lib.custom_dir = custom_dir.clone();
        if custom_dir.is_dir() {
            let mut files: Vec<PathBuf> = std::fs::read_dir(&custom_dir)
                .map_err(|source| Error::Io { path: custom_dir.clone(), source })?
                .filter_map(|e| e.ok().map(|e| e.path()))
                .filter(|p| p.extension().is_some_and(|e| e.eq_ignore_ascii_case("toml")))
                .collect();
            files.sort();
            for path in files {
                let id = path.file_stem().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
                if !valid_id(&id) {
                    return Err(Error::Config(format!(
                        "{}: profile id '{id}' must use a-z, 0-9, '_' or '-'",
                        path.display()
                    )));
                }
                let text = std::fs::read_to_string(&path).map_err(|source| Error::Io { path: path.clone(), source })?;
                let mut table: toml::Table =
                    toml::from_str(&text).map_err(|source| Error::Toml { path: path.clone(), source })?;
                let based_on = take_string(&mut table, "based_on")?;
                let source = take_string(&mut table, "source")?;
                let params = params_from_table(table, &path.display().to_string())?;
                lib.entries.push(LaunchVehicleEntry { id, is_preset: false, based_on, source, params, path: Some(path) });
            }
        }
        lib.check_ids()?;
        Ok(lib)
    }

    fn check_ids(&self) -> Result<()> {
        let mut seen = HashSet::new();
        for e in &self.entries {
            if !seen.insert(e.id.as_str()) {
                return Err(Error::Config(format!(
                    "rocket profile id '{}' is defined twice (custom profiles cannot reuse preset ids)",
                    e.id
                )));
            }
        }
        Ok(())
    }

    pub fn entries(&self) -> &[LaunchVehicleEntry] {
        &self.entries
    }

    pub fn get(&self, id: &str) -> Option<&LaunchVehicleEntry> {
        self.entries.iter().find(|e| e.id == id)
    }

    fn chain(&self, id: &str) -> Result<(LaunchVehicleParams, Vec<&'static str>, Vec<String>)> {
        let mut stack: Vec<&LaunchVehicleEntry> = Vec::new();
        let mut cur = id;
        loop {
            let e = self.get(cur).ok_or_else(|| Error::Config(format!("unknown rocket profile '{cur}'")))?;
            if stack.iter().any(|x| x.id == e.id) {
                return Err(Error::Config(format!("rocket profile '{id}': based_on cycle")));
            }
            stack.push(e);
            match &e.based_on {
                Some(b) => cur = b,
                None => break,
            }
        }
        let mut params = LaunchVehicleParams::default();
        let mut overridden = Vec::new();
        let mut sources = Vec::new();
        for e in stack.iter().rev() {
            if let Some(s) = &e.source {
                sources.push(format!("{}: {s}", e.id));
            }
            if !e.is_preset {
                for f in e.params.set_fields() {
                    if !overridden.contains(&f) {
                        overridden.push(f);
                    }
                }
            }
            params = params.overlay(&e.params);
        }
        Ok((params, overridden, sources))
    }

    /// Собирает ракету: профиль `id` + переопределения этого запуска.
    pub fn resolve(&self, id: &str, overrides: &LaunchVehicleParams) -> Result<ResolvedLaunchVehicle> {
        let entry = self.get(id).ok_or_else(|| {
            Error::Config(format!(
                "unknown rocket '{id}'; available: {}",
                self.entries.iter().map(|e| e.id.as_str()).collect::<Vec<_>>().join(", ")
            ))
        })?;
        let (base, mut overridden, sources) = self.chain(id)?;
        let params = base.overlay(overrides);
        for f in overrides.set_fields() {
            if !overridden.contains(&f) {
                overridden.push(f);
            }
        }
        let base_name = base.name.clone().unwrap_or_else(|| id.to_string());
        let custom = !(entry.is_preset && overrides.is_empty());
        let config = LaunchVehicleConfig::from_params(&params, &base_name, custom)?;
        let status_label = if entry.is_preset && overrides.is_empty() {
            format!("Preset: {}", config.name)
        } else if !overrides.is_empty() {
            format!("Custom (на базе {base_name})")
        } else {
            match &entry.based_on {
                Some(b) => {
                    let b_name = self.chain(b).ok().and_then(|c| c.0.name).unwrap_or_else(|| b.clone());
                    format!("Custom (на базе {b_name})")
                }
                None => "Custom".to_string(),
            }
        };
        let errors = config.validate();
        Ok(ResolvedLaunchVehicle { id: id.to_string(), status_label, config, overridden, sources, errors })
    }

    /// Сохраняет `new_id` = `based_on` + `overrides` в папку своих профилей.
    pub fn save_custom(&self, new_id: &str, based_on: &str, overrides: &LaunchVehicleParams) -> Result<PathBuf> {
        if !valid_id(new_id) {
            return Err(Error::Config(format!("profile id '{new_id}' must use a-z, 0-9, '_' or '-'")));
        }
        if self.get(new_id).is_some_and(|e| e.is_preset) {
            return Err(Error::Config(format!("'{new_id}' is a preset id; choose another name")));
        }
        let base = self.get(based_on).ok_or_else(|| Error::Config(format!("unknown rocket profile '{based_on}'")))?;
        let (based, params) = if new_id == based_on {
            (base.based_on.clone(), base.params.overlay(overrides))
        } else {
            (Some(based_on.to_string()), overrides.clone())
        };
        if self.custom_dir.as_os_str().is_empty() {
            return Err(Error::Config("no folder for custom rocket profiles".into()));
        }
        std::fs::create_dir_all(&self.custom_dir).map_err(|source| Error::Io { path: self.custom_dir.clone(), source })?;
        let mut table = toml::Table::new();
        if let Some(b) = based {
            table.insert("based_on".into(), toml::Value::String(b));
        }
        if let toml::Value::Table(t) =
            toml::Value::try_from(&params).map_err(|e| Error::Config(format!("serialize profile: {e}")))?
        {
            table.extend(t);
        }
        let text = format!(
            "# Пользовательский профиль ракеты-носителя (создан программой).\n\
             # based_on - основа; указаны только изменённые параметры.\n{}",
            toml::to_string(&table).map_err(|e| Error::Config(format!("serialize profile: {e}")))?
        );
        let path = self.custom_dir.join(format!("{new_id}.toml"));
        std::fs::write(&path, text).map_err(|source| Error::Io { path: path.clone(), source })?;
        Ok(path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spaceport::SpaceportCatalog;
    use approx::assert_abs_diff_eq;

    fn none() -> LaunchVehicleParams {
        LaunchVehicleParams::default()
    }

    #[test]
    fn presets_are_valid_and_sourced() {
        let lib = LaunchVehicleLibrary::presets_only().unwrap();
        assert!(lib.entries().len() >= 4);
        let text = std::fs::read_to_string("data/spaceports.toml").unwrap();
        let ports = SpaceportCatalog::parse(&text, Path::new("data/spaceports.toml")).unwrap();
        for e in lib.entries() {
            let r = lib.resolve(&e.id, &none()).unwrap();
            assert!(r.errors.is_empty(), "{}: {:?}", e.id, r.errors);
            assert!(!r.sources.is_empty(), "{} has no source", e.id);
            for p in &r.config.spaceports {
                assert!(ports.get(p).is_some(), "{}: unknown spaceport {p}", e.id);
            }
        }
    }

    #[test]
    fn capability_equals_reference_at_rated_payload_and_grows_when_lighter() {
        let lib = LaunchVehicleLibrary::presets_only().unwrap();
        let r = lib.resolve("long_march_5", &none()).unwrap();
        let rated = r.config.capability(r.config.payload_tli_kg);
        assert_abs_diff_eq!(rated.max_tli_delta_v_m_s, REFERENCE_TLI_DELTA_V_M_S, epsilon = 1e-9);
        assert!(rated.within_rated_payload);
        let light = r.config.capability(4000.0);
        assert!(light.max_tli_delta_v_m_s > REFERENCE_TLI_DELTA_V_M_S + 300.0);
        let heavy = r.config.capability(12_000.0);
        assert!(!heavy.within_rated_payload);
        assert!(heavy.max_tli_delta_v_m_s < REFERENCE_TLI_DELTA_V_M_S);
    }

    #[test]
    fn tli_burn_time_is_physical() {
        let lib = LaunchVehicleLibrary::presets_only().unwrap();
        let r = lib.resolve("sls_block1", &none()).unwrap();
        let b = r.config.tli_burn(27_000.0, 3150.0);
        assert!(b.propellant_used_kg > 0.0, "used {}", b.propellant_used_kg);
        // ICPS выдаёт разгон к Луне за ~20 минут.
        assert!(b.duration_s > 600.0 && b.duration_s < 1500.0, "burn {} s", b.duration_s);
        // С круговой орбиты ICPS на полный разгон 27 т не хватает: она стартует с вытянутой орбиты.
        assert!(b.exceeds_stage_propellant);
        let light = r.config.tli_burn(10_000.0, 3150.0);
        assert!(!light.exceeds_stage_propellant);
    }

    #[test]
    fn override_save_and_reload() {
        let dir = std::env::temp_dir().join(format!("lwp-lv-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let lib = LaunchVehicleLibrary::load(&dir).unwrap();
        let mut o = none();
        o.set_from_str("payload_tli_kg=9400").unwrap();
        assert!(o.set_from_str("spaceports=wenchang,xichang").is_err());
        let r = lib.resolve("long_march_5", &o).unwrap();
        assert_eq!(r.status_label, "Custom (на базе Long March 5)");
        assert!(r.config.any_spaceport, "custom rockets launch from any spaceport");
        assert!(!lib.resolve("long_march_5", &none()).unwrap().config.any_spaceport);
        assert_eq!(r.config.upper_stage_isp_s, 442.6);
        lib.save_custom("lm5_upper", "long_march_5", &o).unwrap();
        let lib2 = LaunchVehicleLibrary::load(&dir).unwrap();
        let r2 = lib2.resolve("lm5_upper", &none()).unwrap();
        assert_eq!(r2.config.payload_tli_kg, 9400.0);
        assert!(r2.config.any_spaceport);
        assert_eq!(r2.status_label, "Custom (на базе Long March 5)");
        assert!(lib2.save_custom("ariane_64", "ariane_64", &o).is_err());
    }

    #[test]
    fn standalone_custom_rocket_needs_no_spaceports() {
        let dir = std::env::temp_dir().join(format!("lwp-lv-solo-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("my_rocket.toml"),
            "name = \"My rocket\"\npayload_tli_kg = 5000.0\nupper_stage_dry_mass_kg = 3000.0\n\
             upper_stage_propellant_kg = 20000.0\nupper_stage_thrust_n = 100000.0\nupper_stage_isp_s = 450.0\n",
        )
        .unwrap();
        let lib = LaunchVehicleLibrary::load(&dir).unwrap();
        let r = lib.resolve("my_rocket", &none()).unwrap();
        assert!(r.errors.is_empty(), "{:?}", r.errors);
        assert!(r.config.any_spaceport);
        assert_eq!(r.status_label, "Custom");
    }

    #[test]
    fn input_errors() {
        let mut o = none();
        assert!(o.set_from_str("isp=1").is_err());
        assert!(o.set_from_str("upper_stage_isp_s=x").is_err());
        assert!(o.set_from_str("payload_tli_kg").is_err());
        let lib = LaunchVehicleLibrary::presets_only().unwrap();
        let mut bad = none();
        bad.set_from_str("upper_stage_isp_s=0").unwrap();
        assert!(!lib.resolve("ariane_64", &bad).unwrap().errors.is_empty());
    }
}
