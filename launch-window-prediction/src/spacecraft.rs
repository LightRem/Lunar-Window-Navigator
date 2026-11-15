//! Конфигурация лунного аппарата (Spacecraft Configuration).
//!
//! Реализует справочник «Подсистема конфигурации лунного аппарата»:
//!
//! - **Шаблоны (presets)** — встроенная библиотека `data/spacecraft/presets.toml`
//!   с характеристиками реальных аппаратов (источник указан у каждого).
//! - **Переопределение (override)** — любой параметр шаблона можно изменить; профиль тогда
//!   получает статус `Custom (на базе <шаблон>)`, остальные значения остаются из шаблона.
//! - **Сохранение** — свои профили хранятся в `data/spacecraft/custom/<id>.toml`
//!   (только `based_on` и изменённые поля) и появляются в списке доступных профилей.
//! - **Уровни параметров** — базовый (массы, тяга, Isp) обязателен, продвинутый
//!   (дросселирование, энергетика, связь) — по желанию.
//! - **Производные величины** — стартовая масса, Δv по Циолковскому, тяговооружённость
//!   на Луне, время работы в тени.
//! - **Проверки** — критические ошибки и предупреждения по разделу 4 справочника.

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

/// Стандартное ускорение свободного падения, м/с².
pub const G0: f64 = 9.80665;
/// Ускорение свободного падения на поверхности Луны, м/с².
pub const G_MOON: f64 = 1.62;
/// Δv, нужная для выхода с поверхности Луны на низкую окололунную орбиту (~100 км)
/// с учётом гравитационных потерь, м/с (справочник, раздел 4).
pub const MIN_ASCENT_DELTA_V_M_S: f64 = 1870.0;

/// Встроенная библиотека шаблонов.
const PRESETS_TOML: &str = include_str!("../data/spacecraft/presets.toml");

// ---------------------------------------------------------------------------
// Параметры
// ---------------------------------------------------------------------------

/// Частичный набор параметров: шаблон, пользовательский профиль или переопределения.
/// Неуказанное поле означает «взять из основы».
#[derive(Debug, Clone, Default, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SpacecraftParams {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    // Базовый уровень.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dry_mass_kg: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub propellant_mass_kg: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub payload_mass_kg: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_thrust_n: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub isp_s: Option<f64>,
    // Продвинутый уровень.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_throttle_percent: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub solar_panel_area_m2: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub solar_panel_efficiency: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub battery_capacity_wh: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub idle_power_consumption_w: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tx_power_w: Option<f64>,
}

/// Имена числовых полей (для ввода `поле=значение` и отчёта о переопределениях).
pub const NUMERIC_FIELDS: [&str; 11] = [
    "dry_mass_kg",
    "propellant_mass_kg",
    "payload_mass_kg",
    "max_thrust_n",
    "isp_s",
    "min_throttle_percent",
    "solar_panel_area_m2",
    "solar_panel_efficiency",
    "battery_capacity_wh",
    "idle_power_consumption_w",
    "tx_power_w",
];

impl SpacecraftParams {
    fn field_mut(&mut self, field: &str) -> Option<&mut Option<f64>> {
        Some(match field {
            "dry_mass_kg" => &mut self.dry_mass_kg,
            "propellant_mass_kg" => &mut self.propellant_mass_kg,
            "payload_mass_kg" => &mut self.payload_mass_kg,
            "max_thrust_n" => &mut self.max_thrust_n,
            "isp_s" => &mut self.isp_s,
            "min_throttle_percent" => &mut self.min_throttle_percent,
            "solar_panel_area_m2" => &mut self.solar_panel_area_m2,
            "solar_panel_efficiency" => &mut self.solar_panel_efficiency,
            "battery_capacity_wh" => &mut self.battery_capacity_wh,
            "idle_power_consumption_w" => &mut self.idle_power_consumption_w,
            "tx_power_w" => &mut self.tx_power_w,
            _ => return None,
        })
    }

    fn field(&self, field: &str) -> Option<f64> {
        self.clone().field_mut(field).and_then(|v| *v)
    }

    /// Применяет строку `поле=значение` (например, `payload_mass_kg=500` или `name=Мой аппарат`).
    pub fn set_from_str(&mut self, assignment: &str) -> Result<()> {
        let (key, value) = assignment
            .split_once('=')
            .ok_or_else(|| Error::Config(format!("spacecraft override '{assignment}' must be field=value")))?;
        let (key, value) = (key.trim(), value.trim());
        if key == "name" {
            self.name = Some(value.to_string());
            return Ok(());
        }
        let slot = self.field_mut(key).ok_or_else(|| {
            Error::Config(format!(
                "unknown spacecraft field '{key}'; known: name, {}",
                NUMERIC_FIELDS.join(", ")
            ))
        })?;
        let v: f64 = value
            .parse()
            .map_err(|_| Error::Config(format!("spacecraft field '{key}': '{value}' is not a number")))?;
        *slot = Some(v);
        Ok(())
    }

    /// Поля, заданные в этом наборе (включая `name`).
    pub fn set_fields(&self) -> Vec<&'static str> {
        let mut out = Vec::new();
        if self.name.is_some() {
            out.push("name");
        }
        out.extend(NUMERIC_FIELDS.iter().copied().filter(|f| self.field(f).is_some()));
        out
    }

    pub fn is_empty(&self) -> bool {
        self.set_fields().is_empty()
    }

    /// Накладывает `over` поверх `self`.
    pub fn overlay(&self, over: &SpacecraftParams) -> SpacecraftParams {
        let mut out = self.clone();
        if over.name.is_some() {
            out.name = over.name.clone();
        }
        for f in NUMERIC_FIELDS {
            if let Some(v) = over.field(f) {
                *out.field_mut(f).expect("known field") = Some(v);
            }
        }
        out
    }
}

/// Полная конфигурация аппарата (базовые поля обязательны).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SpacecraftConfig {
    pub name: String,
    pub dry_mass_kg: f64,
    pub propellant_mass_kg: f64,
    pub payload_mass_kg: f64,
    pub max_thrust_n: f64,
    pub isp_s: f64,
    pub min_throttle_percent: Option<f64>,
    pub solar_panel_area_m2: Option<f64>,
    pub solar_panel_efficiency: Option<f64>,
    pub battery_capacity_wh: Option<f64>,
    pub idle_power_consumption_w: Option<f64>,
    pub tx_power_w: Option<f64>,
}

impl SpacecraftConfig {
    fn from_params(p: &SpacecraftParams, fallback_name: &str) -> Result<Self> {
        let missing: Vec<&str> = ["dry_mass_kg", "propellant_mass_kg", "payload_mass_kg", "max_thrust_n", "isp_s"]
            .into_iter()
            .filter(|f| p.field(f).is_none())
            .collect();
        if !missing.is_empty() {
            return Err(Error::Config(format!(
                "spacecraft '{fallback_name}': missing base parameters: {}",
                missing.join(", ")
            )));
        }
        Ok(Self {
            name: p.name.clone().unwrap_or_else(|| fallback_name.to_string()),
            dry_mass_kg: p.dry_mass_kg.unwrap_or_default(),
            propellant_mass_kg: p.propellant_mass_kg.unwrap_or_default(),
            payload_mass_kg: p.payload_mass_kg.unwrap_or_default(),
            max_thrust_n: p.max_thrust_n.unwrap_or_default(),
            isp_s: p.isp_s.unwrap_or_default(),
            min_throttle_percent: p.min_throttle_percent,
            solar_panel_area_m2: p.solar_panel_area_m2,
            solar_panel_efficiency: p.solar_panel_efficiency,
            battery_capacity_wh: p.battery_capacity_wh,
            idle_power_consumption_w: p.idle_power_consumption_w,
            tx_power_w: p.tx_power_w,
        })
    }

    /// Производные величины (раздел 3 справочника).
    pub fn derived(&self) -> Derived {
        let m0 = self.dry_mass_kg + self.propellant_mass_kg + self.payload_mass_kg;
        let m_final = self.dry_mass_kg + self.payload_mass_kg;
        let delta_v = if m_final > 0.0 && m0 > 0.0 { self.isp_s * G0 * (m0 / m_final).ln() } else { f64::NAN };
        let twr = if m0 > 0.0 { self.max_thrust_n / (m0 * G_MOON) } else { f64::NAN };
        let eclipse = match (self.battery_capacity_wh, self.idle_power_consumption_w) {
            (Some(e), Some(p)) if p > 0.0 => Some(e / p),
            _ => None,
        };
        Derived {
            initial_mass_kg: m0,
            final_mass_kg: m_final,
            delta_v_m_s: delta_v,
            twr_moon: twr,
            eclipse_autonomy_h: eclipse,
        }
    }

    /// Проверки (раздел 4 справочника).
    pub fn validate(&self) -> Vec<Issue> {
        let mut issues = Vec::new();
        let mut err = |field: &'static str, msg: String| {
            issues.push(Issue { severity: Severity::Error, field: Some(field), message: msg })
        };
        let positive = [
            ("dry_mass_kg", self.dry_mass_kg),
            ("propellant_mass_kg", self.propellant_mass_kg),
            ("max_thrust_n", self.max_thrust_n),
            ("isp_s", self.isp_s),
        ];
        for (f, v) in positive {
            if !(v > 0.0 && v.is_finite()) {
                err(f, format!("{f} must be > 0 (got {v})"));
            }
        }
        if !(self.payload_mass_kg >= 0.0 && self.payload_mass_kg.is_finite()) {
            err("payload_mass_kg", format!("payload_mass_kg must be >= 0 (got {})", self.payload_mass_kg));
        }
        if let Some(e) = self.solar_panel_efficiency
            && !(e > 0.0 && e <= 1.0)
        {
            err("solar_panel_efficiency", format!("solar_panel_efficiency must be in (0, 1] (got {e})"));
        }
        if let Some(t) = self.min_throttle_percent
            && !(0.0..=100.0).contains(&t)
        {
            err("min_throttle_percent", format!("min_throttle_percent must be in [0, 100] (got {t})"));
        }
        for (f, v) in [
            ("solar_panel_area_m2", self.solar_panel_area_m2),
            ("battery_capacity_wh", self.battery_capacity_wh),
            ("idle_power_consumption_w", self.idle_power_consumption_w),
            ("tx_power_w", self.tx_power_w),
        ] {
            if let Some(v) = v
                && !(v >= 0.0 && v.is_finite())
            {
                err(f, format!("{f} must be >= 0 (got {v})"));
            }
        }
        // Остальные проверки имеют смысл только при корректных базовых величинах.
        if issues.iter().any(|i| i.severity == Severity::Error) {
            return issues;
        }
        let d = self.derived();
        if d.twr_moon <= 1.0 {
            issues.push(Issue {
                severity: Severity::Error,
                field: Some("max_thrust_n"),
                message: format!(
                    "lunar thrust-to-weight {:.2} <= 1.0: the engine cannot lift off or hover on the Moon",
                    d.twr_moon
                ),
            });
        }
        if d.delta_v_m_s < MIN_ASCENT_DELTA_V_M_S {
            issues.push(Issue {
                severity: Severity::Warning,
                field: Some("propellant_mass_kg"),
                message: format!(
                    "delta-v {:.0} m/s < {:.0} m/s: not enough to reach low lunar orbit (~100 km) from the surface",
                    d.delta_v_m_s, MIN_ASCENT_DELTA_V_M_S
                ),
            });
        }
        issues
    }
}

/// Производные величины.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Derived {
    /// m0 = сухая + топливо + ПН, кг.
    pub initial_mass_kg: f64,
    /// Масса после выработки топлива, кг.
    pub final_mass_kg: f64,
    /// Характеристическая скорость (Циолковский), м/с.
    pub delta_v_m_s: f64,
    /// Тяговооружённость на Луне.
    pub twr_moon: f64,
    /// Время работы в тени на аккумуляторах, ч (если заданы ёмкость и потребление).
    pub eclipse_autonomy_h: Option<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Severity {
    /// Критическая ошибка: аппарат с такими параметрами не годится.
    Error,
    Warning,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Issue {
    pub severity: Severity,
    pub field: Option<&'static str>,
    pub message: String,
}

// ---------------------------------------------------------------------------
// Профили и библиотека
// ---------------------------------------------------------------------------

/// Откуда профиль.
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum ProfileStatus {
    /// Шаблон без изменений.
    Preset { id: String },
    /// Пользовательский профиль (сохранённый или с переопределениями в этом запуске).
    Custom { id: Option<String>, based_on: Option<String> },
}

/// Запись библиотеки: шаблон или сохранённый профиль.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ProfileEntry {
    pub id: String,
    pub is_preset: bool,
    pub based_on: Option<String>,
    pub source: Option<String>,
    pub params: SpacecraftParams,
    pub path: Option<PathBuf>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct PresetFile {
    #[serde(default)]
    preset: Vec<toml::Table>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct CustomHeader {
    #[serde(default)]
    based_on: Option<String>,
    #[serde(default)]
    source: Option<String>,
}

/// Итог: аппарат, готовый к расчётам.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ResolvedSpacecraft {
    pub status: ProfileStatus,
    /// Человекочитаемый статус: `Preset: …`, `Custom (на базе …)`, `Custom`.
    pub status_label: String,
    pub config: SpacecraftConfig,
    /// Поля, изменённые относительно основы (в сохранённом профиле и в этом запуске).
    pub overridden: Vec<&'static str>,
    /// Источники данных по цепочке основ.
    pub sources: Vec<String>,
    pub derived: Derived,
    pub issues: Vec<Issue>,
}

impl ResolvedSpacecraft {
    pub fn has_errors(&self) -> bool {
        self.issues.iter().any(|i| i.severity == Severity::Error)
    }
}

/// Библиотека: встроенные шаблоны + профили из папки.
#[derive(Debug, Clone)]
pub struct SpacecraftLibrary {
    entries: Vec<ProfileEntry>,
    custom_dir: PathBuf,
}

pub(crate) fn params_from_table<T: serde::de::DeserializeOwned>(table: toml::Table, origin: &str) -> Result<T> {
    toml::Value::Table(table)
        .try_into()
        .map_err(|e| Error::Config(format!("{origin}: {e}")))
}

pub(crate) fn valid_id(id: &str) -> bool {
    !id.is_empty() && id.len() <= 64 && id.chars().all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_' || c == '-')
}

impl SpacecraftLibrary {
    /// Только встроенные шаблоны.
    pub fn presets_only() -> Result<Self> {
        let file: PresetFile = toml::from_str(PRESETS_TOML)
            .map_err(|source| Error::Toml { path: PathBuf::from("data/spacecraft/presets.toml"), source })?;
        let mut entries = Vec::new();
        for mut row in file.preset {
            let take_str = |row: &mut toml::Table, key: &str| -> Result<String> {
                match row.remove(key) {
                    Some(toml::Value::String(s)) => Ok(s),
                    _ => Err(Error::Config(format!("preset: '{key}' must be a string"))),
                }
            };
            let id = take_str(&mut row, "id")?;
            let source = take_str(&mut row, "source")?;
            let params = params_from_table(row, &format!("preset '{id}'"))?;
            entries.push(ProfileEntry { id, is_preset: true, based_on: None, source: Some(source), params, path: None });
        }
        let lib = Self { entries, custom_dir: PathBuf::new() };
        lib.check_ids()?;
        Ok(lib)
    }

    /// Шаблоны + профили из `custom_dir` (`*.toml`; папки может не быть).
    /// Синхронно; из async вызывать через `spawn_blocking`.
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
                lib.entries.push(Self::read_custom(&path)?);
            }
        }
        lib.check_ids()?;
        Ok(lib)
    }

    fn read_custom(path: &Path) -> Result<ProfileEntry> {
        let id = path.file_stem().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
        if !valid_id(&id) {
            return Err(Error::Config(format!(
                "{}: profile id '{id}' must use a-z, 0-9, '_' or '-'",
                path.display()
            )));
        }
        let text = std::fs::read_to_string(path).map_err(|source| Error::Io { path: path.to_path_buf(), source })?;
        let mut table: toml::Table = toml::from_str(&text)
            .map_err(|source| Error::Toml { path: path.to_path_buf(), source })?;
        let mut header = toml::Table::new();
        for key in ["based_on", "source"] {
            if let Some(v) = table.remove(key) {
                header.insert(key.into(), v);
            }
        }
        let header: CustomHeader = toml::Value::Table(header)
            .try_into()
            .map_err(|e| Error::Config(format!("{}: {e}", path.display())))?;
        let params = params_from_table(table, &path.display().to_string())?;
        Ok(ProfileEntry {
            id,
            is_preset: false,
            based_on: header.based_on,
            source: header.source,
            params,
            path: Some(path.to_path_buf()),
        })
    }

    fn check_ids(&self) -> Result<()> {
        let mut seen = HashSet::new();
        for e in &self.entries {
            if !seen.insert(e.id.as_str()) {
                return Err(Error::Config(format!(
                    "spacecraft profile id '{}' is defined twice (custom profiles cannot reuse preset ids)",
                    e.id
                )));
            }
        }
        Ok(())
    }

    pub fn entries(&self) -> &[ProfileEntry] {
        &self.entries
    }

    pub fn get(&self, id: &str) -> Option<&ProfileEntry> {
        self.entries.iter().find(|e| e.id == id)
    }

    /// Параметры профиля с учётом цепочки `based_on`, поля, изменённые относительно
    /// шаблона-основы, и источники по цепочке.
    fn chain(&self, id: &str) -> Result<(SpacecraftParams, Vec<&'static str>, Vec<String>, Option<String>)> {
        let mut stack = Vec::new();
        let mut cur = id;
        loop {
            let e = self
                .get(cur)
                .ok_or_else(|| Error::Config(format!("unknown spacecraft profile '{cur}'")))?;
            if stack.iter().any(|x: &&ProfileEntry| x.id == e.id) {
                return Err(Error::Config(format!("spacecraft profile '{id}': based_on cycle")));
            }
            stack.push(e);
            match &e.based_on {
                Some(b) => cur = b,
                None => break,
            }
        }
        let root = stack.last().expect("non-empty");
        let mut params = SpacecraftParams::default();
        let mut overridden: Vec<&'static str> = Vec::new();
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
        let root_preset = root.is_preset.then(|| root.id.clone());
        Ok((params, overridden, sources, root_preset))
    }

    /// Собирает аппарат: профиль `id` (шаблон или сохранённый) плюс переопределения этого запуска.
    pub fn resolve(&self, id: &str, overrides: &SpacecraftParams) -> Result<ResolvedSpacecraft> {
        let entry = self.get(id).ok_or_else(|| {
            Error::Config(format!(
                "unknown spacecraft profile '{id}'; available: {}",
                self.entries.iter().map(|e| e.id.as_str()).collect::<Vec<_>>().join(", ")
            ))
        })?;
        let (base, mut overridden, sources, _root) = self.chain(id)?;
        let params = base.overlay(overrides);
        for f in overrides.set_fields() {
            if !overridden.contains(&f) {
                overridden.push(f);
            }
        }
        let fallback_name = base.name.clone().unwrap_or_else(|| id.to_string());
        let config = SpacecraftConfig::from_params(&params, &fallback_name)?;

        let (status, status_label) = if entry.is_preset && overrides.is_empty() {
            (ProfileStatus::Preset { id: id.to_string() }, format!("Preset: {}", config.name))
        } else if entry.is_preset {
            let base_name = base.name.clone().unwrap_or_else(|| id.to_string());
            (
                ProfileStatus::Custom { id: None, based_on: Some(id.to_string()) },
                format!("Custom (на базе {base_name})"),
            )
        } else {
            let based_on = if overrides.is_empty() { entry.based_on.clone() } else { Some(id.to_string()) };
            let label = match &based_on {
                Some(b) => {
                    let b_name = self
                        .get(b)
                        .and_then(|e| e.params.name.clone())
                        .or_else(|| self.chain(b).ok().and_then(|c| c.0.name))
                        .unwrap_or_else(|| b.clone());
                    format!("Custom (на базе {b_name})")
                }
                None => "Custom".to_string(),
            };
            let custom_id = overrides.is_empty().then(|| id.to_string());
            (ProfileStatus::Custom { id: custom_id, based_on }, label)
        };

        let derived = config.derived();
        let issues = config.validate();
        Ok(ResolvedSpacecraft { status, status_label, config, overridden, sources, derived, issues })
    }

    /// Сохраняет профиль `new_id` = основа `based_on` + `overrides` в папку пользовательских профилей.
    /// Возвращает путь к файлу. Существующий пользовательский профиль перезаписывается,
    /// шаблон — нет.
    pub fn save_custom(&self, new_id: &str, based_on: &str, overrides: &SpacecraftParams) -> Result<PathBuf> {
        if !valid_id(new_id) {
            return Err(Error::Config(format!("profile id '{new_id}' must use a-z, 0-9, '_' or '-'")));
        }
        if self.get(new_id).is_some_and(|e| e.is_preset) {
            return Err(Error::Config(format!("'{new_id}' is a preset id; choose another name")));
        }
        if new_id == based_on {
            // Изменение существующего профиля: основа — его собственная основа.
            let e = self.get(based_on).expect("checked by caller");
            let merged = e.params.overlay(overrides);
            return self.write_custom(new_id, e.based_on.as_deref(), e.source.as_deref(), &merged);
        }
        self.get(based_on)
            .ok_or_else(|| Error::Config(format!("unknown spacecraft profile '{based_on}'")))?;
        self.write_custom(new_id, Some(based_on), None, overrides)
    }

    fn write_custom(
        &self,
        id: &str,
        based_on: Option<&str>,
        source: Option<&str>,
        params: &SpacecraftParams,
    ) -> Result<PathBuf> {
        if self.custom_dir.as_os_str().is_empty() {
            return Err(Error::Config("no folder for custom spacecraft profiles".into()));
        }
        std::fs::create_dir_all(&self.custom_dir)
            .map_err(|source| Error::Io { path: self.custom_dir.clone(), source })?;
        let mut table = toml::Table::new();
        if let Some(b) = based_on {
            table.insert("based_on".into(), toml::Value::String(b.into()));
        }
        if let Some(s) = source {
            table.insert("source".into(), toml::Value::String(s.into()));
        }
        let body = toml::Value::try_from(params).map_err(|e| Error::Config(format!("serialize profile: {e}")))?;
        if let toml::Value::Table(t) = body {
            table.extend(t);
        }
        let text = format!(
            "# Пользовательский профиль лунного аппарата (создан программой).\n\
             # based_on - основа; указаны только изменённые параметры.\n{}",
            toml::to_string(&table).map_err(|e| Error::Config(format!("serialize profile: {e}")))?
        );
        let path = self.custom_dir.join(format!("{id}.toml"));
        std::fs::write(&path, text).map_err(|source| Error::Io { path: path.clone(), source })?;
        Ok(path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn lib_with(files: &[(&str, &str)]) -> (SpacecraftLibrary, PathBuf) {
        let dir = std::env::temp_dir().join(format!(
            "lwp-sc-{}-{}",
            std::process::id(),
            files.first().map(|f| f.0).unwrap_or("none")
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        for (name, text) in files {
            std::fs::write(dir.join(format!("{name}.toml")), text).unwrap();
        }
        (SpacecraftLibrary::load(&dir).unwrap(), dir)
    }

    #[test]
    fn presets_load_and_are_valid() {
        let lib = SpacecraftLibrary::presets_only().unwrap();
        assert!(lib.entries().len() >= 3);
        for e in lib.entries() {
            let r = lib.resolve(&e.id, &SpacecraftParams::default()).unwrap();
            assert!(!r.has_errors(), "{}: {:?}", e.id, r.issues);
            assert!(matches!(r.status, ProfileStatus::Preset { .. }));
            assert!(r.overridden.is_empty());
            assert!(!r.sources.is_empty(), "{} has no source", e.id);
        }
    }

    #[test]
    fn derived_values_follow_the_formulas() {
        let lib = SpacecraftLibrary::presets_only().unwrap();
        let r = lib.resolve("apollo_lm_ascent", &SpacecraftParams::default()).unwrap();
        let c = &r.config;
        let m0 = c.dry_mass_kg + c.propellant_mass_kg + c.payload_mass_kg;
        assert_abs_diff_eq!(r.derived.initial_mass_kg, m0, epsilon = 1e-9);
        let dv = c.isp_s * G0 * (m0 / (c.dry_mass_kg + c.payload_mass_kg)).ln();
        assert_abs_diff_eq!(r.derived.delta_v_m_s, dv, epsilon = 1e-9);
        assert_abs_diff_eq!(r.derived.twr_moon, c.max_thrust_n / (m0 * G_MOON), epsilon = 1e-12);
        // Взлётная ступень Apollo: ~2.2 км/с, тяговооружённость на Луне ~2.
        assert!(r.derived.delta_v_m_s > 2100.0 && r.derived.delta_v_m_s < 2400.0);
        assert!(r.derived.twr_moon > 1.8 && r.derived.twr_moon < 2.5);
        assert!(r.derived.eclipse_autonomy_h.is_none());
    }

    #[test]
    fn override_switches_status_and_keeps_other_values() {
        let lib = SpacecraftLibrary::presets_only().unwrap();
        let mut o = SpacecraftParams::default();
        o.set_from_str("payload_mass_kg=500").unwrap();
        let r = lib.resolve("apollo_lm_ascent", &o).unwrap();
        assert_eq!(r.status_label, "Custom (на базе Apollo LM Ascent Stage)");
        assert_eq!(r.overridden, vec!["payload_mass_kg"]);
        assert_eq!(r.config.payload_mass_kg, 500.0);
        assert_eq!(r.config.isp_s, 311.0);
        assert_eq!(r.config.name, "Apollo LM Ascent Stage");
    }

    #[test]
    fn guardrails_errors_and_warnings() {
        let lib = SpacecraftLibrary::presets_only().unwrap();
        let mut o = SpacecraftParams::default();
        o.set_from_str("max_thrust_n=5000").unwrap(); // TWR < 1
        let r = lib.resolve("apollo_lm_ascent", &o).unwrap();
        assert!(r.has_errors());
        let mut o = SpacecraftParams::default();
        o.set_from_str("propellant_mass_kg=1000").unwrap(); // Δv < 1870
        let r = lib.resolve("apollo_lm_ascent", &o).unwrap();
        assert!(!r.has_errors());
        assert!(r.issues.iter().any(|i| i.severity == Severity::Warning));
        let mut o = SpacecraftParams::default();
        o.set_from_str("solar_panel_efficiency=1.5").unwrap();
        assert!(lib.resolve("apollo_lm_ascent", &o).unwrap().has_errors());
        let mut o = SpacecraftParams::default();
        o.set_from_str("payload_mass_kg=-1").unwrap();
        assert!(lib.resolve("apollo_lm_ascent", &o).unwrap().has_errors());
    }

    #[test]
    fn override_input_errors() {
        let mut o = SpacecraftParams::default();
        assert!(o.set_from_str("isp").is_err());
        assert!(o.set_from_str("isp_sec=300").is_err());
        assert!(o.set_from_str("isp_s=abc").is_err());
        o.set_from_str("name=My lander").unwrap();
        assert_eq!(o.name.as_deref(), Some("My lander"));
    }

    #[test]
    fn eclipse_autonomy_from_battery_and_idle_power() {
        let lib = SpacecraftLibrary::presets_only().unwrap();
        let mut o = SpacecraftParams::default();
        o.set_from_str("battery_capacity_wh=2000").unwrap();
        o.set_from_str("idle_power_consumption_w=100").unwrap();
        let r = lib.resolve("apollo_lm_descent", &o).unwrap();
        assert_abs_diff_eq!(r.derived.eclipse_autonomy_h.unwrap(), 20.0, epsilon = 1e-12);
    }

    #[test]
    fn custom_profiles_chain_and_status() {
        let (lib, _) = lib_with(&[
            ("lm_heavy", "based_on = \"apollo_lm_descent\"\nname = \"Heavy LM\"\npayload_mass_kg = 5000.0\n"),
            ("lm_heavy_hot", "based_on = \"lm_heavy\"\nisp_s = 320.0\n"),
        ]);
        let r = lib.resolve("lm_heavy", &SpacecraftParams::default()).unwrap();
        assert_eq!(r.status_label, "Custom (на базе Apollo LM Descent Stage)");
        assert_eq!(r.config.payload_mass_kg, 5000.0);
        assert_eq!(r.config.dry_mass_kg, 2767.0);
        let r2 = lib.resolve("lm_heavy_hot", &SpacecraftParams::default()).unwrap();
        assert_eq!(r2.config.isp_s, 320.0);
        assert_eq!(r2.config.payload_mass_kg, 5000.0);
        assert_eq!(r2.status_label, "Custom (на базе Heavy LM)");
        assert!(r2.overridden.contains(&"isp_s") && r2.overridden.contains(&"payload_mass_kg"));
    }

    #[test]
    fn standalone_custom_profile_needs_all_base_fields() {
        let (lib, _) = lib_with(&[("mine", "name = \"Mine\"\ndry_mass_kg = 1000.0\n")]);
        assert!(lib.resolve("mine", &SpacecraftParams::default()).is_err());
        let (lib, _) = lib_with(&[(
            "full",
            "name = \"Full\"\ndry_mass_kg = 1000.0\npropellant_mass_kg = 1500.0\npayload_mass_kg = 0.0\nmax_thrust_n = 20000.0\nisp_s = 320.0\n",
        )]);
        let r = lib.resolve("full", &SpacecraftParams::default()).unwrap();
        assert_eq!(r.status_label, "Custom");
    }

    #[test]
    fn save_and_reload_custom_profile() {
        let (lib, dir) = lib_with(&[]);
        let mut o = SpacecraftParams::default();
        o.set_from_str("payload_mass_kg=750").unwrap();
        let path = lib.save_custom("my_ascent", "apollo_lm_ascent", &o).unwrap();
        assert!(path.starts_with(&dir));
        let lib2 = SpacecraftLibrary::load(&dir).unwrap();
        let r = lib2.resolve("my_ascent", &SpacecraftParams::default()).unwrap();
        assert_eq!(r.config.payload_mass_kg, 750.0);
        assert_eq!(r.status_label, "Custom (на базе Apollo LM Ascent Stage)");
        // Только изменённые поля в файле.
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.contains("payload_mass_kg") && !text.contains("isp_s"));
        // Шаблон перезаписать нельзя.
        assert!(lib2.save_custom("apollo_lm_ascent", "apollo_lm_ascent", &o).is_err());
        assert!(lib2.save_custom("Bad Name", "apollo_lm_ascent", &o).is_err());
    }

    #[test]
    fn rejects_cycles_unknown_fields_and_duplicate_ids() {
        let dir = std::env::temp_dir().join(format!("lwp-sc-bad-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("a.toml"), "based_on = \"b\"\n").unwrap();
        std::fs::write(dir.join("b.toml"), "based_on = \"a\"\n").unwrap();
        let lib = SpacecraftLibrary::load(&dir).unwrap();
        assert!(lib.resolve("a", &SpacecraftParams::default()).is_err());

        std::fs::write(dir.join("c.toml"), "isp = 300\n").unwrap();
        assert!(SpacecraftLibrary::load(&dir).is_err());
        std::fs::remove_file(dir.join("c.toml")).unwrap();

        std::fs::write(dir.join("apollo_lm_ascent.toml"), "isp_s = 300\n").unwrap();
        assert!(SpacecraftLibrary::load(&dir).is_err());
    }

    #[test]
    fn bundled_custom_example_resolves() {
        let lib = SpacecraftLibrary::load("data/spacecraft/custom").unwrap();
        let r = lib.resolve("apollo_lm_heavy_payload", &SpacecraftParams::default()).unwrap();
        assert_eq!(r.status_label, "Custom (на базе Apollo LM Descent Stage)");
        assert!(!r.has_errors());
    }
}
