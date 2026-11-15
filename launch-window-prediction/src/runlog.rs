//! Журнал информации о полёте: один JSON-документ, удобный для чтения человеком.
//!
//! Каждый запуск пишет свой файл `logs/flight-<время>.json`. Рядом с тем же временем в имени
//! пишется журнал траектории для 3D (`trajectory-<время>.json`, модуль [`crate::trajlog`]).
//!
//! Устройство документа (сверху вниз):
//! - `format`, `about` — формат и короткое описание разделов;
//! - разделы-сводки, которые собирает программа по ходу расчёта ([`RunLog::section`],
//!   [`RunLog::merge`]): `run`, `summary`, `launch_site`, `landing_site`, `files`;
//! - `log` — все события расчёта по порядку: `{"seq", "wall_time_utc", "event", "data"}`;
//!   в `data` попадают все входные данные, параметры и промежуточные результаты.
//!
//! Файл форматируется с отступами; короткие массивы чисел (векторы) пишутся в одну строку,
//! длинные — по нескольку чисел в строке. Дробные числа округлены до 6 знаков после запятой
//! (очень малые — до 6 значащих цифр), время — до миллисекунд без лишних нулей.
//! Полная точность положений и скоростей — в журнале траектории.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use anise::math::Vector3;
use anise::time::{Duration, Epoch, Unit};
use serde::{Serialize, Serializer};
use serde_json::{Map, Value};

use crate::error::{Error, Result};

/// Идентификатор формата журнала информации о полёте.
pub const FORMAT: &str = "lwp-flight/2";

/// Разделы-сводки в порядке вывода (перед журналом событий).
const SECTIONS: [&str; 5] = ["run", "summary", "launch_site", "landing_site", "files"];

/// Не чаще, чем раз в столько секунд, файл переписывается после очередного события
/// (разделы и завершение пишутся сразу).
const WRITE_INTERVAL_S: f64 = 0.5;

/// Ширина строки при форматировании.
const WIDTH: usize = 100;

/// Журнал одного запуска.
pub struct RunLog {
    path: PathBuf,
    state: Mutex<State>,
    seq: AtomicU64,
}

struct State {
    sections: Vec<(String, Value)>,
    events: Vec<Value>,
    last_write: Option<Instant>,
    dirty: bool,
}

#[derive(Serialize)]
struct Line<'a, T: Serialize> {
    seq: u64,
    wall_time_utc: String,
    event: &'a str,
    data: T,
}

impl RunLog {
    /// Создаёт новый файл журнала в папке `dir` (папка создаётся при необходимости).
    pub fn create(dir: impl AsRef<Path>) -> Result<Self> {
        let dir = dir.as_ref();
        std::fs::create_dir_all(dir).map_err(|source| Error::Io { path: dir.to_path_buf(), source })?;
        let stamp = wall_now()
            .round(1 * Unit::Second)
            .to_string()
            .replace([':', ' '], "")
            .replace("UTC", "Z");
        let mut path = dir.join(format!("flight-{stamp}.json"));
        let mut n = 1;
        while path.exists() || dir.join(format!("trajectory-{}", Self::suffix(&path))).exists() {
            n += 1;
            path = dir.join(format!("flight-{stamp}-{n}.json"));
        }
        let log = Self {
            path,
            state: Mutex::new(State {
                sections: SECTIONS.iter().map(|s| (s.to_string(), Value::Null)).collect(),
                events: Vec::new(),
                last_write: None,
                dirty: true,
            }),
            seq: AtomicU64::new(0),
        };
        // Файл появляется сразу: так имя занято и второй запуск в ту же секунду возьмёт другое.
        let mut st = log.lock();
        log.write_locked(&mut st).map_err(|source| Error::Io { path: log.path.clone(), source })?;
        drop(st);
        Ok(log)
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Часть имени после `flight-` с расширением `.json`.
    fn suffix(path: &Path) -> String {
        let name = path.file_stem().and_then(|s| s.to_str()).unwrap_or("flight");
        format!("{}.json", name.trim_start_matches("flight-"))
    }

    /// Путь журнала траектории с тем же временем в имени.
    pub fn trajectory_path(&self) -> PathBuf {
        self.path.with_file_name(format!("trajectory-{}", Self::suffix(&self.path)))
    }

    fn lock(&self) -> std::sync::MutexGuard<'_, State> {
        self.state.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// Записывает событие. Ошибка записи не прерывает расчёт: она уходит в предупреждение.
    pub fn event<T: Serialize>(&self, event: &str, data: T) {
        let line = Line {
            seq: self.seq.fetch_add(1, Ordering::Relaxed),
            wall_time_utc: iso(&wall_now()),
            event,
            data,
        };
        let value = match serde_json::to_value(&line) {
            Ok(v) => v,
            Err(e) => {
                tracing::warn!(error = %e, event, "failed to serialize run log event");
                return;
            }
        };
        let mut st = self.lock();
        st.events.push(value);
        st.dirty = true;
        let due = st.last_write.is_none_or(|t| t.elapsed().as_secs_f64() >= WRITE_INTERVAL_S);
        if due {
            self.write_reporting(&mut st);
        }
    }

    /// Задаёт раздел-сводку целиком (заменяет прежнее содержимое).
    pub fn section<T: Serialize>(&self, name: &str, data: T) {
        let Some(value) = self.to_value(name, data) else { return };
        let mut st = self.lock();
        match st.sections.iter_mut().find(|(n, _)| n == name) {
            Some((_, v)) => *v = value,
            None => st.sections.push((name.to_string(), value)),
        }
        st.dirty = true;
        self.write_reporting(&mut st);
    }

    /// Добавляет поля в раздел-сводку (поля с теми же именами заменяются).
    pub fn merge<T: Serialize>(&self, name: &str, data: T) {
        let Some(value) = self.to_value(name, data) else { return };
        let mut st = self.lock();
        let slot = match st.sections.iter().position(|(n, _)| n == name) {
            Some(i) => &mut st.sections[i].1,
            None => {
                st.sections.push((name.to_string(), Value::Null));
                &mut st.sections.last_mut().expect("just pushed").1
            }
        };
        match (slot, value) {
            (Value::Object(old), Value::Object(new)) => old.extend(new),
            (slot, value) => *slot = value,
        }
        st.dirty = true;
        self.write_reporting(&mut st);
    }

    /// Дописывает файл, если после последней записи что-то изменилось.
    pub fn flush(&self) {
        let mut st = self.lock();
        if st.dirty {
            self.write_reporting(&mut st);
        }
    }

    fn to_value<T: Serialize>(&self, name: &str, data: T) -> Option<Value> {
        serde_json::to_value(data)
            .map_err(|e| tracing::warn!(error = %e, section = name, "failed to serialize run log section"))
            .ok()
    }

    fn write_reporting(&self, st: &mut State) {
        if let Err(e) = self.write_locked(st) {
            tracing::warn!(error = %e, path = %self.path.display(), "failed to write run log");
        }
    }

    fn write_locked(&self, st: &mut State) -> std::io::Result<()> {
        let mut doc = Map::new();
        doc.insert("format".into(), Value::from(FORMAT));
        doc.insert("about".into(), about());
        for (name, value) in &st.sections {
            if !value.is_null() {
                doc.insert(name.clone(), value.clone());
            }
        }
        doc.insert("log".into(), Value::Array(st.events.clone()));
        let mut text = String::new();
        render(&Value::Object(doc), 0, &mut text);
        text.push('\n');
        std::fs::write(&self.path, text)?;
        st.last_write = Some(Instant::now());
        st.dirty = false;
        Ok(())
    }
}

impl Drop for RunLog {
    fn drop(&mut self) {
        self.flush();
    }
}

fn about() -> Value {
    serde_json::json!({
        "run": "program, version, command line, when it ran and how it ended",
        "summary": "the answer in short: dates of every flight stage, flight times, delta-v, alternatives",
        "launch_site": "everything about the spaceport and the launch from it",
        "landing_site": "everything about the landing point on the Moon: terrain, horizon, the landing itself, Sun and Earth",
        "files": "other files written by this run",
        "log": "every step of the calculation in order, with all inputs and intermediate results",
        "units": "units are in the key names (_km, _m_s, _deg, ...); times are UTC; t_plus is time since launch",
    })
}

fn wall_now() -> Epoch {
    Epoch::now().unwrap_or_else(|_| Epoch::from_gregorian_utc_at_midnight(1970, 1, 1))
}

/// Момент в виде строки UTC с миллисекундами.
pub fn iso(e: &Epoch) -> String {
    e.round(1 * Unit::Millisecond).to_string()
}

/// Момент для человека: `2026-12-29T05:24:12.229 UTC` (без лишних нулей в дробной части).
pub fn iso_short(e: &Epoch) -> String {
    trim_time(&iso(e)).unwrap_or_else(|| iso(e))
}

/// Длительность для человека: `2 d 14 h 32 min`, `1 h 17 min 30 s`, `45.2 s`.
pub fn human_duration(seconds: f64) -> String {
    if !seconds.is_finite() {
        return "?".into();
    }
    let sign = if seconds < 0.0 { "-" } else { "" };
    let s = seconds.abs();
    if s < 60.0 {
        return format!("{sign}{} s", (s * 10.0).round() / 10.0);
    }
    let total = s.round() as u64;
    let (d, h, m, sec) = (total / 86_400, total / 3600 % 24, total / 60 % 60, total % 60);
    let mut parts = Vec::new();
    if d > 0 {
        parts.push(format!("{d} d"));
    }
    if h > 0 {
        parts.push(format!("{h} h"));
    }
    if m > 0 {
        parts.push(format!("{m} min"));
    }
    // Секунды — только у коротких промежутков.
    if sec > 0 && d == 0 {
        parts.push(format!("{sec} s"));
    }
    format!("{sign}{}", parts.join(" "))
}

/// `2026-12-29T05:24:12.229000000 UTC` → `2026-12-29T05:24:12.229 UTC`.
fn trim_time(s: &str) -> Option<String> {
    let body = s.strip_suffix(" UTC")?;
    let (head, frac) = body.split_once('.')?;
    if !head.contains('T') || frac.is_empty() || !frac.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let frac = frac.trim_end_matches('0');
    Some(if frac.is_empty() { format!("{head} UTC") } else { format!("{head}.{frac} UTC") })
}

// ---- Форматирование ----

fn number(n: &serde_json::Number) -> String {
    if n.is_i64() || n.is_u64() {
        return n.to_string();
    }
    let Some(x) = n.as_f64() else { return n.to_string() };
    let ax = x.abs();
    let r = if x == 0.0 || ax >= 1e12 {
        x
    } else if ax >= 1e-3 {
        (x * 1e6).round() / 1e6
    } else {
        format!("{x:.5e}").parse().unwrap_or(x)
    };
    serde_json::to_string(&r).unwrap_or_else(|_| n.to_string())
}

fn scalar(v: &Value) -> String {
    match v {
        Value::Number(n) => number(n),
        Value::String(s) => {
            let s = trim_time(s).map_or_else(|| s.clone(), |t| t);
            serde_json::to_string(&s).unwrap_or_default()
        }
        other => other.to_string(),
    }
}

fn is_scalar(v: &Value) -> bool {
    !matches!(v, Value::Array(_) | Value::Object(_))
}

fn pad(out: &mut String, indent: usize) {
    out.extend(std::iter::repeat_n(' ', indent));
}

fn render(v: &Value, indent: usize, out: &mut String) {
    match v {
        Value::Array(items) if items.is_empty() => out.push_str("[]"),
        Value::Object(map) if map.is_empty() => out.push_str("{}"),
        Value::Array(items) if items.iter().all(is_scalar) => {
            let parts: Vec<String> = items.iter().map(scalar).collect();
            let inline = parts.join(", ");
            if indent + inline.len() + 2 <= WIDTH {
                let _ = write!(out, "[{inline}]");
                return;
            }
            // Длинный ряд чисел — по нескольку в строке.
            out.push_str("[\n");
            let mut line = String::new();
            for (i, p) in parts.iter().enumerate() {
                let piece = if i + 1 < parts.len() { format!("{p},") } else { p.clone() };
                if !line.is_empty() && indent + 2 + line.len() + 1 + piece.len() > WIDTH {
                    pad(out, indent + 2);
                    out.push_str(&line);
                    out.push('\n');
                    line.clear();
                }
                if !line.is_empty() {
                    line.push(' ');
                }
                line.push_str(&piece);
            }
            pad(out, indent + 2);
            out.push_str(&line);
            out.push('\n');
            pad(out, indent);
            out.push(']');
        }
        Value::Array(items) => {
            out.push_str("[\n");
            for (i, item) in items.iter().enumerate() {
                pad(out, indent + 2);
                render(item, indent + 2, out);
                out.push_str(if i + 1 < items.len() { ",\n" } else { "\n" });
            }
            pad(out, indent);
            out.push(']');
        }
        Value::Object(map) => {
            // Маленький объект из простых значений — в одну строку.
            if map.len() <= 4 && map.values().all(is_scalar) {
                let inline: Vec<String> = map
                    .iter()
                    .map(|(k, v)| format!("{}: {}", serde_json::to_string(k).unwrap_or_default(), scalar(v)))
                    .collect();
                let inline = format!("{{ {} }}", inline.join(", "));
                if indent + inline.len() <= WIDTH - 20 {
                    out.push_str(&inline);
                    return;
                }
            }
            out.push_str("{\n");
            for (i, (k, v)) in map.iter().enumerate() {
                pad(out, indent + 2);
                let _ = write!(out, "{}: ", serde_json::to_string(k).unwrap_or_default());
                render(v, indent + 2, out);
                out.push_str(if i + 1 < map.len() { ",\n" } else { "\n" });
            }
            pad(out, indent);
            out.push('}');
        }
        other => out.push_str(&scalar(other)),
    }
}

// ---- Помощники для `#[serde(serialize_with = ...)]` ----

pub fn ser_epoch<S: Serializer>(e: &Epoch, s: S) -> std::result::Result<S::Ok, S::Error> {
    s.serialize_str(&iso(e))
}

pub fn ser_opt_epoch<S: Serializer>(e: &Option<Epoch>, s: S) -> std::result::Result<S::Ok, S::Error> {
    match e {
        Some(e) => s.serialize_some(&iso(e)),
        None => s.serialize_none(),
    }
}

/// Длительность в секундах.
pub fn ser_duration_s<S: Serializer>(d: &Duration, s: S) -> std::result::Result<S::Ok, S::Error> {
    s.serialize_f64(d.to_seconds())
}

pub fn ser_vec3<S: Serializer>(v: &Vector3, s: S) -> std::result::Result<S::Ok, S::Error> {
    [v.x, v.y, v.z].serialize(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_a_readable_json_document() {
        let dir = std::env::temp_dir().join(format!("lwp-runlog-{}", std::process::id()));
        let log = RunLog::create(&dir).unwrap();
        log.event("first", serde_json::json!({ "a": 1 }));
        #[derive(Serialize)]
        struct S {
            #[serde(serialize_with = "ser_epoch")]
            t: Epoch,
            #[serde(serialize_with = "ser_vec3")]
            v: Vector3,
            #[serde(serialize_with = "ser_duration_s")]
            d: Duration,
            x: f64,
            tiny: f64,
            long: Vec<f64>,
        }
        log.event(
            "second",
            S {
                t: Epoch::from_gregorian_utc_hms(2026, 12, 1, 7, 25, 33) + 229 * Unit::Millisecond,
                v: Vector3::new(1.0, 2.0, 3.0),
                d: 90 * Unit::Second,
                x: 4.179197862343359,
                tiny: 6.363123678062213e-05,
                long: (0..100).map(f64::from).collect(),
            },
        );
        log.section("summary", serde_json::json!({ "status": "ok" }));
        log.merge("landing_site", serde_json::json!({ "lat_deg": 5.0 }));
        log.merge("landing_site", serde_json::json!({ "lon_deg": -10.0 }));
        log.flush();
        let text = std::fs::read_to_string(log.path()).unwrap();
        assert!(log.path().extension().is_some_and(|e| e == "json"));
        assert!(text.lines().count() > 10, "pretty-printed:\n{text}");
        assert!(text.contains("[1.0, 2.0, 3.0]"), "vectors on one line:\n{text}");
        let doc: Value = serde_json::from_str(&text).unwrap();
        assert_eq!(doc["format"], FORMAT);
        // Разделы идут перед журналом.
        let keys: Vec<&str> = doc.as_object().unwrap().keys().map(String::as_str).collect();
        assert_eq!(keys, ["format", "about", "summary", "landing_site", "log"]);
        assert_eq!(doc["landing_site"], serde_json::json!({ "lat_deg": 5.0, "lon_deg": -10.0 }));
        let log_lines = doc["log"].as_array().unwrap();
        assert_eq!(log_lines.len(), 2);
        assert_eq!(log_lines[0]["seq"], 0);
        assert_eq!(log_lines[0]["event"], "first");
        let data = &log_lines[1]["data"];
        assert_eq!(data["t"], "2026-12-01T07:25:33.229 UTC");
        assert_eq!(data["v"], serde_json::json!([1.0, 2.0, 3.0]));
        assert_eq!(data["d"], 90.0);
        assert_eq!(data["x"], 4.179198);
        assert_eq!(data["tiny"], 6.36312e-5);
        assert_eq!(data["long"].as_array().unwrap().len(), 100);
        // Второй журнал в той же секунде не затирает первый.
        let other = RunLog::create(&dir).unwrap();
        assert_ne!(other.path(), log.path());
    }

    #[test]
    fn durations_and_times_for_people() {
        assert_eq!(human_duration(45.21), "45.2 s");
        assert_eq!(human_duration(4650.0), "1 h 17 min 30 s");
        assert_eq!(human_duration(225_257.8), "2 d 14 h 34 min");
        assert_eq!(trim_time("2026-12-31T21:16:00 UTC"), None);
        assert_eq!(trim_time("2026-12-31T21:16:00.500000000 UTC").as_deref(), Some("2026-12-31T21:16:00.5 UTC"));
        assert_eq!(trim_time("2026-12-31T21:16:00.000000000 UTC").as_deref(), Some("2026-12-31T21:16:00 UTC"));
    }
}
