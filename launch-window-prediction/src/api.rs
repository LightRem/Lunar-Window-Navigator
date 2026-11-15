//! Машинный интерфейс (API) для других программ: `launch-window-prediction api <команда>`.
//!
//! Обмен — JSON: запрос читается из стандартного ввода (или из файла `--input`), ответ пишется в
//! стандартный вывод одним JSON-объектом. Журнал работы программы идёт в поток ошибок, поэтому
//! стандартный вывод можно сразу разбирать как JSON. Код выхода: 0 — успех, 1 — внутренняя
//! ошибка, 2 — неверный запрос, 3 — задача не найдена, 4 — результат ещё не готов.
//!
//! Расчёт занимает от секунды до минуты, поэтому он идёт задачей в фоне:
//! `submit` создаёт задачу и сразу возвращает её `job_id`, фоновый процесс считает,
//! `status` показывает состояние, `result` отдаёт данные. У каждой задачи своя папка
//! `<jobs-dir>/<job_id>/`: запрос, состояние, журналы расчёта и текстовый вывод.
//!
//! Подробное описание с примерами — `api help` и файл `api.md`.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Stdio;

use clap::Parser;
use launch_window_prediction::launch_vehicle::{LaunchVehicleLibrary, LaunchVehicleParams};
use launch_window_prediction::spacecraft::{SpacecraftLibrary, SpacecraftParams};
use launch_window_prediction::spaceport::SpaceportCatalog;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value, json};

use super::{Cli, execute, parse_date};

/// Версия интерфейса: меняется, если меняются поля запроса или ответа несовместимо.
pub const API_VERSION: &str = "lwp-api/1";

const DEFAULT_JOBS_DIR: &str = "jobs";
const DEFAULT_KERNELS_DIR: &str = "data/kernels";
const DEFAULT_TERRAIN_DIR: &str = "data/terrain";
const DEFAULT_SPACEPORTS: &str = "data/spaceports.toml";
const DEFAULT_ROCKET_DIR: &str = "data/launch_vehicles/custom";
const DEFAULT_SPACECRAFT_DIR: &str = "data/spacecraft/custom";

// ---------------------------------------------------------------------------
// Ошибки и ответ
// ---------------------------------------------------------------------------

#[derive(Debug)]
struct ApiError {
    code: &'static str,
    message: String,
    exit: i32,
}

impl ApiError {
    fn bad(message: impl Into<String>) -> Self {
        Self { code: "bad_request", message: message.into(), exit: 2 }
    }
    fn not_found(message: impl Into<String>) -> Self {
        Self { code: "not_found", message: message.into(), exit: 3 }
    }
    fn not_ready(message: impl Into<String>) -> Self {
        Self { code: "not_ready", message: message.into(), exit: 4 }
    }
    fn internal(message: impl Into<String>) -> Self {
        Self { code: "internal", message: message.into(), exit: 1 }
    }
}

impl From<std::io::Error> for ApiError {
    fn from(e: std::io::Error) -> Self {
        Self::internal(e.to_string())
    }
}

impl From<serde_json::Error> for ApiError {
    fn from(e: serde_json::Error) -> Self {
        Self::internal(e.to_string())
    }
}

type ApiResult<T> = std::result::Result<T, ApiError>;

/// Точка входа: разбирает аргументы, выполняет команду, печатает JSON. Возвращает код выхода.
pub async fn main(args: &[String]) -> i32 {
    let opts = match Opts::parse(args) {
        Ok(o) => o,
        Err(e) => return print_error(&e, false),
    };
    // Фоновый процесс задачи ничего не печатает в формате API: его вывод — текст расчёта.
    if opts.command.as_deref() == Some("run") {
        return match worker(&opts).await {
            Ok(()) => 0,
            Err(e) => {
                eprintln!("{}: {}", e.code, e.message);
                e.exit
            }
        };
    }
    match dispatch(&opts).await {
        Ok(body) => {
            let mut out = Map::new();
            out.insert("ok".into(), Value::Bool(true));
            out.insert("api".into(), Value::from(API_VERSION));
            if let Value::Object(m) = body {
                out.extend(m);
            }
            print_json(&Value::Object(out), opts.pretty);
            0
        }
        Err(e) => print_error(&e, opts.pretty),
    }
}

fn print_json(v: &Value, pretty: bool) {
    let text = if pretty { serde_json::to_string_pretty(v) } else { serde_json::to_string(v) };
    println!("{}", text.unwrap_or_else(|e| format!("{{\"ok\":false,\"error\":{{\"code\":\"internal\",\"message\":\"{e}\"}}}}")));
}

fn print_error(e: &ApiError, pretty: bool) -> i32 {
    print_json(&json!({ "ok": false, "api": API_VERSION, "error": { "code": e.code, "message": e.message } }), pretty);
    e.exit
}

// ---------------------------------------------------------------------------
// Аргументы
// ---------------------------------------------------------------------------

struct Opts {
    jobs_dir: PathBuf,
    input: Option<PathBuf>,
    part: String,
    wait: bool,
    pretty: bool,
    command: Option<String>,
    positional: Vec<String>,
}

impl Opts {
    fn parse(args: &[String]) -> ApiResult<Self> {
        let mut o = Self {
            jobs_dir: PathBuf::from(DEFAULT_JOBS_DIR),
            input: None,
            part: "summary".into(),
            wait: false,
            pretty: false,
            command: None,
            positional: Vec::new(),
        };
        let mut it = args.iter();
        while let Some(a) = it.next() {
            let mut value = |name: &str| it.next().cloned().ok_or_else(|| ApiError::bad(format!("{name} needs a value")));
            match a.as_str() {
                "--jobs-dir" => o.jobs_dir = PathBuf::from(value("--jobs-dir")?),
                "--input" => o.input = Some(PathBuf::from(value("--input")?)),
                "--part" => o.part = value("--part")?,
                "--wait" => o.wait = true,
                "--pretty" => o.pretty = true,
                s if s.starts_with("--") => return Err(ApiError::bad(format!("unknown option '{s}' (see 'api help')"))),
                s if o.command.is_none() => o.command = Some(s.to_string()),
                s => o.positional.push(s.to_string()),
            }
        }
        Ok(o)
    }

    fn job_id(&self) -> ApiResult<&str> {
        let id = self.positional.first().ok_or_else(|| ApiError::bad("job id is required"))?;
        // Только безопасные символы: id — имя папки.
        if id.is_empty() || !id.chars().all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_') {
            return Err(ApiError::bad(format!("invalid job id '{id}'")));
        }
        Ok(id)
    }

    fn job_dir(&self) -> ApiResult<PathBuf> {
        let dir = self.jobs_dir.join(self.job_id()?);
        if !dir.join(STATUS_FILE).exists() {
            return Err(ApiError::not_found(format!("job '{}' not found in {}", self.job_id()?, self.jobs_dir.display())));
        }
        Ok(dir)
    }
}

async fn dispatch(o: &Opts) -> ApiResult<Value> {
    match o.command.as_deref() {
        None | Some("help") => Ok(help()),
        Some("catalog") => catalog().await,
        Some("submit") => submit(o).await,
        Some("status") => read_status(&o.job_dir()?).map(|s| json!({ "job": s })),
        Some("result") => result(o),
        Some("list") => list(o),
        Some("delete") => delete(o),
        Some(other) => Err(ApiError::bad(format!("unknown command '{other}' (see 'api help')"))),
    }
}

// ---------------------------------------------------------------------------
// Запрос на расчёт
// ---------------------------------------------------------------------------

/// Запрос на расчёт (тело команды `submit`). Все поля необязательны; нужна хотя бы дата
/// (`launch_date` или `arrival_date`) или точка посадки.
#[derive(Debug, Clone, Default, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Request {
    /// Сутки старта (UTC), "ГГГГ-ММ-ДД": найти лучший старт в эти сутки.
    pub launch_date: Option<String>,
    /// Сутки прилёта (UTC), "ГГГГ-ММ-ДД". С точкой посадки — сутки посадки в неё.
    pub arrival_date: Option<String>,
    /// Космодром (id из каталога). Без него — все космодромы (или космодромы ракеты).
    pub spaceport: Option<String>,
    /// Ракета-носитель (id из каталога).
    pub rocket: Option<String>,
    /// Переопределения параметров ракеты: { "поле": значение }.
    #[serde(default)]
    pub rocket_overrides: BTreeMap<String, Value>,
    /// Лунный аппарат (id из каталога).
    pub spacecraft: Option<String>,
    /// Переопределения параметров аппарата: { "поле": значение }.
    #[serde(default)]
    pub spacecraft_overrides: BTreeMap<String, Value>,
    /// Точка посадки на Луне.
    pub landing_site: Option<LandingSite>,
    /// Высота поверхности над опорной сферой, м (только без карт рельефа).
    pub height_m: Option<f64>,
    /// Ряд видимости Солнца и Земли из точки посадки.
    pub visibility: Option<Visibility>,
    /// Папка ядер эфемерид; по умолчанию "data/kernels".
    pub kernels_dir: Option<String>,
    /// Скачать ядра из интернета вместо `kernels_dir`.
    #[serde(default)]
    pub download_kernels: bool,
    /// Папка карт рельефа; по умолчанию "data/terrain", если она есть. "" — без рельефа.
    pub terrain_dir: Option<String>,
    /// Свои файлы параметров (по умолчанию — из data/).
    pub transfer_config: Option<String>,
    pub landing_config: Option<String>,
    pub visibility_config: Option<String>,
    /// Не использовать кэш моделирования подъёма.
    #[serde(default)]
    pub no_cache: bool,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct LandingSite {
    pub lat_deg: f64,
    pub lon_deg: f64,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Visibility {
    /// Начало ряда, например "2027-02-14T00:00:00 UTC"; по умолчанию — момент посадки.
    pub start: Option<String>,
    /// Длительность, часы (0 — один момент).
    pub hours: Option<f64>,
    /// Шаг, минуты.
    pub step_min: Option<f64>,
}

fn override_value(v: &Value) -> String {
    match v {
        Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

impl Request {
    /// Аргументы командной строки обычного режима; журналы пишутся в `log_dir`.
    fn to_args(&self, log_dir: &Path) -> Vec<String> {
        let mut a: Vec<String> = vec![env!("CARGO_PKG_NAME").into(), "--log-dir".into(), log_dir.display().to_string()];
        let mut push = |k: &str, v: String| {
            a.push(k.into());
            a.push(v);
        };
        if let Some(d) = &self.launch_date {
            push("--launch-date", d.clone());
        }
        if let Some(d) = &self.arrival_date {
            push("--arrival-date", d.clone());
        }
        if let Some(p) = &self.spaceport {
            push("--spaceport", p.clone());
        }
        if let Some(r) = &self.rocket {
            push("--rocket", r.clone());
        }
        for (k, v) in &self.rocket_overrides {
            push("--rk", format!("{k}={}", override_value(v)));
        }
        if let Some(s) = &self.spacecraft {
            push("--spacecraft", s.clone());
        }
        for (k, v) in &self.spacecraft_overrides {
            push("--sc", format!("{k}={}", override_value(v)));
        }
        if let Some(s) = &self.landing_site {
            push("--lat", s.lat_deg.to_string());
            push("--lon", s.lon_deg.to_string());
        }
        if let Some(h) = self.height_m {
            push("--height-m", h.to_string());
        }
        if let Some(v) = &self.visibility {
            if let Some(s) = &v.start {
                push("--start", s.clone());
            }
            if let Some(h) = v.hours {
                push("--hours", h.to_string());
            }
            if let Some(m) = v.step_min {
                push("--step-min", m.to_string());
            }
        }
        if self.download_kernels {
            a.push("--download-kernels".into());
        } else {
            push("--kernels-dir", self.kernels_dir.clone().unwrap_or_else(|| DEFAULT_KERNELS_DIR.into()));
        }
        let mut push = |k: &str, v: String| {
            a.push(k.into());
            a.push(v);
        };
        match self.terrain_dir.as_deref() {
            Some("") => {}
            Some(d) => push("--terrain-dir", d.into()),
            None if Path::new(DEFAULT_TERRAIN_DIR).is_dir() => push("--terrain-dir", DEFAULT_TERRAIN_DIR.into()),
            None => {}
        }
        for (flag, v) in [
            ("--transfer-config", &self.transfer_config),
            ("--landing-config", &self.landing_config),
            ("--visibility-config", &self.visibility_config),
        ] {
            if let Some(v) = v {
                push(flag, v.clone());
            }
        }
        if self.no_cache {
            a.push("--no-cache".into());
        }
        a
    }

    /// Проверка до запуска: формат, допустимые сочетания, id и переопределения по каталогам.
    async fn validate(&self) -> ApiResult<Cli> {
        if self.launch_date.is_none() && self.arrival_date.is_none() && self.landing_site.is_none() {
            return Err(ApiError::bad("nothing to calculate: give launch_date or arrival_date, and/or landing_site"));
        }
        for d in [&self.launch_date, &self.arrival_date].into_iter().flatten() {
            parse_date(d).map_err(|e| ApiError::bad(e.to_string()))?;
        }
        if let Some(s) = &self.landing_site {
            if !(-90.0..=90.0).contains(&s.lat_deg) || !(-180.0..=360.0).contains(&s.lon_deg) {
                return Err(ApiError::bad("landing_site: lat_deg must be in [-90, 90], lon_deg in [-180, 360]"));
            }
        }
        let cli = Cli::try_parse_from(self.to_args(Path::new("."))).map_err(|e| {
            ApiError::bad(e.to_string().lines().next().unwrap_or("invalid request").trim_start_matches("error: ").to_string())
        })?;
        if let Some(dir) = &cli.kernels_dir {
            if !dir.is_dir() {
                return Err(ApiError::bad(format!("kernels_dir '{}' does not exist", dir.display())));
            }
        }
        if let Some(dir) = &cli.terrain_dir {
            if !dir.is_dir() {
                return Err(ApiError::bad(format!("terrain_dir '{}' does not exist", dir.display())));
            }
        }
        let catalog = SpaceportCatalog::load(&cli.spaceports).await.map_err(|e| ApiError::internal(e.to_string()))?;
        if let Some(p) = &self.spaceport {
            if catalog.get(p).is_none() {
                let ids: Vec<&str> = catalog.iter().map(|s| s.id.as_str()).collect();
                return Err(ApiError::bad(format!("unknown spaceport '{p}' (known: {})", ids.join(", "))));
            }
        }
        let (sd, rd) = (cli.spacecraft_dir.clone(), cli.rocket_dir.clone());
        let (library, rockets) = tokio::task::spawn_blocking(move || {
            Ok::<_, launch_window_prediction::Error>((SpacecraftLibrary::load(sd)?, LaunchVehicleLibrary::load(rd)?))
        })
        .await
        .map_err(|e| ApiError::internal(e.to_string()))?
        .map_err(|e| ApiError::internal(e.to_string()))?;
        if let Some(id) = &self.spacecraft {
            let mut ov = SpacecraftParams::default();
            for a in &cli.spacecraft_overrides {
                ov.set_from_str(a).map_err(|e| ApiError::bad(format!("spacecraft_overrides: {e}")))?;
            }
            library.resolve(id, &ov).map_err(|e| ApiError::bad(format!("spacecraft: {e}")))?;
        }
        if let Some(id) = &self.rocket {
            let mut ov = LaunchVehicleParams::default();
            for a in &cli.rocket_overrides {
                ov.set_from_str(a).map_err(|e| ApiError::bad(format!("rocket_overrides: {e}")))?;
            }
            let r = rockets.resolve(id, &ov).map_err(|e| ApiError::bad(format!("rocket: {e}")))?;
            if let Some(p) = &self.spaceport {
                if !r.config.any_spaceport && !r.config.spaceports.contains(p) {
                    return Err(ApiError::bad(format!(
                        "rocket '{id}' does not launch from '{p}' (its spaceports: {})",
                        r.config.spaceports.join(", ")
                    )));
                }
            }
        }
        Ok(cli)
    }
}

// ---------------------------------------------------------------------------
// Задачи
// ---------------------------------------------------------------------------

const REQUEST_FILE: &str = "request.json";
const STATUS_FILE: &str = "status.json";
const OUTPUT_FILE: &str = "output.txt";
const WORKER_LOG: &str = "worker.log";

/// Состояние задачи (`status.json` в папке задачи).
#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct JobStatus {
    pub job_id: String,
    /// queued — создана, фоновый процесс запускается; running — считается; done — готово
    /// (смотрите `outcome`: расчёт мог и остановиться, например на непригодной площадке);
    /// failed — ошибка (см. `error`).
    pub status: String,
    pub created_utc: String,
    pub started_utc: Option<String>,
    pub finished_utc: Option<String>,
    pub elapsed_s: Option<f64>,
    /// Итог расчёта: "done" или "stopped: …" / "done: …" с причиной.
    pub outcome: Option<String>,
    pub error: Option<String>,
    /// Процесс, который считает задачу.
    pub pid: Option<u32>,
    /// Абсолютные пути файлов результата (после завершения).
    pub flight_log: Option<String>,
    pub trajectory_log: Option<String>,
    pub job_dir: String,
}

fn now_utc() -> String {
    launch_window_prediction::runlog::iso_short(&anise::time::Epoch::now().unwrap_or_default())
}

fn new_job_id() -> String {
    let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap_or_default();
    let stamp = anise::time::Epoch::now()
        .unwrap_or_default()
        .round(anise::time::Duration::from_seconds(1.0))
        .to_string();
    // "2026-11-15T14:05:12 UTC" -> "20261115-140512"
    let digits: String = stamp.chars().filter(|c| c.is_ascii_digit()).take(14).collect();
    let (date, time) = digits.split_at(digits.len().min(8));
    let salt = (now.subsec_nanos() ^ std::process::id().rotate_left(16)) & 0xffff;
    format!("{date}-{time}-{salt:04x}")
}

fn write_status(dir: &Path, s: &JobStatus) -> ApiResult<()> {
    let tmp = dir.join(format!("{STATUS_FILE}.tmp"));
    std::fs::write(&tmp, serde_json::to_string_pretty(s)?)?;
    std::fs::rename(&tmp, dir.join(STATUS_FILE))?;
    Ok(())
}

fn read_status(dir: &Path) -> ApiResult<JobStatus> {
    let text = std::fs::read_to_string(dir.join(STATUS_FILE))?;
    Ok(serde_json::from_str(&text)?)
}

fn absolute(p: &Path) -> String {
    std::path::absolute(p).unwrap_or_else(|_| p.to_path_buf()).display().to_string()
}

fn read_request(o: &Opts) -> ApiResult<Request> {
    let text = match &o.input {
        Some(p) => std::fs::read_to_string(p).map_err(|e| ApiError::bad(format!("cannot read {}: {e}", p.display())))?,
        None => std::io::read_to_string(std::io::stdin())?,
    };
    if text.trim().is_empty() {
        return Err(ApiError::bad("empty request: send a JSON object on standard input or use --input FILE"));
    }
    serde_json::from_str(&text).map_err(|e| ApiError::bad(format!("request JSON: {e}")))
}

async fn submit(o: &Opts) -> ApiResult<Value> {
    let req = read_request(o)?;
    req.validate().await?;
    std::fs::create_dir_all(&o.jobs_dir)?;
    let id = new_job_id();
    let dir = o.jobs_dir.join(&id);
    std::fs::create_dir(&dir)?;
    std::fs::write(dir.join(REQUEST_FILE), serde_json::to_string_pretty(&req)?)?;
    let mut status = JobStatus {
        job_id: id.clone(),
        status: "queued".into(),
        created_utc: now_utc(),
        job_dir: absolute(&dir),
        ..Default::default()
    };
    write_status(&dir, &status)?;

    // Фоновый процесс: та же программа, команда `api run <id>`; текст расчёта — в output.txt.
    let exe = std::env::current_exe()?;
    let mut cmd = std::process::Command::new(exe);
    cmd.arg("api")
        .arg("--jobs-dir")
        .arg(&o.jobs_dir)
        .arg("run")
        .arg(&id)
        .stdin(Stdio::null())
        .stdout(std::fs::File::create(dir.join(OUTPUT_FILE))?)
        .stderr(std::fs::File::create(dir.join(WORKER_LOG))?);
    #[cfg(windows)]
    {
        use std::os::windows::process::CommandExt;
        // Без окна консоли и в своей группе процессов (Ctrl+C у вызывающего его не прервёт).
        const CREATE_NO_WINDOW: u32 = 0x0800_0000;
        const CREATE_NEW_PROCESS_GROUP: u32 = 0x0000_0200;
        cmd.creation_flags(CREATE_NO_WINDOW | CREATE_NEW_PROCESS_GROUP);
    }
    let mut child = match cmd.spawn() {
        Ok(c) => c,
        Err(e) => {
            status.status = "failed".into();
            status.error = Some(format!("cannot start the worker process: {e}"));
            write_status(&dir, &status)?;
            return Err(ApiError::internal(status.error.unwrap_or_default()));
        }
    };
    // Состояние могло уже обновиться фоновым процессом: pid дописываем, только пока оно queued.
    if let Ok(mut s) = read_status(&dir) {
        if s.status == "queued" {
            s.pid = Some(child.id());
            write_status(&dir, &s)?;
        }
    }
    if o.wait {
        child.wait()?;
        let s = read_status(&dir)?;
        let mut body = json!({ "job": s });
        if s.status == "done" {
            body["summary"] = summary_of(&s)?;
        }
        return Ok(body);
    }
    Ok(json!({ "job": read_status(&dir)? }))
}

/// Фоновый процесс задачи (`api run <id>`).
async fn worker(o: &Opts) -> ApiResult<()> {
    let dir = o.job_dir()?;
    let req: Request = serde_json::from_str(&std::fs::read_to_string(dir.join(REQUEST_FILE))?)?;
    let mut status = read_status(&dir)?;
    status.status = "running".into();
    status.started_utc = Some(now_utc());
    status.pid = Some(std::process::id());
    write_status(&dir, &status)?;

    // Паника не должна оставить задачу «running» навсегда.
    {
        let dir = dir.clone();
        let base = status.clone();
        let prev = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            let mut s = base.clone();
            s.status = "failed".into();
            s.finished_utc = Some(now_utc());
            s.error = Some(format!("internal error: {info}"));
            let _ = write_status(&dir, &s);
            prev(info);
        }));
    }

    let started = std::time::Instant::now();
    let finish = |status: &mut JobStatus| {
        status.finished_utc = Some(now_utc());
        status.elapsed_s = Some((started.elapsed().as_secs_f64() * 1000.0).round() / 1000.0);
    };
    let cli = match Cli::try_parse_from(req.to_args(&dir)) {
        Ok(c) => c,
        Err(e) => {
            status.status = "failed".into();
            status.error = Some(e.to_string());
            finish(&mut status);
            write_status(&dir, &status)?;
            return Err(ApiError::bad(e.to_string()));
        }
    };
    match execute(&cli).await {
        Ok((result, flight)) => {
            status.flight_log = Some(absolute(&flight));
            let traj = flight.with_file_name(
                flight.file_name().map(|n| n.to_string_lossy().replacen("flight-", "trajectory-", 1)).unwrap_or_default(),
            );
            status.trajectory_log = traj.exists().then(|| absolute(&traj));
            match result {
                Ok(outcome) => {
                    status.status = "done".into();
                    status.outcome = Some(outcome.to_string());
                }
                Err(e) => {
                    status.status = "failed".into();
                    status.error = Some(format!("{e:#}"));
                }
            }
        }
        Err(e) => {
            status.status = "failed".into();
            status.error = Some(format!("{e:#}"));
        }
    }
    finish(&mut status);
    write_status(&dir, &status)?;
    Ok(())
}

fn read_json(path: &str) -> ApiResult<Value> {
    let text = std::fs::read_to_string(path).map_err(|e| ApiError::internal(format!("cannot read {path}: {e}")))?;
    Ok(serde_json::from_str(&text)?)
}

/// Сводка из журнала полёта: разделы run, summary, launch_site, landing_site, files.
fn summary_of(s: &JobStatus) -> ApiResult<Value> {
    let path = s.flight_log.as_deref().ok_or_else(|| ApiError::not_ready("the job has no flight log yet"))?;
    let log = read_json(path)?;
    let mut out = Map::new();
    for k in ["run", "summary", "launch_site", "landing_site", "files"] {
        if let Some(v) = log.get(k) {
            out.insert(k.into(), v.clone());
        }
    }
    Ok(Value::Object(out))
}

fn result(o: &Opts) -> ApiResult<Value> {
    let dir = o.job_dir()?;
    let s = read_status(&dir)?;
    let part = o.part.as_str();
    // Текстовый вывод и полный журнал доступны и у неудавшейся задачи (для разбора ошибки).
    let diagnostic = matches!(part, "output" | "log");
    match s.status.as_str() {
        "done" => {}
        "failed" if diagnostic => {}
        "failed" => {
            return Err(ApiError {
                code: "job_failed",
                message: s.error.clone().unwrap_or_else(|| "the job failed".into()),
                exit: 1,
            });
        }
        other => return Err(ApiError::not_ready(format!("job is {other}; poll 'status' until it is done"))),
    }
    let data = match part {
        "summary" => summary_of(&s)?,
        "trajectory" => match &s.trajectory_log {
            Some(p) => read_json(p)?,
            None => Value::Null,
        },
        "log" => match &s.flight_log {
            Some(p) => read_json(p)?,
            None => return Err(ApiError::not_ready("the job has no flight log")),
        },
        "output" => Value::String(std::fs::read_to_string(dir.join(OUTPUT_FILE)).unwrap_or_default()),
        "all" => json!({
            "summary": summary_of(&s)?,
            "trajectory": match &s.trajectory_log { Some(p) => read_json(p)?, None => Value::Null },
            "log": match &s.flight_log { Some(p) => read_json(p)?, None => Value::Null },
        }),
        other => {
            return Err(ApiError::bad(format!("unknown part '{other}' (summary, trajectory, log, output, all)")));
        }
    };
    Ok(json!({ "job": s, "part": part, "data": data }))
}

fn list(o: &Opts) -> ApiResult<Value> {
    let mut jobs = Vec::new();
    if let Ok(rd) = std::fs::read_dir(&o.jobs_dir) {
        for e in rd.flatten() {
            if let Ok(s) = read_status(&e.path()) {
                jobs.push(s);
            }
        }
    }
    jobs.sort_by(|a, b| b.created_utc.cmp(&a.created_utc));
    Ok(json!({ "jobs": jobs }))
}

fn delete(o: &Opts) -> ApiResult<Value> {
    let dir = o.job_dir()?;
    let s = read_status(&dir)?;
    if matches!(s.status.as_str(), "queued" | "running") {
        return Err(ApiError::bad(format!("job is {}; wait until it finishes", s.status)));
    }
    std::fs::remove_dir_all(&dir)?;
    Ok(json!({ "deleted": s.job_id }))
}

// ---------------------------------------------------------------------------
// Справочники и описание
// ---------------------------------------------------------------------------

async fn catalog() -> ApiResult<Value> {
    let catalog = SpaceportCatalog::load(DEFAULT_SPACEPORTS).await.map_err(|e| ApiError::internal(e.to_string()))?;
    let (library, rockets) = tokio::task::spawn_blocking(|| {
        Ok::<_, launch_window_prediction::Error>((
            SpacecraftLibrary::load(DEFAULT_SPACECRAFT_DIR)?,
            LaunchVehicleLibrary::load(DEFAULT_ROCKET_DIR)?,
        ))
    })
    .await
    .map_err(|e| ApiError::internal(e.to_string()))?
    .map_err(|e| ApiError::internal(e.to_string()))?;
    let rockets: Vec<Value> = rockets
        .entries()
        .iter()
        .map(|e| match rockets.resolve(&e.id, &LaunchVehicleParams::default()) {
            Ok(r) => json!({ "id": e.id, "is_preset": e.is_preset, "based_on": e.based_on, "resolved": r }),
            Err(err) => json!({ "id": e.id, "is_preset": e.is_preset, "error": err.to_string() }),
        })
        .collect();
    let spacecraft: Vec<Value> = library
        .entries()
        .iter()
        .map(|e| match library.resolve(&e.id, &SpacecraftParams::default()) {
            Ok(r) => json!({ "id": e.id, "entry": e, "resolved": r }),
            Err(err) => json!({ "id": e.id, "entry": e, "error": err.to_string() }),
        })
        .collect();
    Ok(json!({
        "spaceports": catalog.iter().collect::<Vec<_>>(),
        "rockets": rockets,
        "spacecraft": spacecraft,
    }))
}

fn help() -> Value {
    json!({
        "usage": "launch-window-prediction api [--jobs-dir DIR] [--pretty] <command> [args]",
        "commands": {
            "help": "this description",
            "catalog": "spaceports, rockets and spacecraft with all parameters (ids for requests)",
            "submit [--input FILE] [--wait]": "start a calculation; the request JSON comes on standard input (or from FILE). Returns the job at once (status queued); with --wait, waits and returns the finished job and its summary",
            "status <job_id>": "job state: queued, running, done, failed",
            "result <job_id> [--part summary|trajectory|log|output|all]": "results of a finished job: summary (dates, delta-v, launch and landing sites), trajectory (3D trajectory, format lwp-trajectory/1), log (full flight log, format lwp-flight/2), output (text output), all",
            "list": "all jobs, newest first",
            "delete <job_id>": "remove a finished job and its files",
        },
        "exit_codes": { "0": "ok", "1": "internal error or the job failed", "2": "bad request", "3": "job not found", "4": "result not ready yet" },
        "response": "one JSON object: {\"ok\": true, \"api\": \"lwp-api/1\", ...} or {\"ok\": false, \"error\": {\"code\", \"message\"}}",
        "request_example": {
            "arrival_date": "2027-02-14",
            "spaceport": "baikonur",
            "rocket": "soyuz_21b_fregat",
            "landing_site": { "lat_deg": -89.0, "lon_deg": 45.0 },
        },
        "request_fields": {
            "launch_date | arrival_date": "UTC day YYYY-MM-DD (one of them); with landing_site, arrival_date is the landing day",
            "spaceport": "spaceport id (see catalog); omit to search all",
            "rocket, rocket_overrides": "rocket id and {field: value} overrides",
            "spacecraft, spacecraft_overrides": "spacecraft id and {field: value} overrides",
            "landing_site": "{lat_deg, lon_deg} on the Moon",
            "height_m": "surface height without terrain maps",
            "visibility": "{start, hours, step_min}: Sun/Earth visibility series from the landing site",
            "kernels_dir, download_kernels, terrain_dir": "data folders (defaults data/kernels, data/terrain; terrain_dir \"\" disables terrain)",
            "transfer_config, landing_config, visibility_config": "own parameter files",
            "no_cache": "do not use the ascent cache",
        },
    })
}
