use std::path::PathBuf;
use std::str::FromStr;
use std::sync::Arc;

use anise::prelude::Almanac;
use anise::time::{Duration, Epoch, Unit};
use clap::Parser;
use launch_window_prediction::ephemeris::{AlmanacSource, load_almanac};
use launch_window_prediction::landing::{self, LandingConfig, Verdict};
use launch_window_prediction::lunar_landing::{DescentModel, DescentVehicle, LandingTarget};
use launch_window_prediction::refine::refine;
use launch_window_prediction::sky::MOON_RADIUS_KM;
use launch_window_prediction::launch_vehicle::{
    LaunchCapability, LaunchVehicleLibrary, LaunchVehicleParams, ResolvedLaunchVehicle,
};
use launch_window_prediction::flightlog;
use launch_window_prediction::runlog::{RunLog, iso_short};
use launch_window_prediction::sky::sky_series;
use launch_window_prediction::spacecraft::{ResolvedSpacecraft, Severity, SpacecraftLibrary, SpacecraftParams};
use launch_window_prediction::spaceport::{Spaceport, SpaceportCatalog};
use launch_window_prediction::terrain::TerrainSet;
use launch_window_prediction::trajectory::{
    DateConstraint, DeltaVBudget, Maneuver, RefineAttempt, ScanPoint, TrajectoryKind, TransferConfig, TransferReport,
    TransferSolution,
};
use launch_window_prediction::visibility::{
    SiteVisibility, VisibilityConfig, VisibilitySample, VisibilitySummary,
};
use serde_json::json;

mod api;

/// Файлы параметров по умолчанию; если файла нет, берутся встроенные значения.
const DEFAULT_VISIBILITY_CONFIG: &str = "data/visibility.toml";
const DEFAULT_LANDING_CONFIG: &str = "data/landing.toml";
const DEFAULT_TRANSFER_CONFIG: &str = "data/transfer.toml";
/// Сколько вариантов прилёта проверить численным уточнением, если лучший не подошёл.
const MAX_REFINE_ATTEMPTS: usize = 8;

/// Lunch Window Prediction: даты старта и прилёта Земля -> Луна, проверка площадки,
/// видимость Солнца и Земли.
///
/// Порядок расчёта: площадка (если задана точка и карты) -> перелёт (если задана дата)
/// -> видимость с момента прилёта (или с --start).
#[derive(Debug, Parser)]
struct Cli {
    /// Файл с космодромами.
    #[arg(long, default_value = "data/spaceports.toml")]
    spaceports: PathBuf,

    /// Космодром (id из файла космодромов). Без него перелёт считается со всех.
    #[arg(long)]
    spaceport: Option<String>,

    /// Дата старта (UTC), ГГГГ-ММ-ДД: найти лучший старт в эти сутки и дату прилёта.
    #[arg(long, conflicts_with = "arrival_date")]
    launch_date: Option<String>,

    /// Искать старт только после этого времени суток (UTC), ЧЧ:ММ. Только с --launch-date.
    #[arg(long, requires = "launch_date", value_name = "HH:MM")]
    launch_from_utc: Option<String>,

    /// Искать старт только до этого времени суток (UTC), ЧЧ:ММ (24:00 — конец суток). Только с --launch-date.
    #[arg(long, requires = "launch_date", value_name = "HH:MM")]
    launch_to_utc: Option<String>,

    /// Дата прилёта (UTC), ГГГГ-ММ-ДД. С точкой посадки (--lat/--lon) — сутки посадки в эту точку,
    /// без неё — сутки прилёта к Луне. Находит лучший вариант и дату старта.
    #[arg(long)]
    arrival_date: Option<String>,

    /// Папка с локальными ядрами (.bsp/.bpc/.pca/.epa).
    #[arg(long, conflicts_with = "download_kernels")]
    kernels_dir: Option<PathBuf>,

    /// Скачать актуальные ядра через MetaAlmanac (нужен интернет).
    #[arg(long)]
    download_kernels: bool,

    /// Широта точки посадки на Луне, градусы (+ север).
    #[arg(long, requires = "lon", allow_negative_numbers = true)]
    lat: Option<f64>,

    /// Долгота точки посадки на Луне, градусы (+ восток).
    #[arg(long, requires = "lat", allow_negative_numbers = true)]
    lon: Option<f64>,

    /// Папка с картами высот LOLA (.LBL + .IMG). Нужна для проверки площадки;
    /// без неё горизонт считается ровным.
    #[arg(long)]
    terrain_dir: Option<PathBuf>,

    /// Высота поверхности над опорной сферой, м (только без --terrain-dir).
    #[arg(long, default_value_t = 0.0, allow_negative_numbers = true)]
    height_m: f64,

    /// Файл с требованиями к площадке.
    #[arg(long)]
    landing_config: Option<PathBuf>,

    /// Файл с параметрами видимости.
    #[arg(long)]
    visibility_config: Option<PathBuf>,

    /// Файл с параметрами перелёта.
    #[arg(long)]
    transfer_config: Option<PathBuf>,

    /// Тип траектории перелёта: direct (прямой перелёт, задача Ламберта) или hohmann
    /// (перелёт Гомана: разгон по касательной, апогей у Луны, 180°).
    #[arg(long, default_value = "direct", value_name = "direct|hohmann")]
    trajectory: TrajectoryKind,

    /// Начало расчёта видимости, например "2026-11-14T12:00:00 UTC".
    /// По умолчанию — момент прилёта найденного перелёта.
    #[arg(long, requires = "lat")]
    start: Option<String>,

    /// Длительность расчёта видимости, часы (0 = один момент времени).
    #[arg(long, default_value_t = 0.0)]
    hours: f64,

    /// Шаг по времени для видимости, минуты.
    #[arg(long, default_value_t = 60.0)]
    step_min: f64,

    /// Папка для JSON-журналов расчёта (файлы flight-<время>.json и trajectory-<время>.json на каждый запуск).
    #[arg(long, default_value = "logs")]
    log_dir: PathBuf,

    /// Папка кэша моделирования подъёма ракеты (ускоряет повторные запуски с той же ракетой).
    #[arg(long, default_value = "cache")]
    cache_dir: PathBuf,

    /// Не использовать кэш подъёма (всё считать заново).
    #[arg(long)]
    no_cache: bool,

    /// Лунный аппарат: id шаблона или сохранённого профиля (см. --list-spacecraft).
    /// С ним перелёт учитывает Δv аппарата.
    #[arg(long)]
    spacecraft: Option<String>,

    /// Переопределить параметр аппарата: поле=значение (можно несколько раз),
    /// например --sc payload_mass_kg=500 --sc isp_s=320.
    #[arg(long = "sc", value_name = "FIELD=VALUE", requires = "spacecraft")]
    spacecraft_overrides: Vec<String>,

    /// Сохранить аппарат (основа + переопределения) как пользовательский профиль с этим id.
    #[arg(long, value_name = "ID", requires = "spacecraft")]
    save_spacecraft: Option<String>,

    /// Показать каталог лунных аппаратов и ракет-носителей (шаблоны и свои профили) и выйти.
    #[arg(long, visible_alias = "list-vehicles")]
    list_spacecraft: bool,

    /// Ракета-носитель: id из каталога или своего профиля. С ней перелёт учитывает её
    /// грузоподъёмность к Луне, разгонную ступень и космодромы.
    #[arg(long)]
    rocket: Option<String>,

    /// Переопределить параметр ракеты: поле=значение (можно несколько раз),
    /// например --rk payload_tli_kg=9000. Ракета с переопределениями (как и свой профиль)
    /// стартует с любого космодрома.
    #[arg(long = "rk", value_name = "FIELD=VALUE", requires = "rocket")]
    rocket_overrides: Vec<String>,

    /// Сохранить ракету (основа + переопределения) как свой профиль с этим id.
    #[arg(long, value_name = "ID", requires = "rocket")]
    save_rocket: Option<String>,

    /// Папка своих профилей ракет.
    #[arg(long, default_value = "data/launch_vehicles/custom")]
    rocket_dir: PathBuf,

    /// Папка пользовательских профилей аппаратов.
    #[arg(long, default_value = "data/spacecraft/custom")]
    spacecraft_dir: PathBuf,
}

/// Время суток ЧЧ:ММ -> часы (24:00 допустимо).
fn parse_hhmm(s: &str) -> anyhow::Result<f64> {
    let (h, m) = s.split_once(':').ok_or_else(|| anyhow::anyhow!("time of day must be HH:MM, got {s:?}"))?;
    let (h, m): (u32, u32) = (h.trim().parse()?, m.trim().parse()?);
    anyhow::ensure!(m < 60 && (h < 24 || (h == 24 && m == 0)), "time of day out of range: {s:?}");
    Ok(f64::from(h) + f64::from(m) / 60.0)
}

/// Явно указанный файл обязателен; файл по умолчанию — только если существует.
async fn config_path(explicit: Option<PathBuf>, default: &str) -> Option<PathBuf> {
    match explicit {
        Some(path) => Some(path),
        None => {
            let default = PathBuf::from(default);
            if tokio::fs::try_exists(&default).await.unwrap_or(false) {
                Some(default)
            } else {
                tracing::info!("{} not found, using built-in defaults", default.display());
                None
            }
        }
    }
}

fn parse_date(text: &str) -> anyhow::Result<Epoch> {
    let parts: Vec<&str> = text.trim().split('-').collect();
    let bad = || anyhow::anyhow!("date '{text}' must be YYYY-MM-DD");
    anyhow::ensure!(parts.len() == 3, bad());
    let y: i32 = parts[0].parse().map_err(|_| bad())?;
    let m: u8 = parts[1].parse().map_err(|_| bad())?;
    let d: u8 = parts[2].parse().map_err(|_| bad())?;
    anyhow::ensure!((1..=12).contains(&m) && (1..=31).contains(&d), bad());
    Ok(Epoch::from_gregorian_utc_at_midnight(y, m, d))
}

fn hours(d: Duration) -> f64 {
    d.to_seconds() / 3600.0
}

/// Имя файла без папки.
fn file_name(path: &std::path::Path) -> String {
    path.file_name().map_or_else(|| path.display().to_string(), |n| n.to_string_lossy().into_owned())
}

/// Момент, округлённый до секунды, для вывода.
fn ts(e: Epoch) -> String {
    format!("{}", e.round(1 * Unit::Second))
}

fn almanac_source(cli: &Cli) -> Option<AlmanacSource> {
    match (&cli.kernels_dir, cli.download_kernels) {
        (Some(dir), _) => Some(AlmanacSource::LocalDir(dir.clone())),
        (None, true) => Some(AlmanacSource::Latest),
        (None, false) => None,
    }
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    // Машинный интерфейс: `launch-window-prediction api <команда>` — JSON на входе и выходе.
    if args.get(1).map(String::as_str) == Some("api") {
        init_tracing(true);
        let code = api::main(&args[2..]).await;
        std::process::exit(code);
    }
    init_tracing(false);
    let cli = Cli::parse();
    let (result, _) = execute(&cli).await?;
    result.map(|_| ())
}

/// Журнал работы программы (tracing). В режиме API — в поток ошибок, чтобы в выводе был только JSON.
fn init_tracing(to_stderr: bool) {
    let filter = tracing_subscriber::EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info"));
    if to_stderr {
        tracing_subscriber::fmt().with_env_filter(filter).with_writer(std::io::stderr).with_ansi(false).init();
    } else {
        tracing_subscriber::fmt().with_env_filter(filter).init();
    }
}

/// Один расчёт целиком: журнал, расчёт, итог в журнале. Возвращает итог расчёта и путь
/// журнала информации о полёте (журнал траектории лежит рядом, см. [`RunLog::trajectory_path`]).
async fn execute(cli: &Cli) -> anyhow::Result<(anyhow::Result<&'static str>, PathBuf)> {
    let log = Arc::new(RunLog::create(&cli.log_dir)?);
    let started = std::time::Instant::now();
    let run_info = json!({
        "program": env!("CARGO_PKG_NAME"),
        "version": env!("CARGO_PKG_VERSION"),
        "args": std::env::args().collect::<Vec<_>>(),
        "cwd": std::env::current_dir().ok(),
    });
    log.event("run_start", &run_info);
    log.section("run", &run_info);
    log.merge("run", json!({ "started_utc": iso_short(&Epoch::now().unwrap_or_default()) }));
    log.merge("files", json!({ "flight_log": file_name(log.path()) }));
    let result = run(cli, log.clone()).await;
    let elapsed = started.elapsed().as_secs_f64();
    let end = match &result {
        Ok(outcome) => json!({ "ok": true, "outcome": outcome, "elapsed_s": elapsed }),
        Err(e) => json!({ "ok": false, "error": format!("{e:#}"), "elapsed_s": elapsed }),
    };
    log.event("run_end", &end);
    log.merge("run", json!({ "finished_utc": iso_short(&Epoch::now().unwrap_or_default()) }));
    log.merge("run", end);
    log.flush();
    println!("Flight log: {}", log.path().display());
    Ok((result, log.path().to_path_buf()))
}

async fn run(cli: &Cli, log: Arc<RunLog>) -> anyhow::Result<&'static str> {
    let catalog = SpaceportCatalog::load(&cli.spaceports).await?;
    tracing::info!(count = catalog.len(), "spaceports loaded");
    log.event(
        "spaceports",
        json!({ "path": cli.spaceports, "spaceports": catalog.iter().collect::<Vec<_>>() }),
    );
    if let Some(id) = &cli.spaceport {
        anyhow::ensure!(catalog.get(id).is_some(), "unknown spaceport '{id}'");
    }

    // 0. Лунный аппарат: шаблон/профиль + переопределения, проверки, сохранение.
    let (library, rockets) = {
        let (sd, rd) = (cli.spacecraft_dir.clone(), cli.rocket_dir.clone());
        tokio::task::spawn_blocking(move || {
            Ok::<_, launch_window_prediction::Error>((SpacecraftLibrary::load(sd)?, LaunchVehicleLibrary::load(rd)?))
        })
        .await??
    };
    if cli.list_spacecraft {
        print_spacecraft_list(&library, &log);
        print_rocket_list(&rockets, &log);
        return Ok("done: vehicle catalog");
    }
    let spacecraft = match &cli.spacecraft {
        Some(id) => {
            let mut overrides = SpacecraftParams::default();
            for a in &cli.spacecraft_overrides {
                overrides.set_from_str(a)?;
            }
            let resolved = library.resolve(id, &overrides)?;
            print_spacecraft(&resolved);
            log.event(
                "spacecraft",
                json!({ "requested_id": id, "overrides": overrides, "resolved": resolved }),
            );
            if let Some(new_id) = &cli.save_spacecraft {
                let path = library.save_custom(new_id, id, &overrides)?;
                println!("Spacecraft profile saved: {} (id '{new_id}')", path.display());
                log.event("spacecraft_saved", json!({ "id": new_id, "based_on": id, "overrides": overrides, "path": path }));
            }
            if resolved.has_errors() {
                println!("Spacecraft has critical errors; further calculations are skipped.");
                return Ok("stopped: spacecraft configuration has critical errors");
            }
            Some(resolved)
        }
        None => None,
    };

    // 0b. Ракета-носитель: каталог + переопределения, грузоподъёмность для этого аппарата.
    let rocket = match &cli.rocket {
        Some(id) => {
            let mut overrides = LaunchVehicleParams::default();
            for a in &cli.rocket_overrides {
                overrides.set_from_str(a)?;
            }
            let resolved = rockets.resolve(id, &overrides)?;
            if !resolved.config.any_spaceport {
                for p in &resolved.config.spaceports {
                    anyhow::ensure!(catalog.get(p).is_some(), "rocket '{id}': unknown spaceport '{p}'");
                }
            }
            let payload = spacecraft.as_ref().map_or(resolved.config.payload_tli_kg, |s| s.derived.initial_mass_kg);
            let cap = resolved.config.capability(payload);
            print_rocket(&resolved, &cap, spacecraft.is_some());
            log.event(
                "launch_vehicle",
                json!({ "requested_id": id, "overrides": overrides, "resolved": resolved, "capability": cap }),
            );
            if let Some(new_id) = &cli.save_rocket {
                let path = rockets.save_custom(new_id, id, &overrides)?;
                println!("Rocket profile saved: {} (id '{new_id}')", path.display());
                log.event("launch_vehicle_saved", json!({ "id": new_id, "based_on": id, "overrides": overrides, "path": path }));
            }
            if !resolved.errors.is_empty() {
                println!("Rocket has critical errors; further calculations are skipped.");
                return Ok("stopped: rocket configuration has critical errors");
            }
            if !cap.within_rated_payload {
                println!(
                    "The spacecraft ({:.0} kg) is heavier than the rocket can send to the Moon ({:.0} kg); \
                     further calculations are skipped.",
                    cap.payload_kg, cap.payload_tli_kg
                );
                return Ok("stopped: spacecraft heavier than rocket TLI payload");
            }
            Some((resolved, cap))
        }
        None => None,
    };
    if let (Some((r, _)), Some(id)) = (&rocket, &cli.spaceport) {
        anyhow::ensure!(
            r.config.any_spaceport || r.config.spaceports.contains(id),
            "rocket '{}' does not launch from '{id}' (its spaceports: {})",
            r.config.name,
            r.config.spaceports.join(", ")
        );
    }

    let date = match (&cli.launch_date, &cli.arrival_date) {
        (Some(d), _) => Some(DateConstraint::LaunchDay(parse_date(d)?)),
        (None, Some(d)) => Some(DateConstraint::ArrivalDay(parse_date(d)?)),
        (None, None) => None,
    };
    let site = cli.lat.zip(cli.lon);
    log.event(
        "request",
        json!({
            "spaceport": cli.spaceport,
            "date_constraint": date,
            "landing_site": site.map(|(lat, lon)| json!({ "lat_deg": lat, "lon_deg": lon })),
            "terrain_dir": cli.terrain_dir,
            "kernels_dir": cli.kernels_dir,
            "download_kernels": cli.download_kernels,
            "visibility_start": cli.start,
            "visibility_hours": cli.hours,
            "visibility_step_min": cli.step_min,
            "flat_height_m": cli.height_m,
        }),
    );
    if let Some((lat, lon)) = site {
        log.merge("landing_site", flightlog::landing_request(lat, lon));
    }
    if let Some(p) = cli.spaceport.as_deref().and_then(|id| catalog.get(id)) {
        log.section("launch_site", flightlog::spaceport_info(p));
    }

    // Карты рельефа: чтение — блокирующий I/O, уводим из async-рантайма.
    let terrain = match cli.terrain_dir.clone() {
        Some(dir) => {
            let t = Arc::new(tokio::task::spawn_blocking(move || TerrainSet::open_dir(dir)).await??);
            log.event(
                "terrain_maps",
                json!({
                    "dir": cli.terrain_dir,
                    "maps": t.maps().iter().map(|m| json!({ "name": m.name(), "meta": m.meta() })).collect::<Vec<_>>(),
                }),
            );
            Some(t)
        }
        None => None,
    };

    // 1. Пригодность площадки. Если сесть нельзя — дальше не считаем.
    if let Some((lat, lon)) = site {
        println!("== Landing site: lat {lat:.4} deg, lon {lon:.4} deg");
        if check_landing(cli, &log, terrain.clone(), lat, lon).await? == Some(Verdict::Unsuitable) {
            println!("Landing here is not possible; further calculations are skipped.");
            log.section(
                "summary",
                json!({ "status": "landing site unsuitable", "details": "see landing_site.terrain_check" }),
            );
            return Ok("stopped: landing site unsuitable");
        }
    }

    let needs_kernels = date.is_some() || cli.start.is_some();
    if !needs_kernels {
        if let Some(source) = almanac_source(cli) {
            load_almanac(source).await?;
            tracing::info!("almanac loaded");
            log.event("almanac_loaded", json!({ "source": format!("{:?}", almanac_source(cli)) }));
        }
        return Ok("done");
    }
    let Some(source) = almanac_source(cli) else {
        anyhow::bail!("trajectory and visibility need kernels: add --download-kernels or --kernels-dir");
    };
    let almanac = Arc::new(load_almanac(source.clone()).await?);
    tracing::info!("almanac loaded");
    log.event("almanac_loaded", json!({ "source": format!("{source:?}") }));

    // Точка посадки с рельефом и правилами видимости: нужна и перелёту (условия посадки),
    // и расчёту видимости.
    let site_vis = match site {
        Some((lat, lon)) => {
            let vis_path = config_path(cli.visibility_config.clone(), DEFAULT_VISIBILITY_CONFIG).await;
            let vis_cfg = match &vis_path {
                Some(path) => VisibilityConfig::load(path).await?,
                None => VisibilityConfig::default(),
            };
            log.event(
                "config_visibility",
                json!({ "path": vis_path, "built_in_defaults": vis_path.is_none(), "config": vis_cfg }),
            );
            let flat_height_m = cli.height_m;
            let t = terrain.clone();
            let sv = tokio::task::spawn_blocking(move || {
                SiteVisibility::new(t.as_deref(), lat, lon, flat_height_m, vis_cfg)
            })
            .await??;
            log.merge("landing_site", flightlog::landing_surface(&sv));
            Some(sv)
        }
        None => None,
    };

    // 2. Перелёт (с точкой посадки — до посадки в неё).
    let mut arrival = None;
    if let Some(constraint) = date {
        let ports: Vec<Spaceport> = match (&cli.spaceport, &rocket) {
            (Some(id), _) => vec![catalog.get(id).cloned().expect("checked above")],
            (None, Some((r, _))) if r.config.any_spaceport => catalog.iter().cloned().collect(),
            (None, Some((r, _))) => r.config.spaceports.iter().filter_map(|p| catalog.get(p).cloned()).collect(),
            (None, None) => catalog.iter().cloned().collect(),
        };
        let reports = solve_transfers(
            cli,
            &log,
            almanac.clone(),
            ports,
            constraint,
            spacecraft.as_ref(),
            rocket.as_ref(),
            site_vis.as_ref(),
        )
        .await?;
        print_transfers(&reports, constraint);
        if let Some(best) = reports.iter().find(|r| r.solution.is_some()).or(reports.first()) {
            log.section("launch_site", flightlog::launch_site(best, Some(&almanac)));
            if site.is_some() {
                log.merge("landing_site", flightlog::landing_result(best));
            }
        }
        log.section(
            "summary",
            flightlog::summary(
                &reports,
                spacecraft.as_ref().map(|sc| sc.config.name.as_str()),
                rocket.as_ref().map(|(r, _)| r.config.name.as_str()),
            ),
        );
        if let (Some((r, cap)), Some(s)) = (&rocket, reports.first().and_then(|r| r.solution.as_ref())) {
            // Мгновенный эквивалент разгона: уточнённый, если есть, иначе по формулам.
            let dv = s.refined.as_ref().and_then(|n| n.tli_burn.as_ref()).map_or(s.departure.delta_v_magnitude_km_s, |b| b.impulsive_delta_v_km_s) * 1000.0;
            if let Some(b) = s.refined.as_ref().and_then(|n| n.tli_burn.as_ref()) {
                println!(
                    "  Rocket {}: TLI {:.0} m/s (impulsive equivalent) of {:.0} m/s rated for {:.0} kg; simulated burn {:.1} min, {:.0} m/s from the engine with gravity loss, {:.0} of {:.0} kg propellant left in the upper stage.",
                    r.config.name,
                    dv,
                    cap.max_tli_delta_v_m_s,
                    cap.payload_kg,
                    b.duration_s / 60.0,
                    b.delta_v_km_s * 1000.0,
                    b.propellant_used_kg,
                    b.engine.propellant_kg,
                );
            }
            let burn = r.config.tli_burn(cap.payload_kg, dv);
            let simulated = s.refined.as_ref().is_some_and(|n| n.tli_burn.is_some());
            if !simulated { println!(
                "  Rocket {}: TLI {:.0} m/s of {:.0} m/s possible for {:.0} kg, burn {:.1} min, {:.0} of {:.0} kg upper-stage propellant.",
                r.config.name,
                dv,
                cap.max_tli_delta_v_m_s,
                cap.payload_kg,
                burn.duration_s / 60.0,
                burn.propellant_used_kg,
                r.config.upper_stage_propellant_kg,
            ); }
            if !simulated && burn.exceeds_stage_propellant {
                println!(
                    "  NOTE: from a circular parking orbit the upper stage would need more propellant than it holds; \
                     the rated TLI payload implies the burn starts from a higher (elliptical) orbit reached by the \
                     lower stages. Burn time and propellant are estimates."
                );
            }
            log.event(
                "launch_vehicle_tli",
                json!({
                    "rocket": r.config.name,
                    "payload_kg": cap.payload_kg,
                    "tli_delta_v_m_s": dv,
                    "max_tli_delta_v_m_s": cap.max_tli_delta_v_m_s,
                    "burn_estimate": burn,
                    "burn_simulated": s.refined.as_ref().and_then(|n| n.tli_burn.as_ref()),
                    "upper_stage_propellant_kg": r.config.upper_stage_propellant_kg,
                }),
            );
        }
        // Видимость считается с момента посадки (или прилёта, если точки нет).
        arrival = reports
            .iter()
            .find_map(|r| r.solution.as_ref())
            .map(|s| {
                // Уточнённое касание, иначе касание по формулам, иначе прилёт.
                let refined = s.refined.as_ref();
                refined
                    .and_then(|n| n.landing.as_ref().map(|l| l.touchdown_epoch))
                    .or(s.landing.map(|l| l.epoch))
                    .or(refined.map(|n| n.arrival_epoch))
                    .unwrap_or(s.arrival.epoch)
            });
    }

    // 3. Видимость Солнца и Земли.
    let (Some((lat, lon)), Some(site_vis)) = (site, site_vis) else {
        return Ok("done");
    };
    let start = match cli.start.as_deref() {
        Some(text) => Epoch::from_str(text)
            .map_err(|e| anyhow::anyhow!("cannot parse --start '{text}': {e}"))?,
        None => match arrival {
            Some(t) => t,
            None => return Ok("done: no transfer found, visibility skipped"),
        },
    };
    anyhow::ensure!(cli.hours >= 0.0, "--hours must be >= 0");
    anyhow::ensure!(cli.step_min > 0.0, "--step-min must be > 0");
    let end = start + cli.hours * Unit::Hour;
    let step = cli.step_min * Unit::Minute;

    let flat_height_m = cli.height_m;
    let observer = site_vis.site;
    let a = almanac.clone();
    let sky = tokio::task::spawn_blocking(move || sky_series(&a, &observer, start, end, step)).await??;
    let vis = site_vis.evaluate_all(&sky);
    let summary = VisibilitySummary::from_samples(&vis, step);
    let mask = site_vis.horizon_mask_deg();
    log.event(
        "visibility_site",
        json!({
            "lat_deg": lat,
            "lon_deg": lon,
            "observer": site_vis.site,
            "terrain": site_vis.terrain,
            "config": site_vis.config,
            "start": launch_window_prediction::runlog::iso(&start),
            "end": launch_window_prediction::runlog::iso(&end),
            "step_s": step.to_seconds(),
            "horizon_mask_azimuth_step_deg": mask.map(|m| 360.0 / m.len() as f64),
            "horizon_mask_deg": mask,
        }),
    );
    for v in &vis {
        log.event("visibility_sample", v);
    }
    log.event("visibility_summary", &summary);
    log.merge("landing_site", flightlog::landing_visibility(start, end, step.to_seconds(), &summary));
    println!("== Visibility from {}", ts(start));
    print_visibility(&site_vis, flat_height_m, &vis, &summary);
    if let Some(sc) = &spacecraft {
        let gap_h = hours(summary.longest_sun_gap);
        let verdict = match sc.derived.eclipse_autonomy_h {
            Some(a) if a >= gap_h => {
                println!("Spacecraft eclipse autonomy {a:.1} h covers the longest sunless gap {gap_h:.1} h.");
                "covers"
            }
            Some(a) => {
                println!("WARNING: spacecraft eclipse autonomy {a:.1} h is shorter than the longest sunless gap {gap_h:.1} h.");
                "insufficient"
            }
            None => {
                println!("Spacecraft eclipse autonomy unknown (set battery_capacity_wh and idle_power_consumption_w).");
                "unknown"
            }
        };
        let check = json!({
            "eclipse_autonomy_h": sc.derived.eclipse_autonomy_h,
            "longest_sun_gap_h": gap_h,
            "verdict": verdict,
        });
        log.event("spacecraft_eclipse_check", &check);
        log.merge("landing_site", json!({ "spacecraft_eclipse_check": check }));
    }
    Ok("done")
}

async fn check_landing(
    cli: &Cli,
    log: &RunLog,
    terrain: Option<Arc<TerrainSet>>,
    lat: f64,
    lon: f64,
) -> anyhow::Result<Option<Verdict>> {
    let Some(terrain) = terrain else {
        println!("Landing: not checked (no --terrain-dir)");
        log.event("landing", json!({ "lat_deg": lat, "lon_deg": lon, "checked": false, "reason": "no --terrain-dir" }));
        log.merge("landing_site", json!({ "terrain_check": { "checked": false, "reason": "no --terrain-dir" } }));
        return Ok(None);
    };
    let path = config_path(cli.landing_config.clone(), DEFAULT_LANDING_CONFIG).await;
    let cfg = match &path {
        Some(path) => LandingConfig::load(path).await?,
        None => LandingConfig::default(),
    };
    log.event("config_landing", json!({ "path": path, "built_in_defaults": path.is_none(), "config": cfg }));
    let c = cfg.clone();
    let a = tokio::task::spawn_blocking(move || landing::assess(&terrain, lat, lon, &c)).await??;
    log.event(
        "landing",
        json!({ "lat_deg": lat, "lon_deg": lon, "checked": true, "explanation": a.explain(&cfg), "assessment": a }),
    );
    log.merge("landing_site", flightlog::landing_terrain(&a, &cfg));
    println!("Landing: {}", a.explain(&cfg));
    if !matches!(a.verdict, Verdict::NotChecked { .. }) {
        println!(
            "  area radius {:.0} m, {} points, map {} ({:.0} m/pix), slope baseline {:.0} m",
            cfg.check_radius_m,
            a.points,
            a.map_name.as_deref().unwrap_or("?"),
            a.resolution_m,
            a.baseline_m,
        );
        println!(
            "  slope max {:.1} deg, mean {:.1} deg, relief {:.1} m{}",
            a.max_slope_deg,
            a.mean_slope_deg,
            a.relief_m,
            a.steepest_at
                .map(|(la, lo)| format!(", steepest at ({la:.5}, {lo:.5})"))
                .unwrap_or_default(),
        );
    }
    Ok(Some(a.verdict))
}

/// Считает перелёт для каждого космодрома; результат отсортирован по суммарному Δv.
async fn solve_transfers(
    cli: &Cli,
    log: &RunLog,
    almanac: Arc<Almanac>,
    ports: Vec<Spaceport>,
    constraint: DateConstraint,
    spacecraft: Option<&ResolvedSpacecraft>,
    rocket: Option<&(ResolvedLaunchVehicle, LaunchCapability)>,
    landing_site: Option<&SiteVisibility>,
) -> anyhow::Result<Vec<TransferReport>> {
    let path = config_path(cli.transfer_config.clone(), DEFAULT_TRANSFER_CONFIG).await;
    let cfg = match &path {
        Some(path) => TransferConfig::load(path).await?,
        None => TransferConfig::default(),
    };
    log.event("config_transfer", json!({ "path": path, "built_in_defaults": path.is_none(), "config": cfg }));
    let mut cfg = cfg;
    if let Some(t) = &cli.launch_from_utc {
        cfg.transfer.launch_from_utc_h = parse_hhmm(t)?;
    }
    if let Some(t) = &cli.launch_to_utc {
        cfg.transfer.launch_to_utc_h = parse_hhmm(t)?;
    }
    cfg.transfer.validate()?;
    if let Some((r, cap)) = rocket {
        // Параметры ракеты заменяют допущения из transfer.toml.
        let dt = &mut cfg.transfer;
        dt.max_tli_delta_v_km_s = cap.max_tli_delta_v_m_s / 1000.0;
        anyhow::ensure!(
            dt.max_tli_delta_v_km_s >= 3.0,
            "rocket {} can give at most {:.0} m/s for TLI with {:.0} kg: not enough to reach the Moon",
            r.config.name,
            cap.max_tli_delta_v_m_s,
            cap.payload_kg
        );
        if let Some(h) = r.config.parking_altitude_km {
            dt.parking_altitude_km = h;
        }
        if let Some(a) = r.config.ascent_duration_s {
            dt.coast_min_min = dt.coast_min_min.max(a / 60.0);
            if let Some(c) = r.config.max_coast_s {
                dt.coast_max_min = dt.coast_max_min.min((a + c) / 60.0);
            }
        }
        dt.validate()?;
        log.event("config_transfer_effective", json!({ "rocket": r.config.name, "config": cfg }));
    }
    let mut budget = spacecraft
        .map(|sc| DeltaVBudget::new(&sc.config.name, sc.derived.delta_v_m_s, &cfg.mission))
        .transpose()?;
    let descent = DescentModel {
        config: cfg.descent.clone(),
        vehicle: spacecraft.map(|sc| DescentVehicle {
            initial_mass_kg: sc.derived.initial_mass_kg,
            max_thrust_n: sc.config.max_thrust_n,
            isp_s: sc.config.isp_s,
        }),
    };
    // Δv спуска для бюджета: моделирование торможения по параметрам аппарата
    // (если в transfer.toml не задано число). До спуска аппарат уже потратил топливо на
    // торможение у Луны (типично ~850 м/с), если этот манёвр его.
    if let Some(b) = budget.as_mut() {
        if cfg.mission.descent_delta_v_m_s.is_none() {
            let r_llo = MOON_RADIUS_KM + cfg.transfer.lunar_orbit_altitude_km;
            let site_radius = landing_site.map_or(MOON_RADIUS_KM, |s| s.site.position_km().norm());
            let prior = if b.maneuvers.contains(&Maneuver::Loi) { 850.0 } else { 0.0 };
            match descent.plan(4902.8, r_llo, prior, site_radius)? {
                Ok(plan) => {
                    b.descent_delta_v_m_s = plan.total_delta_v_m_s;
                    log.event("descent_plan_for_budget", json!({ "prior_delta_v_m_s": prior, "plan": plan }));
                }
                Err(reason) => anyhow::bail!("spacecraft {}: {reason}", b.spacecraft_name),
            }
        }
        log.event("delta_v_budget", &*b);
    }
    if landing_site.is_some() {
        log.event("descent_model", &descent);
    }
    let target = landing_site.map(|s| LandingTarget { site: s.clone(), window: cfg.landing_window.clone(), descent });
    // Ракета для моделирования подъёма по космодромам: заданная (--rocket) или первая ракета
    // каталога с данными ступеней, которая летает с этого космодрома.
    let mut ascent = std::collections::HashMap::new();
    let presets = launch_window_prediction::launch_vehicle::LaunchVehicleLibrary::presets_only()?;
    for port in &ports {
        let chosen = match rocket {
            Some((rk, _)) => Some((rk.id.clone(), rk.config.clone())),
            None => presets.entries().iter().find_map(|e| {
                let rk = presets.resolve(&e.id, &Default::default()).ok()?;
                (rk.config.ascent.is_some() && rk.config.spaceports.iter().any(|p| p == &port.id))
                    .then(|| (rk.id.clone(), rk.config))
            }),
        };
        if let Some((id, cfg_rk)) = chosen {
            if let Some(vehicle) = cfg_rk.ascent.clone() {
                // Груз: лунный аппарат, а без него — паспортная масса к Луне.
                let payload_kg = spacecraft.map_or(cfg_rk.payload_tli_kg, |sc| sc.derived.initial_mass_kg);
                ascent.insert(
                    port.id.clone(),
                    launch_window_prediction::trajectory::AscentSetup {
                        rocket_id: id,
                        rocket_name: cfg_rk.name.clone(),
                        vehicle,
                        upper: cfg_rk.upper_stage(),
                        payload_kg,
                    },
                );
            }
        }
    }
    log.event("ascent_setup", &ascent);
    log.event(
        "trajectory_type",
        json!({ "trajectory": cli.trajectory, "name": cli.trajectory.report_name(), "model": cli.trajectory.description() }),
    );
    let trajectory = launch_window_prediction::trajectory::build(
        cli.trajectory,
        &cfg,
        ascent.clone(),
        (!cli.no_cache).then(|| cli.cache_dir.clone()),
    )?;
    let numerical = cfg.numerical.clone();
    // Двигатель разгона к Луне: аппарат, если разгон — его манёвр, иначе разгонная ступень
    // ракеты космодрома (с остатком топлива после подъёма). Без данных — мгновенный импульс.
    let sc_engine = spacecraft.map(|sc| {
        (sc.config.name.clone(), sc.config.max_thrust_n, sc.config.isp_s, sc.derived.initial_mass_kg, sc.config.propellant_mass_kg)
    });
    let sc_does_tli = budget.as_ref().is_some_and(|b| b.maneuvers.contains(&Maneuver::Tli));
    let tli_engine = move |s: &TransferSolution| -> Option<launch_window_prediction::refine::Engine> {
        use launch_window_prediction::refine::Engine;
        if sc_does_tli {
            return sc_engine.clone().map(|(name, thrust_n, isp_s, mass_kg, propellant_kg)| Engine {
                name: format!("{name} (spacecraft)"),
                thrust_n,
                isp_s,
                mass_kg,
                propellant_kg,
            });
        }
        let st = ascent.get(&s.spaceport_id)?;
        let left = s.ascent.simulation.as_ref().map_or(st.upper.propellant_kg, |x| x.upper_stage_propellant_left_kg);
        Some(Engine {
            name: format!("upper stage of {}", st.rocket_name),
            thrust_n: st.upper.thrust_n,
            isp_s: st.upper.isp_s,
            mass_kg: st.upper.dry_mass_kg + left + st.payload_kg,
            propellant_kg: left,
        })
    };
    let mut reports = tokio::task::spawn_blocking(move || {
        use rayon::prelude::*;
        // Космодромы считаются параллельно (внутри каждого перебор тоже параллельный).
        let mut reports = ports
            .par_iter()
            .map(|p| trajectory.solve(&almanac, p, constraint, budget.as_ref(), target.as_ref()))
            .collect::<launch_window_prediction::Result<Vec<_>>>()?;
        if numerical.enabled {
            let fastest_goal =
                budget.as_ref().is_some_and(|b| b.objective == launch_window_prediction::trajectory::Objective::Fastest);
            // Уточняет одно решение: true — уточнение прошло (подошло оно или нет, см. `verdict`).
            let refine_solution = |s: &mut TransferSolution| -> launch_window_prediction::Result<bool> {
                let engine = tli_engine(s);
                match refine(&almanac, &numerical, s, target.as_ref(), budget.as_ref(), engine.as_ref())? {
                    Ok(refined) => {
                        s.refined = Some(refined);
                        Ok(true)
                    }
                    Err(reason) => {
                        s.refine_failure = Some(reason);
                        Ok(false)
                    }
                }
            };
            // Подходит ли уточнённое решение: условия посадки выполнены и касание в заданные сутки.
            let verdict = |s: &TransferSolution| -> Option<String> {
                let Some(n) = &s.refined else {
                    return Some(s.refine_failure.clone().unwrap_or_else(|| "not refined".into()));
                };
                if let Some(b) = &n.tli_burn
                    && b.exceeds_propellant
                    && !b.impulsive_exceeds_propellant
                {
                    return Some(format!(
                        "{} runs out of propellant on the finite trans-lunar injection burn ({:.0} m/s needed with {:.0} m/s gravity loss, {:.0} m/s available)",
                        b.engine.name,
                        b.delta_v_km_s * 1000.0,
                        b.gravity_loss_m_s,
                        b.max_delta_v_km_s * 1000.0
                    ));
                }
                let Some(l) = &n.landing else { return None };
                if !l.conditions_met {
                    return Some(l.unmet_reason.unwrap_or("landing conditions are not met").to_string());
                }
                if let DateConstraint::ArrivalDay(d) = constraint {
                    if l.touchdown_epoch < d || l.touchdown_epoch >= d + 1 * Unit::Day {
                        return Some("the refined touchdown falls outside the requested day".into());
                    }
                }
                None
            };
            let refine_one = |r: &mut TransferReport| -> launch_window_prediction::Result<bool> {
                let Some(s) = r.solution.as_mut() else { return Ok(false) };
                let ok = refine_solution(s)?;
                let Some(why) = verdict(s) else { return Ok(ok) };
                if target.is_none() {
                    return Ok(ok);
                }
                // Уточнение отвергло вариант (обычно медленный перелёт: Луна притягивает аппарат
                // раньше, чем считают формулы, и прилёт смещается). Проверяем другие моменты
                // прилёта из перебора — по тому же порядку, что и выбор лучшего.
                let mut attempts = vec![RefineAttempt {
                    arrival: s.arrival.epoch,
                    launch: Some(s.launch.epoch),
                    accepted: false,
                    reason: Some(why),
                }];
                let mut tried = vec![s.arrival.epoch];
                // Насколько уточнение сдвигает прилёт (у соседних вариантов примерно так же) и
                // сколько длится спуск от схода с орбиты до касания: по ним прогнозируется
                // ожидание на орбите у каждого варианта, и заведомо не влезающие в лимит пропускаются.
                let shift_of = |s: &TransferSolution| s.refined.as_ref().map(|n| (n.arrival_epoch - s.arrival.epoch).to_seconds());
                let mut shift_s = shift_of(s);
                let descent_s = s.landing.map_or(0.0, |l| (l.epoch - l.doi_epoch).to_seconds());
                let max_wait_s = target.as_ref().map_or(f64::MAX, |t| t.window.max_wait_in_orbit_h * 3600.0);
                let cands: Vec<ScanPoint> =
                    r.coarse_scan.iter().chain(&r.fine_scan).filter(|p| p.feasible).cloned().collect();
                let key = |p: &ScanPoint| {
                    if fastest_goal { p.time_of_flight_h.unwrap_or(f64::MAX) } else { p.delta_v_total_km_s.unwrap_or(f64::MAX) }
                };
                let spacing = 30 * Unit::Minute;
                while attempts.len() <= MAX_REFINE_ATTEMPTS {
                    let fits = |p: &ScanPoint, shift: Option<f64>| match (shift, p.landing) {
                        // Запас 1 ч на разброс сдвига между вариантами.
                        (Some(sh), Some(l)) => (l - p.arrival).to_seconds() - descent_s - sh <= max_wait_s - 3600.0,
                        _ => true,
                    };
                    let fresh = |p: &&ScanPoint| !tried.iter().any(|t| (p.arrival - *t).abs() < spacing);
                    let pick = cands
                        .iter()
                        .filter(fresh)
                        .filter(|p| fits(p, shift_s))
                        .min_by(|a, b| key(a).total_cmp(&key(b)))
                        .or_else(|| cands.iter().filter(fresh).min_by(|a, b| key(a).total_cmp(&key(b))));
                    let Some(p) = pick.cloned() else { break };
                    tried.push(p.arrival);
                    let Some(mut alt) =
                        trajectory.solve_at(&almanac, &r.spaceport, constraint, budget.as_ref(), target.as_ref(), p.arrival)?
                    else {
                        attempts.push(RefineAttempt {
                            arrival: p.arrival,
                            launch: p.launch,
                            accepted: false,
                            reason: Some("no transfer at this arrival".into()),
                        });
                        continue;
                    };
                    refine_solution(&mut alt)?;
                    shift_s = shift_of(&alt).or(shift_s);
                    let why = verdict(&alt);
                    attempts.push(RefineAttempt {
                        arrival: p.arrival,
                        launch: Some(alt.launch.epoch),
                        accepted: why.is_none(),
                        reason: why.clone(),
                    });
                    if why.is_none() {
                        r.solution = Some(alt);
                        r.refine_attempts = attempts;
                        return Ok(true);
                    }
                }
                // Ни один не подошёл: остаётся лучший вариант, с причиной в уточнении.
                r.refine_attempts = attempts;
                Ok(ok)
            };
            if numerical.refine_all {
                // Лучший вариант каждого космодрома уточняется численным интегрированием.
                reports.par_iter_mut().try_for_each(|r| refine_one(r).map(|_| ()))?;
            } else {
                // Уточняется только лучший вариант из всех космодромов; если уточнить его
                // не удалось — следующий по порядку.
                let fastest = budget.as_ref().is_some_and(|b| b.objective == launch_window_prediction::trajectory::Objective::Fastest);
                let mut order: Vec<usize> = (0..reports.len()).filter(|&i| reports[i].solution.is_some()).collect();
                order.sort_by(|&a, &b| ranking_key(&reports[a], fastest).total_cmp(&ranking_key(&reports[b], fastest)));
                let mut done = None;
                for &i in &order {
                    if refine_one(&mut reports[i])? {
                        done = Some(i);
                        break;
                    }
                }
                for &i in &order {
                    if Some(i) == done {
                        continue;
                    }
                    if let Some(s) = reports[i].solution.as_mut()
                        && s.refined.is_none()
                        && s.refine_failure.is_none()
                    {
                        s.refine_skipped = Some(
                            "not refined: only the best variant is refined numerically ([numerical] refine_all = false)".into(),
                        );
                    }
                }
            }
        }
        Ok::<_, launch_window_prediction::Error>(reports)
    })
    .await??;
    for r in &reports {
        log.event("transfer", r);
    }
    // Журнал траектории для 3D — отдельный файл.
    let traj_path = log.trajectory_path();
    let site = landing_site.map(|s| (s.site.lat_deg, s.site.lon_deg, s.site.height_km));
    match launch_window_prediction::trajlog::write(&traj_path, log.path(), &reports, site) {
        Ok(()) => {
            log.event("trajectory_log", json!({ "path": traj_path }));
            log.merge("files", json!({ "trajectory_log": file_name(&traj_path) }));
            println!("Trajectory log (3D): {}", traj_path.display());
        }
        Err(e) => {
            log.event("trajectory_log", json!({ "error": e.to_string() }));
            eprintln!("warning: could not write the trajectory log: {e}");
        }
    }
    let fastest = reports
        .first()
        .and_then(|r| r.budget.as_ref())
        .is_some_and(|b| b.objective == launch_window_prediction::trajectory::Objective::Fastest);
    reports.sort_by(|a, b| ranking_key(a, fastest).total_cmp(&ranking_key(b, fastest)));
    log.event(
        "transfer_ranking",
        reports
            .iter()
            .map(|r| json!({
                "spaceport": r.spaceport.id,
                "feasible": r.solution.is_some(),
                "delta_v_total_km_s": r.solution.as_ref().map(|s| s.delta_v_total_km_s),
            }))
            .collect::<Vec<_>>(),
    );
    Ok(reports)
}

/// Порядок вариантов: по времени перелёта (цель «быстрее») или по суммарному Δv.
fn ranking_key(r: &TransferReport, fastest: bool) -> f64 {
    r.solution.as_ref().map_or(f64::INFINITY, |s| {
        if fastest { s.departure.time_of_flight_s.to_seconds() } else { s.delta_v_total_km_s }
    })
}

/// Тип траектории для заголовков вывода.
fn trajectory_label(reports: &[TransferReport]) -> &'static str {
    match reports.first().and_then(|r| r.trajectory.parse::<TrajectoryKind>().ok()).unwrap_or_default() {
        TrajectoryKind::Direct => "direct transfer",
        TrajectoryKind::Hohmann => "Hohmann transfer",
    }
}

fn print_transfers(reports: &[TransferReport], constraint: DateConstraint) {
    match constraint {
        DateConstraint::LaunchDay(d) => println!("== Transfer ({}): best launch on {d} (UTC day)", trajectory_label(reports)),
        DateConstraint::ArrivalDay(d) => println!(
            "== Transfer ({}): best {} on {d} (UTC day)",
            trajectory_label(reports),
            if reports.iter().any(|r| r.landing_site.is_some()) { "landing at the site" } else { "arrival" }
        ),
    }
    println!(
        "{:<12} {:<24} {:<24} {:<24} {:>7} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6}",
        "spaceport", "launch (UTC)", "arrival (UTC)", "landing (UTC)", "tof_h", "incl", "az", "tli", "loi", "total", "rot"
    );
    for r in reports {
        match &r.solution {
            Some(s) => println!(
                "{:<12} {:<24} {:<24} {:<24} {:>7.1} {:>6.1} {:>6.1} {:>6.3} {:>6.3} {:>6.3} {:>6.3}",
                r.spaceport.id,
                ts(s.launch.epoch),
                ts(s.arrival.epoch),
                s.landing.map_or("-".to_string(), |l| ts(l.epoch)),
                hours(s.departure.time_of_flight_s),
                s.parking_orbit.inclination_deg,
                s.launch.azimuth_deg,
                s.departure.delta_v_magnitude_km_s,
                s.arrival.delta_v_loi_km_s,
                s.delta_v_total_km_s,
                s.launch.earth_rotation_gain_km_s,
            ),
            None => println!(
                "{:<12} no transfer this day: {}",
                r.spaceport.id,
                r.no_solution_reason.unwrap_or("no feasible transfer")
            ),
        }
    }
    if let Some((r, s)) = reports.first().and_then(|r| r.solution.as_ref().map(|s| (r, s))) {
        println!(
            "Best: {} — launch {}, parking orbit insertion {} (ascent {:.0} min, {:.0} km downrange), TLI {} (after {:.0} min in parking orbit), arrival {} \
             ({:.1} h of flight), Δv {:.3} km/s (TLI {:.3} + lunar orbit insertion {:.3}).",
            r.spaceport.name,
            ts(s.launch.epoch),
            ts(s.ascent.insertion_epoch),
            s.ascent.duration_s / 60.0,
            s.ascent.downrange_km,
            ts(s.departure.epoch),
            (s.departure.coast_s.to_seconds() - s.ascent.duration_s) / 60.0,
            ts(s.arrival.epoch),
            hours(s.departure.time_of_flight_s),
            s.delta_v_total_km_s,
            s.departure.delta_v_magnitude_km_s,
            s.arrival.delta_v_loi_km_s,
        );
        if let Some(h) = &s.hohmann {
            println!(
                "  Hohmann: half-ellipse {:.0} x {:.0} km (a = {:.0} km), tangential TLI to {:.4} km/s, apogee speed {:.4} km/s; \
                 nominal flight time {:.2} h, timing off by {:.1} s (tolerance {:.0} s), apogee {:+.1} km from the Moon's centre distance {:.0} km.",
                h.perigee_radius_km,
                h.apogee_radius_km,
                h.semi_major_axis_km,
                h.perigee_speed_km_s,
                h.apogee_speed_km_s,
                h.nominal_time_of_flight_s / 3600.0,
                h.timing_error_s,
                h.timing_tolerance_s,
                h.apogee_miss_km,
                h.moon_distance_km,
            );
        }
        let a = &s.ascent;
        match &a.simulation {
            Some(sim) => {
                let stage_events: Vec<String> = sim
                    .events
                    .iter()
                    .filter(|e| !matches!(e.name, "liftoff" | "pitch_kick" | "guidance" | "insertion"))
                    .map(|e| format!("{} T+{:.0}s {:.0} km", e.name, e.t_s, e.altitude_km))
                    .collect();
                println!(
                    "  Ascent (simulated, {}): liftoff mass {:.0} t, T/W {:.2}, pitch-over {:.1} deg; {}; orbit {:.0} x {:.0} km at T+{:.0}s, {:.0} km downrange.",
                    a.rocket.as_deref().unwrap_or("rocket"),
                    sim.liftoff_mass_kg / 1000.0,
                    sim.liftoff_thrust_to_weight,
                    sim.pitch_kick_deg,
                    stage_events.join(", "),
                    sim.perigee_km,
                    sim.apogee_km,
                    sim.duration_s,
                    sim.downrange_km,
                );
                println!(
                    "    max-Q {:.1} kPa at T+{:.0}s{}, max {:.1} g; losses: gravity {:.0}, drag {:.0}, steering {:.0}, back-pressure {:.0} m/s; upper stage {}: {:.0} kg propellant left = {:.0} m/s for TLI (needs {:.0}){}.",
                    sim.max_q_pa / 1000.0,
                    sim.max_q_t_s,
                    if sim.max_q_above_limit { " (above the rocket's limit)" } else { "" },
                    sim.max_acceleration_g,
                    sim.gravity_loss_m_s,
                    sim.drag_loss_m_s,
                    sim.steering_loss_m_s,
                    sim.pressure_loss_m_s,
                    if sim.upper_stage_used { "helped reach orbit" } else { "not used in ascent" },
                    sim.upper_stage_propellant_left_kg,
                    sim.upper_stage_delta_v_left_m_s,
                    s.refined.as_ref().map_or(s.departure.delta_v_magnitude_km_s, |n| n.tli_delta_v_km_s) * 1000.0,
                    if sim.upper_stage_delta_v_left_m_s < s.refined.as_ref().map_or(s.departure.delta_v_magnitude_km_s, |n| n.tli_delta_v_km_s) * 1000.0 {
                        " — NOT ENOUGH"
                    } else {
                        ""
                    },
                );
            }
            None => println!(
                "  Ascent (simplified): {:.0} min, {:.0} km downrange — {}.",
                a.duration_s / 60.0,
                a.downrange_km,
                a.model_note.as_deref().unwrap_or("")
            ),
        }
        println!(
            "  Moon at arrival: {:.0} km, declination {:.1} deg; transfer angle {:.1} deg; v_inf {:.3} km/s; C3 {:.2} km2/s2.",
            s.arrival.moon_distance_km,
            s.arrival.moon_declination_deg,
            s.departure.transfer_angle_deg,
            s.arrival.v_infinity_magnitude_km_s,
            s.departure.c3_km2_s2,
        );
        if let Some(l) = &s.landing {
            let d = &l.descent;
            println!(
                "  Landing at the site (touchdown): {} — orbit inclination {:.1} deg to the lunar equator, plane turned {:+.2} deg by lunar J2 since LOI;",
                ts(l.epoch),
                l.lunar_orbit_inclination_deg,
                l.nodal_precession_deg,
            );
            println!(
                "    descent: DOI {} ({:.1} h in orbit, pass {}, {:.0} m/s) -> coast {:.0} min -> PDI {} ({:.2} km/s) -> braking {:.1} min over {:.0} km ({:.0} m/s at {:.0}% thrust, {} engine) -> final approach {:.1} min (+{:.0} m/s).",
                ts(l.doi_epoch),
                l.wait_after_arrival_s / 3600.0,
                l.orbit_pass,
                d.doi_delta_v_m_s,
                d.coast_s / 60.0,
                ts(l.pdi_epoch),
                d.pdi_speed_km_s,
                d.braking_duration_s / 60.0,
                d.powered_descent_range_km,
                d.powered_descent_delta_v_m_s,
                d.throttle * 100.0,
                d.vehicle_source,
                (d.powered_descent_duration_s - d.braking_duration_s) / 60.0,
                d.terminal_allowance_m_s,
            );
            println!(
                "    at landing: Sun {:.1} deg (terrain {:.1} deg, {}), Earth {:.1} deg (terrain {:.1} deg, {}).",
                l.sun.elevation_deg,
                l.sun.horizon_deg,
                if l.sun.visible { "visible" } else { "hidden" },
                l.earth.elevation_deg,
                l.earth.horizon_deg,
                if l.earth.visible { "visible" } else { "hidden" },
            );
        }
        if let Some(b) = &s.budget {
            println!(
                "  Spacecraft {}: needs {:.0} m/s ({}) + {:.0}% reserve = {:.0} m/s of {:.0} m/s available, margin {:.0} m/s{}.",
                b.spacecraft_name,
                b.required_m_s,
                [
                    b.tli_m_s.map(|v| format!("TLI {v:.0}")),
                    b.loi_m_s.map(|v| format!("LOI {v:.0}")),
                    b.descent_m_s.map(|v| format!("descent {v:.0}")),
                ]
                .into_iter()
                .flatten()
                .collect::<Vec<_>>()
                .join(" + "),
                b.reserve_percent,
                b.required_with_reserve_m_s,
                b.available_m_s,
                b.margin_m_s,
                if b.time_of_flight_at_lower_limit {
                    "; flight time is at the tof_min_days limit (faster is not allowed by transfer.toml)"
                } else {
                    ""
                },
            );
        }
        if let Some(n) = &s.refined {
            let c = &n.config;
            println!(
                "  REFINED by numerical integration (Earth{}, Moon gravity degree {}{}): arrival {} ({:+.1} min vs formulas), TLI {:.3} km/s ({:+.1} m/s), LOI {:.3} km/s ({:+.1} m/s); parking-orbit point moved {:.0} km; {} iterations{}.",
                if c.earth_j2 { " with J2" } else { "" },
                c.lunar_gravity_degree,
                if c.sun_gravity { ", Sun" } else { "" },
                ts(n.arrival_epoch),
                n.arrival_shift_s / 60.0,
                n.tli_delta_v_km_s,
                n.tli_delta_v_change_m_s,
                n.loi_delta_v_km_s,
                n.loi_delta_v_change_m_s,
                n.parking_position_shift_km,
                n.iterations,
                if n.converged { "" } else { ", NOT CONVERGED" },
            );
            match &n.tli_burn {
                Some(b) => println!(
                    "    TLI burn ({}, {:.0} kN, Isp {:.0} s): {} -> {} ({:.1} min, {:.0} deg along the orbit), engine delta-v {:.0} m/s = impulsive {:.0} + gravity loss {:.0} m/s; propellant {:.0} of {:.0} kg{}.",
                    b.engine.name,
                    b.engine.thrust_n / 1000.0,
                    b.engine.isp_s,
                    ts(b.start_epoch),
                    ts(b.end_epoch),
                    b.duration_s / 60.0,
                    b.burn_arc_deg,
                    b.delta_v_km_s * 1000.0,
                    b.impulsive_delta_v_km_s * 1000.0,
                    b.gravity_loss_m_s,
                    b.propellant_used_kg,
                    b.engine.propellant_kg,
                    if b.exceeds_propellant { " — NOT ENOUGH PROPELLANT" } else { "" },
                ),
                None => println!("    TLI burn: impulsive (no engine data for a finite burn)."),
            }
            if let Some(l) = &n.landing {
                println!(
                    "    lunar orbit {:.1}..{:.1} km while waiting (mascons); DOI {} ({:.0} m/s{}) -> PDI {} at {:.1} km ({:.1} km above the site) -> braking {:.1} min over {:.0} km ({:.0} m/s) -> touchdown {} ({:+.1} min vs formulas), miss {:.2} km.",
                    l.orbit_min_altitude_km,
                    l.orbit_max_altitude_km,
                    ts(l.doi_epoch),
                    l.doi_delta_v_m_s,
                    if l.extra_orbits > 0 { format!(", {} orbit(s) later than by formulas", l.extra_orbits) } else { String::new() },
                    ts(l.pdi_epoch),
                    l.pdi_altitude_km,
                    l.pdi_height_above_site_km,
                    l.braking.duration_s / 60.0,
                    l.braking.range_km,
                    l.braking.delta_v_m_s,
                    ts(l.touchdown_epoch),
                    l.touchdown_shift_s / 60.0,
                    l.miss_km,
                );
                println!(
                    "    at touchdown: Sun {:.1} deg ({}), Earth {:.1} deg ({}) — landing conditions {}.",
                    l.sun.elevation_deg,
                    if l.sun.visible { "visible" } else { "hidden" },
                    l.earth.elevation_deg,
                    if l.earth.visible { "visible" } else { "hidden" },
                    match l.unmet_reason {
                        None => "met".to_string(),
                        Some(reason) => format!("NOT MET: {reason}"),
                    },
                );
            }
            if let Some(b) = &n.budget {
                println!(
                    "    spacecraft budget (refined): needs {:.0} m/s ({}) + {:.0}% reserve = {:.0} of {:.0} m/s, margin {:.0} m/s{}.",
                    b.required_m_s,
                    [
                        b.tli_m_s.map(|v| format!("TLI {v:.0}")),
                        b.loi_m_s.map(|v| format!("LOI {v:.0}")),
                        b.descent_m_s.map(|v| format!("descent {v:.0}")),
                    ]
                    .into_iter()
                    .flatten()
                    .collect::<Vec<_>>()
                    .join(" + "),
                    b.reserve_percent,
                    b.required_with_reserve_m_s,
                    b.available_m_s,
                    b.margin_m_s,
                    if b.margin_m_s < 0.0 { " — NOT ENOUGH" } else { "" },
                );
            }
        }
        if let Some(reason) = &s.refine_failure {
            println!("  Numerical refinement failed: {reason}. The figures above are from formulas only.");
        }
        if s.refine_skipped.is_some() {
            println!("  Not refined numerically (only the best variant is; set [numerical] refine_all = true for all). Figures are from formulas.");
        }
        if !r.refine_attempts.is_empty() {
            let accepted = r.refine_attempts.iter().any(|a| a.accepted);
            println!(
                "  Variants checked by numerical refinement ({}):",
                if accepted {
                    "the best formula variant did not hold up, the first one that did is shown"
                } else {
                    "none met the landing conditions; the best formula variant is shown"
                }
            );
            for a in &r.refine_attempts {
                println!(
                    "    arrival {} (launch {}): {}",
                    ts(a.arrival),
                    a.launch.map_or_else(|| "-".into(), ts),
                    if a.accepted { "OK".to_string() } else { a.reason.clone().unwrap_or_default() },
                );
            }
        }
        if let Some(src) = &r.spaceport.limits_source {
            println!("  Spaceport limits applied: {src}");
        } else {
            println!("  Spaceport limits: none known (any eastward launch allowed).");
        }
    }
    println!(
        "Δv in km/s from the parking orbit (tli) and into lunar orbit (loi); rot = speed gained from \
         Earth's rotation at launch (reference only, not part of Δv)."
    );
    let kind = reports.first().and_then(|r| r.trajectory.parse::<TrajectoryKind>().ok()).unwrap_or_default();
    println!("Model: {}.", kind.description());
}

fn print_visibility(
    site: &SiteVisibility,
    flat_height_m: f64,
    vis: &[VisibilitySample],
    summary: &VisibilitySummary,
) {
    match &site.terrain {
        Some(t) => {
            println!(
                "Terrain: map {} ({:.0} m/pix here), surface height {:.1} m",
                t.map_name, t.resolution_m, t.surface_height_m,
            );
            if t.truncated_rays > 0 {
                println!(
                    "WARNING: {} horizon rays left the maps before {} km; far terrain there is not counted.",
                    t.truncated_rays, site.config.horizon.max_range_km
                );
            }
        }
        None => println!("Terrain: not used (flat horizon), surface height {flat_height_m:.1} m"),
    }
    println!(
        "Visibility rules: observer +{:.1} m, sun {:?} +{:.1} deg, earth {:?} +{:.1} deg",
        site.config.horizon.observer_height_m,
        site.config.sun.criterion,
        site.config.sun.min_clearance_deg,
        site.config.earth.criterion,
        site.config.earth.min_clearance_deg,
    );
    println!(
        "{:<24} {:>8} {:>7} {:>7} {:>5} {:>4}   {:>8} {:>7} {:>7} {:>5}",
        "epoch", "sun_az", "sun_el", "horiz", "lit%", "sun", "earth_az", "earth_el", "horiz", "earth"
    );
    for v in vis {
        println!(
            "{:<24} {:>8.3} {:>7.3} {:>7.3} {:>5.0} {:>4}   {:>8.3} {:>7.3} {:>7.3} {:>5}",
            format!("{}", v.epoch),
            v.sun.azimuth_deg,
            v.sun.elevation_deg,
            v.sun.horizon_deg,
            v.sun.visible_fraction * 100.0,
            if v.sun.visible { "up" } else { "down" },
            v.earth.azimuth_deg,
            v.earth.elevation_deg,
            v.earth.horizon_deg,
            if v.earth.visible { "up" } else { "down" },
        );
    }
    println!(
        "Summary over {} samples: sun {:.1}%, earth {:.1}%, both {:.1}%",
        summary.samples, summary.sun_pct, summary.earth_pct, summary.both_pct
    );
    println!(
        "Longest gaps: sun {:.1} h, earth {:.1} h, both {:.1} h (a gap at the window edge may continue beyond it)",
        hours(summary.longest_sun_gap),
        hours(summary.longest_earth_gap),
        hours(summary.longest_both_gap),
    );
}

fn print_spacecraft(sc: &ResolvedSpacecraft) {
    let c = &sc.config;
    let mark = |f: &str| if sc.overridden.contains(&f) { "*" } else { "" };
    let opt = |v: Option<f64>| v.map_or("-".to_string(), |v| format!("{v}"));
    println!("== Spacecraft: {} [{}]", c.name, sc.status_label);
    println!(
        "  base: dry {}{} kg, propellant {}{} kg, payload {}{} kg, thrust {}{} N, Isp {}{} s",
        c.dry_mass_kg, mark("dry_mass_kg"),
        c.propellant_mass_kg, mark("propellant_mass_kg"),
        c.payload_mass_kg, mark("payload_mass_kg"),
        c.max_thrust_n, mark("max_thrust_n"),
        c.isp_s, mark("isp_s"),
    );
    println!(
        "  advanced: min throttle {}{} %, solar {}{} m2 x {}{}, battery {}{} Wh, idle {}{} W, tx {}{} W",
        opt(c.min_throttle_percent), mark("min_throttle_percent"),
        opt(c.solar_panel_area_m2), mark("solar_panel_area_m2"),
        opt(c.solar_panel_efficiency), mark("solar_panel_efficiency"),
        opt(c.battery_capacity_wh), mark("battery_capacity_wh"),
        opt(c.idle_power_consumption_w), mark("idle_power_consumption_w"),
        opt(c.tx_power_w), mark("tx_power_w"),
    );
    let d = &sc.derived;
    println!(
        "  derived: m0 {:.0} kg, delta-v {:.0} m/s, lunar TWR {:.2}, eclipse autonomy {}",
        d.initial_mass_kg,
        d.delta_v_m_s,
        d.twr_moon,
        d.eclipse_autonomy_h.map_or("unknown".to_string(), |h| format!("{h:.1} h")),
    );
    if !sc.overridden.is_empty() {
        println!("  overridden (*): {}", sc.overridden.join(", "));
    }
    for i in &sc.issues {
        let tag = match i.severity {
            Severity::Error => "ERROR",
            Severity::Warning => "WARNING",
        };
        println!("  {tag}: {}", i.message);
    }
}

fn print_spacecraft_list(library: &SpacecraftLibrary, log: &RunLog) {
    println!("== Lunar spacecraft");
    println!("{:<26} {:<34} {:>9} {:>6}  {}", "id", "status", "dv_m_s", "twr", "issues");
    let mut rows = Vec::new();
    for e in library.entries() {
        match library.resolve(&e.id, &SpacecraftParams::default()) {
            Ok(r) => {
                let issues = r.issues.iter().map(|i| i.message.as_str()).collect::<Vec<_>>().join("; ");
                println!(
                    "{:<26} {:<34} {:>9.0} {:>6.2}  {}",
                    e.id, r.status_label, r.derived.delta_v_m_s, r.derived.twr_moon,
                    if issues.is_empty() { "-" } else { &issues }
                );
                rows.push(json!({ "id": e.id, "entry": e, "resolved": r }));
            }
            Err(err) => {
                println!("{:<26} cannot resolve: {err}", e.id);
                rows.push(json!({ "id": e.id, "entry": e, "error": err.to_string() }));
            }
        }
    }
    log.event("spacecraft_library", rows);
}

fn print_rocket(r: &ResolvedLaunchVehicle, cap: &LaunchCapability, payload_is_spacecraft: bool) {
    let c = &r.config;
    let mark = |f: &str| if r.overridden.contains(&f) { "*" } else { "" };
    println!("== Rocket: {} [{}]", c.name, r.status_label);
    println!(
        "  to the Moon (TLI) {}{} kg, spaceports: {}",
        c.payload_tli_kg,
        mark("payload_tli_kg"),
        if c.any_spaceport { "any (custom rocket)".to_string() } else { c.spaceports.join(", ") },
    );
    println!(
        "  upper stage {}: dry {}{} kg, propellant {}{} kg, thrust {}{} N, Isp {}{} s",
        c.upper_stage_name.as_deref().unwrap_or("-"),
        c.upper_stage_dry_mass_kg, mark("upper_stage_dry_mass_kg"),
        c.upper_stage_propellant_kg, mark("upper_stage_propellant_kg"),
        c.upper_stage_thrust_n, mark("upper_stage_thrust_n"),
        c.upper_stage_isp_s, mark("upper_stage_isp_s"),
    );
    let what = if payload_is_spacecraft { "spacecraft initial mass" } else { "rated TLI payload" };
    if cap.within_rated_payload {
        println!(
            "  payload {:.0} kg ({what}): TLI delta-v up to {:.0} m/s (rated payload corresponds to {:.0} m/s)",
            cap.payload_kg, cap.max_tli_delta_v_m_s, cap.reference_tli_delta_v_m_s,
        );
    } else {
        println!(
            "  payload {:.0} kg ({what}) is HEAVIER than the rated {:.0} kg to the Moon",
            cap.payload_kg, cap.payload_tli_kg
        );
    }
    if !r.overridden.is_empty() {
        println!("  overridden (*): {}", r.overridden.join(", "));
    }
    for e in &r.errors {
        println!("  ERROR: {e}");
    }
}

fn print_rocket_list(rockets: &LaunchVehicleLibrary, log: &RunLog) {
    println!("== Rockets");
    println!("{:<26} {:<34} {:>9}  {}", "id", "status", "tli_kg", "spaceports");
    let mut rows = Vec::new();
    for e in rockets.entries() {
        match rockets.resolve(&e.id, &LaunchVehicleParams::default()) {
            Ok(r) => {
                println!(
                    "{:<26} {:<34} {:>9.0}  {}",
                    e.id,
                    r.status_label,
                    r.config.payload_tli_kg,
                    if r.config.any_spaceport { "any (custom)".to_string() } else { r.config.spaceports.join(", ") }
                );
                rows.push(json!({ "id": e.id, "entry": e, "resolved": r }));
            }
            Err(err) => {
                println!("{:<26} cannot resolve: {err}", e.id);
                rows.push(json!({ "id": e.id, "entry": e, "error": err.to_string() }));
            }
        }
    }
    log.event("launch_vehicle_library", rows);
}
