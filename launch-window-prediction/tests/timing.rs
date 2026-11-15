//! Замер времени работы программы на разных вводных данных.
//!
//! Запускает собранную программу как обычный пользователь (через командную строку) на наборе
//! сценариев: с точкой посадки и без, по дате старта и по дате посадки, один космодром и все,
//! короткий и длинный перебор, с численным уточнением и без. Каждый сценарий запускается
//! несколько раз; печатается таблица (минимум, медиана, максимум) и итог.
//!
//! Нужны ядра в `data/kernels` и рельеф в `data/terrain`, поэтому тест помечен `#[ignore]`.
//! Запуск (время имеет смысл только в release):
//!
//! ```text
//! cargo test --release --test timing -- --ignored --nocapture
//! ```
//!
//! Моделирование подъёма ракеты кэшируется на диске (папка `--cache-dir`), поэтому у каждого
//! сценария первый запуск — без кэша («макс»), следующие — с кэшем («мин»). Кэш — во временной
//! папке теста; при запуске теста он пустой.
//!
//! Переменные окружения: `LWP_TIMING_RUNS` — число запусков каждого сценария (по умолчанию 3),
//! `LWP_TIMING_LIMIT_S` — предел для самого долгого сценария, с (по умолчанию 120).

use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Duration, Instant};

const BIN: &str = env!("CARGO_BIN_EXE_launch-window-prediction");

struct Scenario {
    name: &'static str,
    args: Vec<String>,
    /// Ключи `[секция] ключ = значение` поверх `data/transfer.toml` (пусто — файл как есть).
    overrides: &'static [(&'static str, &'static str, &'static str)],
}

fn args(s: &str) -> Vec<String> {
    s.split_whitespace().map(str::to_string).collect()
}

fn scenarios() -> Vec<Scenario> {
    let common = "--kernels-dir data/kernels --terrain-dir data/terrain";
    let sc = "--rocket sls_block1 --spacecraft apollo_lm_descent --sc payload_mass_kg=0";
    // Без ракеты: SLS летает только с мыса Канаверал, а здесь нужны все космодромы.
    let sc_any = "--spacecraft apollo_lm_descent --sc payload_mass_kg=0";
    let off: &[_] = &[("numerical", "enabled", "false")];
    let smooth: &[_] = &[("numerical", "lunar_gravity_degree", "0")];
    vec![
        Scenario {
            name: "1 порт, дата старта, без точки",
            args: args(&format!("{common} --launch-date 2026-12-10 --spaceport ksc")),
            overrides: &[],
        },
        Scenario {
            name: "все 11 портов, дата старта, без точки",
            args: args(&format!("{common} --launch-date 2026-12-10")),
            overrides: &[],
        },
        Scenario {
            name: "1 порт, точка, дата посадки (Новый год)",
            args: args(&format!("{common} --arrival-date 2027-01-01 --spaceport ksc {sc} --lat 5 --lon -10")),
            overrides: &[],
        },
        Scenario {
            name: "  то же, без численного уточнения",
            args: args(&format!("{common} --arrival-date 2027-01-01 --spaceport ksc {sc} --lat 5 --lon -10")),
            overrides: off,
        },
        Scenario {
            name: "  то же, Луна без масконов (степень 0)",
            args: args(&format!("{common} --arrival-date 2027-01-01 --spaceport ksc {sc} --lat 5 --lon -10")),
            overrides: smooth,
        },
        Scenario {
            name: "1 порт, точка, дата старта (Море Спокойствия)",
            args: args(&format!("{common} --launch-date 2026-12-10 --spaceport ksc {sc} --lat 0.674 --lon 23.473")),
            overrides: &[],
        },
        Scenario {
            name: "1 порт, южный полюс (нет посадки)",
            args: args(&format!("{common} --launch-date 2026-12-05 --spaceport ksc {sc} --lat -89.484 --lon 29.38")),
            overrides: &[],
        },
        Scenario {
            name: "все 11 портов, точка, дата посадки",
            args: args(&format!("{common} --arrival-date 2027-01-01 {sc_any} --lat 5 --lon -10")),
            overrides: &[],
        },
        Scenario {
            name: "1 порт, точка + видимость 72 ч шаг 6 ч",
            args: args(&format!(
                "{common} --arrival-date 2027-01-01 --spaceport ksc {sc} --lat 5 --lon -10 --hours 72 --step-min 360"
            )),
            overrides: &[],
        },
        Scenario {
            name: "1 порт, точка + видимость 30 сут шаг 1 ч",
            args: args(&format!(
                "{common} --arrival-date 2027-01-01 --spaceport ksc {sc} --lat 5 --lon -10 --hours 720 --step-min 60"
            )),
            overrides: &[],
        },
        Scenario {
            name: "все 11 портов, точка, без численного уточнения",
            args: args(&format!("{common} --arrival-date 2027-01-01 {sc_any} --lat 5 --lon -10")),
            overrides: off,
        },
    ]
}

/// Копия `data/transfer.toml` с заменёнными ключами.
fn transfer_config(dir: &Path, name: &str, overrides: &[(&str, &str, &str)]) -> PathBuf {
    let mut table: toml::Table = std::fs::read_to_string("data/transfer.toml").unwrap().parse().unwrap();
    for (section, key, value) in overrides {
        let v: toml::Value = format!("x = {value}").parse::<toml::Table>().unwrap()["x"].clone();
        table
            .entry(section.to_string())
            .or_insert_with(|| toml::Value::Table(Default::default()))
            .as_table_mut()
            .unwrap()
            .insert(key.to_string(), v);
    }
    let path = dir.join(format!("{name}.toml"));
    std::fs::write(&path, toml::to_string(&table).unwrap()).unwrap();
    path
}

fn run(s: &Scenario, config: Option<&Path>, log_dir: &Path, cache_dir: &Path) -> (Duration, String) {
    let mut cmd = Command::new(BIN);
    cmd.args(&s.args).arg("--log-dir").arg(log_dir).arg("--cache-dir").arg(cache_dir);
    if let Some(c) = config {
        cmd.arg("--transfer-config").arg(c);
    }
    let t = Instant::now();
    let out = cmd.output().expect("run the program");
    let elapsed = t.elapsed();
    let stdout = String::from_utf8_lossy(&out.stdout).to_string();
    assert!(
        out.status.success(),
        "{}: exit {:?}\n{}\n{}",
        s.name,
        out.status.code(),
        stdout,
        String::from_utf8_lossy(&out.stderr)
    );
    (elapsed, stdout)
}

fn secs(d: Duration) -> f64 {
    d.as_secs_f64()
}

#[test]
#[ignore = "needs data/kernels and data/terrain; run with --release"]
fn program_run_time() {
    if cfg!(debug_assertions) {
        eprintln!("ВНИМАНИЕ: отладочная сборка — время будет в 10-50 раз больше. Запускайте с --release.");
    }
    let runs: usize = std::env::var("LWP_TIMING_RUNS").ok().and_then(|v| v.parse().ok()).unwrap_or(3).max(1);
    let limit_s: f64 = std::env::var("LWP_TIMING_LIMIT_S").ok().and_then(|v| v.parse().ok()).unwrap_or(120.0);
    let tmp = std::env::temp_dir().join(format!("lwp-timing-{}", std::process::id()));
    std::fs::create_dir_all(&tmp).unwrap();
    let log_dir = tmp.join("logs");
    let cache_dir = tmp.join("cache");

    // Прогрев: первый запуск читает ядра и карты с диска, дальше они в кэше системы.
    // Кэш подъёма после прогрева очищается, чтобы сценарии начинали без него.
    let list = scenarios();
    let (cold, _) = run(&list[0], None, &log_dir, &cache_dir);
    let _ = std::fs::remove_dir_all(&cache_dir);

    println!();
    println!("Замер времени работы программы ({runs} запуск(а) на сценарий, {} ядер процессора)", num_cpus());
    println!("Первый (холодный) запуск: {:.2} с — чтение ядер и карт с диска", secs(cold));
    println!();
    println!("(макс — первый запуск сценария, без кэша подъёма ракеты; мин — с кэшем)");
    println!("{:<46} {:>8} {:>8} {:>8}   результат", "сценарий", "мин, с", "медиана", "макс, с");
    println!("{}", "-".repeat(110));
    let mut total = Duration::ZERO;
    let mut worst = (0.0f64, "");
    for (i, s) in list.iter().enumerate() {
        let config = (!s.overrides.is_empty()).then(|| transfer_config(&tmp, &format!("cfg{i}"), s.overrides));
        let mut times = Vec::with_capacity(runs);
        let mut last = String::new();
        for _ in 0..runs {
            let (t, out) = run(s, config.as_deref(), &log_dir, &cache_dir);
            total += t;
            times.push(t);
            last = out;
        }
        times.sort();
        let (min, med, max) = (secs(times[0]), secs(times[runs / 2]), secs(times[runs - 1]));
        if med > worst.0 {
            worst = (med, s.name);
        }
        // Строки таблицы космодромов: «id  дата-старта ...» (с вариантом) или «id  no transfer ...».
        let ports_ok = last
            .lines()
            .filter(|l| l.split_whitespace().nth(1).is_some_and(|w| w.starts_with("20") && w.contains('T')))
            .count();
        let ports_none = last.lines().filter(|l| l.contains("no transfer")).count();
        let refined = last.matches("REFINED").count();
        let failed = last.matches("refinement failed").count();
        let summary = format!(
            "портов с перелётом {ports_ok}/{}, уточнено {refined}{}",
            ports_ok + ports_none,
            if failed > 0 { format!(", не уточнено {failed}") } else { String::new() }
        );
        println!("{:<46} {:>8.2} {:>8.2} {:>8.2}   {summary}", s.name, min, med, max);
    }
    println!("{}", "-".repeat(110));
    println!("Всего {} запусков за {:.1} с; самый долгий сценарий: {} — {:.2} с", list.len() * runs, secs(total), worst.1, worst.0);
    println!();
    let _ = std::fs::remove_dir_all(&tmp);
    assert!(worst.0 < limit_s, "сценарий «{}» работает {:.1} с — дольше предела {limit_s} с", worst.1, worst.0);
}

fn num_cpus() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
}
