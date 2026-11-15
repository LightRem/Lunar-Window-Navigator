//! Замер времени расчёта в трёх режимах: быстрый (`fast`), экономный (`economy`) и окно старта
//! (`--launch-window`).
//!
//! Каждая группа сценариев считает одну и ту же задачу во всех режимах, после чего печатается
//! таблица: время, во сколько раз режим медленнее быстрого и что он нашёл. Тест проверяет:
//! 1. каждый режим укладывается в абсолютный предел (`LWP_TIMING_LIMIT_S`);
//! 2. экономный режим и окно старта не медленнее быстрого больше чем в `LWP_TIMING_RATIO` раз
//!    (экономный просматривает 2N+1 суток, перебирает наклонения и допускает долгий перелёт и
//!    ожидание у Луны, поэтому он ожидаемо в 5-25 раз медленнее; предел с запасом);
//! 3. режим сделал то, что от него ждут (в выводе есть блок экономного режима или окно старта).
//!
//! Нужны ядра в `data/kernels`; рельеф `data/terrain` используется, если папка есть. Тест помечен
//! `#[ignore]`; время имеет смысл только в release:
//!
//! ```text
//! cargo test --release --test timing_modes -- --ignored --nocapture
//! ```
//!
//! Переменные окружения: `LWP_TIMING_RUNS` — запусков каждого сценария (по умолчанию 2; первый —
//! без кэша подъёма ракеты, остальные с ним), `LWP_TIMING_LIMIT_S` — предел одного режима, с
//! (по умолчанию 300), `LWP_TIMING_RATIO` — наибольшее отношение к быстрому режиму (по умолчанию 40).

use std::path::Path;
use std::process::Command;
use std::time::{Duration, Instant};

const BIN: &str = env!("CARGO_BIN_EXE_launch-window-prediction");

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Mode {
    Fast,
    Economy,
    Window,
}

impl Mode {
    const ALL: [Mode; 3] = [Mode::Fast, Mode::Economy, Mode::Window];

    fn name(self) -> &'static str {
        match self {
            Mode::Fast => "fast",
            Mode::Economy => "economy",
            Mode::Window => "launch-window",
        }
    }

    fn flags(self) -> &'static [&'static str] {
        match self {
            Mode::Fast => &[],
            Mode::Economy => &["--mode", "economy"],
            Mode::Window => &["--launch-window"],
        }
    }
}

/// Одна задача, которая считается во всех режимах.
struct Group {
    name: &'static str,
    /// Общие аргументы (без режима).
    args: String,
    /// Режимы, которые применимы к задаче (окно старта — только с датой посадки).
    modes: &'static [Mode],
}

fn groups() -> Vec<Group> {
    let ksc_rocket = "--rocket sls_block1 --spacecraft apollo_lm_descent --sc payload_mass_kg=0";
    let any_port = "--spacecraft apollo_lm_descent --sc payload_mass_kg=0";
    let site = "--lat 0.674 --lon 23.473";
    vec![
        Group {
            name: "KSC, прямой, посадка в точку (дата посадки)",
            args: format!("--arrival-date 2027-02-13 --spaceport ksc {ksc_rocket} {site}"),
            modes: &Mode::ALL,
        },
        Group {
            name: "KSC, прямой, без точки (дата прилёта)",
            args: "--arrival-date 2027-02-13 --spaceport ksc".to_string(),
            modes: &Mode::ALL,
        },
        Group {
            name: "KSC, Гоман, посадка в точку (дата посадки)",
            args: format!("--arrival-date 2027-02-13 --spaceport ksc --trajectory hohmann {ksc_rocket} {site}"),
            modes: &Mode::ALL,
        },
        Group {
            name: "KSC, прямой, дата старта (окно старта не ищется)",
            args: "--launch-date 2027-02-06 --spaceport ksc".to_string(),
            modes: &[Mode::Fast, Mode::Economy],
        },
        Group {
            name: "все 11 портов, прямой, посадка в точку",
            args: format!("--arrival-date 2027-02-13 {any_port} {site}"),
            modes: &Mode::ALL,
        },
    ]
}

/// Что режим обязан показать в выводе. `None` — всё в порядке, иначе описание проблемы.
fn check_output(mode: Mode, out: &str) -> Option<String> {
    let has = |s: &str| out.contains(s);
    match mode {
        Mode::Fast if has("== Mode: fast") => None,
        Mode::Fast => Some("нет строки «== Mode: fast»".into()),
        Mode::Economy if has("== Mode: economy") && has("Economy: chosen day per spaceport") => None,
        Mode::Economy => Some("нет блока экономного режима".into()),
        Mode::Window if has("== Launch window") && has("Launch window per spaceport") && has("opens ") => None,
        Mode::Window if has("no launch window") => None, // окна нет — допустимый ответ, время всё равно замеряется
        Mode::Window => Some("нет блока окна старта".into()),
    }
}

fn run(group: &Group, mode: Mode, log_dir: &Path, cache_dir: &Path) -> (Duration, String) {
    let mut cmd = Command::new(BIN);
    cmd.args(["--kernels-dir", "data/kernels"]);
    if Path::new("data/terrain").is_dir() {
        cmd.args(["--terrain-dir", "data/terrain"]);
    }
    cmd.args(group.args.split_whitespace())
        .args(mode.flags())
        .arg("--log-dir")
        .arg(log_dir)
        .arg("--cache-dir")
        .arg(cache_dir);
    let t = Instant::now();
    let out = cmd.output().expect("run the program");
    let elapsed = t.elapsed();
    let stdout = String::from_utf8_lossy(&out.stdout).to_string();
    assert!(
        out.status.success(),
        "{} / {}: exit {:?}\n{}\n{}",
        group.name,
        mode.name(),
        out.status.code(),
        stdout,
        String::from_utf8_lossy(&out.stderr)
    );
    (elapsed, stdout)
}

/// Короткий итог из вывода: сколько космодромов с вариантом и (для окна) его размер.
fn result_note(mode: Mode, out: &str) -> String {
    // Основная таблица: от строки заголовка до первого следующего блока.
    let rows: Vec<&str> = out
        .lines()
        .skip_while(|l| !l.starts_with("spaceport "))
        .skip(1)
        .take_while(|l| !(l.starts_with("Economy:") || l.starts_with("Launch window per") || l.starts_with("Best:") || l.starts_with("Δv")))
        .collect();
    let ok = rows.iter().filter(|l| l.split_whitespace().nth(1).is_some_and(|w| w.starts_with("20") && w.contains('T'))).count();
    let mut note = format!("портов с вариантом {ok}/{}", rows.len());
    if mode == Mode::Window {
        // «  ksc  opens … closes …  (49.8 h from first to last, 6 interval(s), 11 min in total); least …»
        if let Some(l) = out.lines().find(|l| l.contains("interval(s)")) {
            if let (Some(i), Some(j)) = (l.find('('), l.find("); ")) {
                note.push_str(&format!(", окно {}", &l[i + 1..j]));
            }
        }
    }
    if mode == Mode::Economy {
        let shifted = out.lines().filter(|l| l.contains("saves ")).count();
        note.push_str(&format!(", со сдвигом даты ради экономии: {shifted}"));
    }
    note
}

fn secs(d: Duration) -> f64 {
    d.as_secs_f64()
}

#[test]
#[ignore = "needs data/kernels; run with --release"]
fn mode_run_time() {
    if cfg!(debug_assertions) {
        eprintln!("ВНИМАНИЕ: отладочная сборка — время будет в 3-20 раз больше. Запускайте с --release.");
    }
    let env_f = |k: &str, d: f64| std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d);
    let runs = env_f("LWP_TIMING_RUNS", 2.0).max(1.0) as usize;
    let limit_s = env_f("LWP_TIMING_LIMIT_S", 300.0);
    let ratio_limit = env_f("LWP_TIMING_RATIO", 40.0);

    let tmp = std::env::temp_dir().join(format!("lwp-timing-modes-{}", std::process::id()));
    std::fs::create_dir_all(&tmp).unwrap();
    let log_dir = tmp.join("logs");
    let cache_dir = tmp.join("cache");

    // Прогрев: первый запуск читает ядра с диска, дальше они в кэше системы.
    let list = groups();
    let (cold, _) = run(&list[3], Mode::Fast, &log_dir, &cache_dir);
    println!();
    println!(
        "Замер времени по режимам ({runs} запуск(а) на сценарий, {} ядер процессора, рельеф: {})",
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1),
        if Path::new("data/terrain").is_dir() { "есть" } else { "нет" }
    );
    println!("Первый (холодный) запуск: {:.2} с — чтение ядер с диска", secs(cold));
    println!("(макс — первый запуск сценария, без кэша подъёма; мин — с кэшем)");
    println!();
    println!("{:<52} {:<14} {:>7} {:>7} {:>7} {:>9}  результат", "сценарий", "режим", "мин, с", "медиана", "макс, с", "× быстр.");
    println!("{}", "-".repeat(140));

    let mut problems: Vec<String> = Vec::new();
    for g in &list {
        let _ = std::fs::remove_dir_all(&cache_dir);
        let mut fast_median = None;
        for &mode in g.modes {
            let _ = std::fs::remove_dir_all(&cache_dir);
            let mut times = Vec::with_capacity(runs);
            let mut last = String::new();
            for _ in 0..runs {
                let (t, out) = run(g, mode, &log_dir, &cache_dir);
                times.push(t);
                last = out;
            }
            times.sort();
            let (min, med, max) = (secs(times[0]), secs(times[runs / 2]), secs(times[runs - 1]));
            if mode == Mode::Fast {
                fast_median = Some(med);
            }
            let ratio = fast_median.filter(|f| *f > 0.0).map(|f| med / f);
            println!(
                "{:<52} {:<14} {:>7.2} {:>7.2} {:>7.2} {:>9}  {}",
                g.name,
                mode.name(),
                min,
                med,
                max,
                ratio.map_or("-".to_string(), |r| format!("{r:.1}")),
                result_note(mode, &last)
            );
            if let Some(p) = check_output(mode, &last) {
                problems.push(format!("«{}» / {}: {p}", g.name, mode.name()));
            }
            if med > limit_s {
                problems.push(format!("«{}» / {}: {med:.1} с — дольше предела {limit_s} с", g.name, mode.name()));
            }
            if mode != Mode::Fast {
                if let Some(r) = ratio.filter(|r| *r > ratio_limit) {
                    problems.push(format!(
                        "«{}» / {}: в {r:.1} раза медленнее быстрого — больше предела {ratio_limit}",
                        g.name,
                        mode.name()
                    ));
                }
            }
        }
    }
    println!("{}", "-".repeat(140));
    println!("Пределы: один режим — {limit_s} с, отношение к быстрому — {ratio_limit}.");
    println!();
    let _ = std::fs::remove_dir_all(&tmp);
    assert!(problems.is_empty(), "замечания:\n{}", problems.join("\n"));
}
