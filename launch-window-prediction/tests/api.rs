//! Проверка машинного интерфейса (`launch-window-prediction api ...`) так, как его вызывает
//! другая программа: запуск процесса, JSON на входе, JSON на выходе, коды выхода.
//!
//! Тесты без расчёта идут всегда. Полный цикл задачи (submit → status → result) требует ядер в
//! `data/kernels` и рельефа в `data/terrain`, поэтому помечен `#[ignore]`:
//!
//! ```text
//! cargo test --release --test api -- --include-ignored
//! ```

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

use serde_json::Value;

const BIN: &str = env!("CARGO_BIN_EXE_launch-window-prediction");

fn jobs_dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("lwp-api-test-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    d
}

/// Запускает `api <args>` с `stdin` на входе; возвращает (код выхода, JSON ответа).
fn api(args: &[&str], stdin: &str) -> (i32, Value) {
    let mut child = Command::new(BIN)
        .arg("api")
        .args(args)
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("run the program");
    child.stdin.take().unwrap().write_all(stdin.as_bytes()).unwrap();
    let out = child.wait_with_output().unwrap();
    let text = String::from_utf8(out.stdout).unwrap();
    let v: Value = serde_json::from_str(&text).unwrap_or_else(|e| panic!("stdout is not JSON ({e}): {text}"));
    (out.status.code().unwrap_or(-1), v)
}

#[test]
fn help_and_errors_are_json_with_exit_codes() {
    let (code, v) = api(&["help"], "");
    assert_eq!(code, 0);
    assert_eq!(v["ok"], true);
    assert_eq!(v["api"], "lwp-api/1");
    assert!(v["commands"].as_object().unwrap().len() >= 7);

    let (code, v) = api(&["frobnicate"], "");
    assert_eq!((code, v["ok"].clone(), v["error"]["code"].clone()), (2, Value::Bool(false), Value::from("bad_request")));

    let dir = jobs_dir("errors");
    let d = dir.to_str().unwrap();
    for (body, needle) in [
        ("", "empty request"),
        ("{\"foo\": 1}", "unknown field"),
        ("{\"arrival_date\": \"2027-13-01\"}", "YYYY-MM-DD"),
        ("{\"spaceport\": \"ksc\"}", "nothing to calculate"),
        ("{\"launch_date\": \"2027-01-01\", \"arrival_date\": \"2027-01-05\"}", "cannot be used with"),
        ("{\"launch_date\": \"2027-01-01\", \"spaceport\": \"mars\"}", "unknown spaceport"),
        ("{\"launch_date\": \"2027-01-01\", \"mode\": \"slow\"}", "unknown mode"),
        ("{\"launch_date\": \"2027-01-01\", \"launch_window\": true}", "only for a landing day"),
        ("{\"arrival_date\": \"2027-02-14\", \"launch_window\": true, \"mode\": \"economy\"}", "replaces mode"),
        ("{\"arrival_date\": \"2027-02-14\", \"launch_window_margin_m_s\": 30}", "needs"),
        ("{\"launch_date\": \"2027-01-01\", \"rocket\": \"sls_block1\", \"spaceport\": \"baikonur\"}", "does not launch from"),
        ("{\"launch_date\": \"2027-01-01\", \"kernels_dir\": \"no/such/dir\"}", "does not exist"),
    ] {
        let (code, v) = api(&["--jobs-dir", d, "submit"], body);
        assert_eq!(code, 2, "{body}: {v}");
        assert_eq!(v["error"]["code"], "bad_request");
        assert!(v["error"]["message"].as_str().unwrap().contains(needle), "{body}: {v}");
    }
    // Неверный запрос не создаёт задач.
    let (_, v) = api(&["--jobs-dir", d, "list"], "");
    assert_eq!(v["jobs"].as_array().unwrap().len(), 0);
    let (code, v) = api(&["--jobs-dir", d, "status", "20000101-000000-0000"], "");
    assert_eq!((code, v["error"]["code"].clone()), (3, Value::from("not_found")));
    let (code, _) = api(&["--jobs-dir", d, "status", "../etc"], "");
    assert_eq!(code, 2);
}

#[test]
fn catalog_lists_the_reference_data() {
    let (code, v) = api(&["catalog"], "");
    assert_eq!(code, 0, "{v}");
    assert_eq!(v["spaceports"].as_array().unwrap().len(), 11);
    let ids = |k: &str| v[k].as_array().unwrap().iter().map(|e| e["id"].as_str().unwrap().to_string()).collect::<Vec<_>>();
    assert!(ids("rockets").contains(&"soyuz_21b_fregat".to_string()));
    assert!(ids("spacecraft").contains(&"apollo_lm_descent".to_string()));
    assert!(v["rockets"][0]["resolved"]["config"]["upper_stage_thrust_n"].is_number());
}

#[test]
#[ignore = "needs kernels in data/kernels and terrain in data/terrain"]
fn job_runs_in_the_background_and_returns_all_parts() {
    let dir = jobs_dir("job");
    let d = dir.to_str().unwrap();
    let request = r#"{
        "arrival_date": "2027-01-01", "spaceport": "ksc", "rocket": "sls_block1",
        "spacecraft": "apollo_lm_descent", "spacecraft_overrides": {"payload_mass_kg": 0},
        "landing_site": {"lat_deg": 5, "lon_deg": -10}
    }"#;
    let (code, v) = api(&["--jobs-dir", d, "submit"], request);
    assert_eq!(code, 0, "{v}");
    let id = v["job"]["job_id"].as_str().unwrap().to_string();
    assert!(matches!(v["job"]["status"].as_str(), Some("queued" | "running")));

    // Опрос до готовности.
    let mut status = String::new();
    for _ in 0..240 {
        let (_, s) = api(&["--jobs-dir", d, "status", &id], "");
        status = s["job"]["status"].as_str().unwrap().to_string();
        if status == "done" || status == "failed" {
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(250));
    }
    assert_eq!(status, "done");

    let (code, v) = api(&["--jobs-dir", d, "result", &id], "");
    assert_eq!(code, 0);
    let s = &v["data"]["summary"];
    assert_eq!(s["status"], "flight and landing found");
    assert!(s["touchdown_utc"].as_str().unwrap().starts_with("2027-01-01"));
    assert!(v["data"]["launch_site"]["launch"]["possible"].as_bool().unwrap());
    assert!(v["data"]["landing_site"]["landing"]["possible"].as_bool().unwrap());

    let (_, v) = api(&["--jobs-dir", d, "result", &id, "--part", "trajectory"], "");
    assert_eq!(v["data"]["format"], "lwp-trajectory/1");
    let (_, v) = api(&["--jobs-dir", d, "result", &id, "--part", "log"], "");
    assert_eq!(v["data"]["format"], "lwp-flight/2");
    let (_, v) = api(&["--jobs-dir", d, "result", &id, "--part", "output"], "");
    assert!(v["data"].as_str().unwrap().contains("Flight log"));

    // Синхронный вариант и задача, которая останавливается на непригодной площадке.
    let (code, v) = api(
        &["--jobs-dir", d, "submit", "--wait"],
        r#"{"arrival_date": "2027-02-14", "spaceport": "baikonur", "rocket": "soyuz_21b_fregat", "landing_site": {"lat_deg": -89.5, "lon_deg": 0}}"#,
    );
    assert_eq!(code, 0, "{v}");
    assert_eq!(v["job"]["status"], "done");
    assert_eq!(v["job"]["outcome"], "stopped: landing site unsuitable");
    assert_eq!(v["summary"]["summary"]["status"], "landing site unsuitable");

    let (_, v) = api(&["--jobs-dir", d, "list"], "");
    assert_eq!(v["jobs"].as_array().unwrap().len(), 2);
    let (code, _) = api(&["--jobs-dir", d, "delete", &id], "");
    assert_eq!(code, 0);
    let _ = std::fs::remove_dir_all(&dir);
}
