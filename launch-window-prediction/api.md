# API программы Lunch Window Prediction (`lwp-api/1`)

Программа считает перелёт Земля — Луна: даты старта, прилёта и посадки, траекторию, место старта
и место посадки. Другая программа получает эти данные через **командную строку с JSON**: она
запускает `launch-window-prediction api <команда>`, передаёт запрос в стандартный ввод и читает
ответ из стандартного вывода.

- Ответ — **ровно один JSON-объект** в стандартном выводе. Журнал работы программы идёт в поток
  ошибок (stderr), его можно игнорировать.
- Расчёт длится от ~1 с до ~1 мин, поэтому он идёт **задачей в фоне**: `submit` сразу
  возвращает `job_id`, дальше вы опрашиваете `status` и забираете `result`.
- Работает локально, без ключей. Запускайте из папки проекта: пути к данным (`data/…`)
  относительные.

## Команды

```text
launch-window-prediction api [--jobs-dir DIR] [--pretty] <команда> [аргументы]
```

| Команда | Что делает |
|---|---|
| `help` | описание команд и полей запроса (JSON) |
| `catalog` | справочники: космодромы, ракеты, аппараты со всеми параметрами |
| `submit [--input FILE] [--wait]` | запустить расчёт; запрос — JSON в stdin (или из файла). Сразу возвращает задачу. С `--wait` ждёт конца и возвращает задачу и сводку |
| `status <job_id>` | состояние задачи |
| `result <job_id> [--part P]` | результат готовой задачи; `P`: `summary` (по умолчанию), `trajectory`, `log`, `output`, `all` |
| `list` | все задачи, новые первыми |
| `delete <job_id>` | удалить завершённую задачу и её файлы |

Общие параметры: `--jobs-dir DIR` — папка задач (по умолчанию `jobs`), `--pretty` — ответ с
отступами (удобно смотреть глазами).

## Формат ответа

Успех:

```json
{ "ok": true, "api": "lwp-api/1", ... поля команды ... }
```

Ошибка:

```json
{ "ok": false, "api": "lwp-api/1", "error": { "code": "bad_request", "message": "unknown spaceport 'mars' (known: ksc, ...)" } }
```

| Код выхода | `error.code` | Когда |
|---|---|---|
| 0 | — | успех |
| 1 | `internal`, `job_failed` | внутренняя ошибка; `result` у упавшей задачи |
| 2 | `bad_request` | неверный запрос или команда (сообщение объясняет, что не так) |
| 3 | `not_found` | нет такой задачи |
| 4 | `not_ready` | задача ещё считается |

## Запрос на расчёт (`submit`)

Все поля необязательны, но нужна хотя бы дата или точка посадки. Неизвестные поля — ошибка
(защита от опечаток).

```json
{
  "arrival_date": "2027-02-14",
  "trajectory": "hohmann",
  "mode": "economy",
  "spaceport": "baikonur",
  "rocket": "soyuz_21b_fregat",
  "landing_site": { "lat_deg": -89.0, "lon_deg": 45.0 }
}
```

| Поле | Тип | Смысл |
|---|---|---|
| `launch_date` | `"ГГГГ-ММ-ДД"` | сутки старта (UTC): найти лучший старт в эти сутки |
| `arrival_date` | `"ГГГГ-ММ-ДД"` | сутки прилёта (UTC); с `landing_site` — **сутки посадки** в точку. Только одно из двух полей дат |
| `trajectory` | строка | тип траектории: `"direct"` (по умолчанию; прямой перелёт, задача Ламберта) или `"hohmann"` (перелёт Гомана: разгон по касательной, апогей у Луны, 180°) |
| `mode` | строка | режим полёта: `"fast"` (по умолчанию; лучший вариант в заданные сутки, как раньше) или `"economy"` (наименьший суммарный Δv миссии — разгон + торможение + спуск; поиск в пределах ±N суток от даты, секция `[economy]` в `transfer.toml`). Экономный режим считается дольше (по умолчанию 7 суток вместо одних) |
| `launch_window` | логическое | `true` — рассчитать окно старта вместо быстрого/экономного режима: все моменты старта, при которых посадка в сутки `arrival_date` обходится не дороже наименьшего Δv миссии плюс запас. Только с `arrival_date`; с `launch_date` и `mode` — ошибка `bad_request`. По умолчанию `false` |
| `launch_window_margin_m_s` | число | запас окна старта над наименьшим Δv миссии, м/с (по умолчанию `[launch_window] margin_m_s` = 50). Только вместе с `"launch_window": true` |
| `spaceport` | строка | id космодрома (см. `catalog`); без него — все подходящие |
| `rocket` | строка | id ракеты |
| `rocket_overrides` | объект | переопределения параметров ракеты, `{ "payload_tli_kg": 9000 }` |
| `spacecraft` | строка | id лунного аппарата |
| `spacecraft_overrides` | объект | переопределения параметров аппарата, `{ "payload_mass_kg": 0 }` |
| `landing_site` | `{lat_deg, lon_deg}` | точка посадки на Луне, градусы (широта + север, долгота + восток) |
| `height_m` | число | высота поверхности без карт рельефа |
| `visibility` | `{start?, hours?, step_min?}` | ряд видимости Солнца и Земли из точки посадки; `start` по умолчанию — момент посадки |
| `kernels_dir` | строка | папка ядер эфемерид (по умолчанию `data/kernels`) |
| `download_kernels` | bool | скачать ядра из интернета |
| `terrain_dir` | строка | папка карт рельефа (по умолчанию `data/terrain`, если есть; `""` — без рельефа) |
| `transfer_config`, `landing_config`, `visibility_config` | строка | свои файлы параметров |
| `no_cache` | bool | не использовать кэш моделирования подъёма |

Запрос проверяется **до** запуска: формат дат, сочетания полей, существование id космодрома,
ракеты и аппарата, имена и значения переопределений, летает ли ракета с этого космодрома, есть
ли папки данных. Ошибка возвращается сразу, задача не создаётся.

## Задача и её состояние

Ответ `submit` и `status`:

```json
{
  "ok": true, "api": "lwp-api/1",
  "job": {
    "job_id": "20261115-140841-feb5",
    "status": "done",
    "created_utc": "2026-11-15T14:08:40.975 UTC",
    "started_utc": "2026-11-15T14:08:40.977 UTC",
    "finished_utc": "2026-11-15T14:08:42.285 UTC",
    "elapsed_s": 1.307,
    "outcome": "done",
    "error": null,
    "pid": 421,
    "flight_log": "C:\\...\\jobs\\20261115-140841-feb5\\flight-2026-11-15T140841Z.json",
    "trajectory_log": "C:\\...\\jobs\\20261115-140841-feb5\\trajectory-2026-11-15T140841Z.json",
    "job_dir": "C:\\...\\jobs\\20261115-140841-feb5"
  }
}
```

| `status` | Значение |
|---|---|
| `queued` | задача создана, фоновый процесс запускается |
| `running` | считается |
| `done` | готово. **Смотрите `outcome`**: `"done"` — расчёт прошёл до конца; `"stopped: …"` — остановился по делу (например, `"stopped: landing site unsuitable"`); причина и подробности — в сводке |
| `failed` | ошибка, текст в `error`; `result --part output` и `--part log` помогут разобраться |

Опрашивайте `status` раз в 0,5–1 с. Если процесс задачи убить снаружи, задача останется
`running` — такую можно считать зависшей, если она не завершилась за несколько минут.

«Перелёт не найден» — **не ошибка**: задача `done`, а в сводке `summary.status = "no flight found"`
и `reason`.

## Результат (`result`)

```json
{ "ok": true, "api": "lwp-api/1", "job": { ... }, "part": "summary", "data": { ... } }
```

| `--part` | `data` | Размер |
|---|---|---|
| `summary` | разделы журнала полёта: `run`, `summary`, `launch_site`, `landing_site`, `files` | ~10 КБ |
| `trajectory` | траектория для 3D, формат `lwp-trajectory/1` (см. `trajectory-format.md`) | ~1 МБ |
| `log` | полный журнал расчёта, формат `lwp-flight/2` (разделы + все шаги в `log`) | ~100–200 КБ |
| `output` | текстовый вывод программы (строка) | ~5 КБ |
| `all` | `{ "summary", "trajectory", "log" }` | ~1,2 МБ |

Те же данные лежат файлами в папке задачи (пути — в `job.flight_log` и `job.trajectory_log`);
большие части можно читать прямо оттуда.

### Что внутри `summary`

- `summary` — итог коротко: `status`, `trajectory` (`"direct_transfer"` или `"hohmann"`), `mode` (`"fast"`, `"economy"` или
  `"launch_window"`), `economy` (только в экономном режиме: `requested` и `chosen` — заданные и выбранные сутки,
  `chosen_shift_days`, `saving_vs_requested_m_s`, `days` — Δv миссии по каждым просмотренным суткам),
  `launch_window` (только при `"launch_window": true`: `opens_utc` и `closes_utc` — первый и последний
  подходящий старт, `span_h`, `open_total_min`, `margin_m_s`, `min_mission_delta_v_m_s` и
  `limit_mission_delta_v_m_s`, `optimum` — самый экономный старт, `segments` — непрерывные интервалы с лучшим
  стартом в каждом, `points_checked`/`points_feasible`/`points_in_window`), космодром, ракета, аппарат, `launch_utc`,
  `lunar_orbit_insertion_utc`, `touchdown_utc`, `timeline` (все этапы с UTC и временем от
  старта), `durations`, `delta_v` (включая `mission_total_m_s` — разгон + торможение + спуск, и `tli_burn`: начало и конец разгона, потери, топливо),
  `spaceports_compared`, `variants_checked_by_refinement`; для перелёта Гомана — `hohmann`
  (номинальное время полёта, расхождение времени и допуск, полуось, перигей и апогей, промах
  апогея относительно расстояния до Луны, скорости в перигее и апогее).
- `launch_site` — космодром: координаты, положение WGS84, ограничения трасс; старт: время UTC
  и местное, Солнце над космодромом, азимут, подъём ракеты, опорная орбита, точка отрыва.
- `landing_site` — точка посадки: сторона Луны, проверка рельефа (уклоны, вердикт), высота и
  горизонт; посадка: выход на орбиту, сход, торможение, касание, промах; Солнце и Земля в момент
  касания; видимость после посадки.
- `run` — команда, время работы, итог; `files` — имена файлов журналов.

Единицы — в именах полей (`_km`, `_m_s`, `_deg`, `_utc`, …), время — UTC.

## Пример клиента (Python)

```python
import json, subprocess, time

EXE = r"target\release\launch-window-prediction.exe"   # на Linux/macOS: target/release/launch-window-prediction
CWD = r"C:\path\to\launch-window-prediction"            # папка проекта

def api(*args, request=None):
    p = subprocess.run([EXE, "api", *args], cwd=CWD, input=json.dumps(request) if request else "",
                       capture_output=True, text=True, encoding="utf-8")
    reply = json.loads(p.stdout)
    if not reply["ok"]:
        raise RuntimeError(f'{reply["error"]["code"]}: {reply["error"]["message"]}')
    return reply

# справочники
cat = api("catalog")
print([s["id"] for s in cat["spaceports"]])

# расчёт
job = api("submit", request={
    "arrival_date": "2027-02-14",
    "spaceport": "baikonur",
    "rocket": "soyuz_21b_fregat",
    "landing_site": {"lat_deg": -89.0, "lon_deg": 45.0},
})["job"]

while job["status"] in ("queued", "running"):
    time.sleep(0.5)
    job = api("status", job["job_id"])["job"]

if job["status"] == "failed":
    raise RuntimeError(job["error"])

s = api("result", job["job_id"])["data"]          # сводка
print(job["outcome"], s["summary"]["status"])
print("старт", s["summary"].get("launch_utc"), "касание", s["summary"].get("touchdown_utc"))

traj = api("result", job["job_id"], "--part", "trajectory")["data"]   # траектория для 3D
# или прочитать файл напрямую: json.load(open(job["trajectory_log"], encoding="utf-8"))
```

Синхронный вариант одним вызовом: `api("submit", "--wait", request=...)` — вернёт `job` и
`summary` после окончания расчёта.

## Пример из командной строки (Windows PowerShell)

```powershell
'{"arrival_date":"2027-01-01","spaceport":"ksc","rocket":"sls_block1","landing_site":{"lat_deg":5,"lon_deg":-10}}' |
  .\target\release\launch-window-prediction.exe api submit --wait --pretty
```

## Версия интерфейса

Поле `api` в каждом ответе: `"lwp-api/1"`. При несовместимых изменениях полей запроса или ответа
версия увеличивается. Новые поля в ответах могут добавляться без смены версии — не считайте
неизвестные поля ошибкой.
