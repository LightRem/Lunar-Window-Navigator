#!/usr/bin/env bash
# Скачивает карты высот LOLA (LDEM_*.IMG) в data/terrain/.
# Источник: NASA PDS Geosciences Node, набор LRO-L-LOLA-4-GDR-V1.0; запасной сервер - imbrium.mit.edu.
#
# Запуск из любой папки:
#   scripts/download_terrain.sh                 # все карты (~1.2 ГБ)
#   scripts/download_terrain.sh global south    # только выбранные
# Карты: global (вся Луна, 474 м), south (от 80° ю. ш., 80 м), north (от 80° с. ш., 80 м),
#        south10m (от 87.5° ю. ш., 10 м), all.
# Папку можно сменить: TERRAIN_DIR=/путь scripts/download_terrain.sh
# Уже скачанные карты с верной контрольной суммой пропускаются; прерванная загрузка продолжается.
set -euo pipefail

PDS=https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr
MIT=https://imbrium.mit.edu/DATA/LOLA_GDR

# ключ | файл | размер, байт | sha256 | путь на PDS (без расширения) | путь на MIT (без расширения)
MAPS="
global|LDEM_64|530841600|1c4958699d4cffd7e777d51421b6044bba300aa0bef286a544309efb082cdbd6|cylindrical/img/ldem_64|CYLINDRICAL/IMG/LDEM_64
south|LDEM_80S_80M|115520000|e6591eeea35402ca5e124ca00f5bf0774381a2d4fa919144b455d82b8ff65fdc|polar/img/ldem_80s_80m|POLAR/IMG/LDEM_80S_80M
north|LDEM_80N_80M|115520000|dfaa75f2d1337258eaf5e13c0195fcad0d70794911b2daf689026ed12dce2233|polar/img/ldem_80n_80m|POLAR/IMG/LDEM_80N_80M
south10m|LDEM_875S_10M|460136448|3c774fd3d8fcab8a9f2e40160821ec9ca83d44ea616f5b7cd13fa78a89f39c53|polar/img/ldem_875s_10m|POLAR/IMG/LDEM_875S_10M
"

DIR="${TERRAIN_DIR:-$(cd "$(dirname "$0")/.." && pwd)/data/terrain}"
mkdir -p "$DIR"

command -v curl >/dev/null || { echo "нужен curl" >&2; exit 1; }
if command -v sha256sum >/dev/null; then sha() { sha256sum "$1" | cut -d' ' -f1; }
elif command -v shasum >/dev/null; then sha() { shasum -a 256 "$1" | cut -d' ' -f1; }
else echo "нужен sha256sum или shasum" >&2; exit 1; fi
size() { wc -c < "$1" | tr -d ' '; }

want=("$@"); [ ${#want[@]} -eq 0 ] && want=(all)
for w in "${want[@]}"; do
  case "$w" in global|south|north|south10m|all) ;; *) echo "неизвестная карта '$w' (global, south, north, south10m, all)" >&2; exit 2;; esac
done
selected() { for w in "${want[@]}"; do [ "$w" = all ] || [ "$w" = "$1" ] && return 0; done; return 1; }

# скачать url в файл, продолжая прерванную загрузку; если докачка не вышла - заново с нуля
fetch() {
  curl -fL --retry 3 --retry-delay 5 -C - -o "$2" "$1" && return 0
  rm -f "$2"
  curl -fL --retry 3 --retry-delay 5 -o "$2" "$1"
}

failed=0
while IFS='|' read -r key name bytes hash pds mit; do
  [ -z "$key" ] && continue
  selected "$key" || continue
  img="$DIR/$name.IMG"; lbl="$DIR/$name.LBL"; part="$img.part"

  if [ ! -f "$lbl" ]; then
    echo "== $name.LBL: скачиваю описание"
    fetch "$PDS/$pds.lbl" "$lbl" || fetch "$MIT/$mit.LBL" "$lbl" || { echo "!! не удалось скачать $name.LBL" >&2; failed=1; }
  fi

  if [ -f "$img" ] && [ "$(size "$img")" = "$bytes" ]; then
    echo "== $name.IMG: уже есть, проверяю контрольную сумму"
    if [ "$(sha "$img")" = "$hash" ]; then echo "   в порядке"; continue; fi
    echo "   сумма не совпала, скачиваю заново"
    rm -f "$img"
  fi

  echo "== $name.IMG: скачиваю $((bytes / 1048576)) МБ"
  # .part больше нужного размера - испорчен, начинаем заново
  if [ -f "$part" ] && [ "$(size "$part")" -gt "$bytes" ]; then rm -f "$part"; fi
  if ! fetch "$PDS/$pds.img" "$part"; then
    echo "   PDS недоступен, пробую imbrium.mit.edu"
    fetch "$MIT/$mit.IMG" "$part" || { echo "!! не удалось скачать $name.IMG" >&2; failed=1; continue; }
  fi
  if [ "$(size "$part")" != "$bytes" ] || [ "$(sha "$part")" != "$hash" ]; then
    echo "!! $name.IMG: размер или контрольная сумма не совпали; удалите $part и запустите снова" >&2
    failed=1; continue
  fi
  mv -f "$part" "$img"
  echo "   готово: $img"
done <<< "$MAPS"

[ $failed -eq 0 ] && echo "Карты в $DIR готовы." || { echo "Часть карт не скачалась, см. сообщения выше." >&2; exit 1; }
