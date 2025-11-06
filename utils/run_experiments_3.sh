#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

PROG="task_3"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
LAB_SRC="$DIR/../lab1/${PROG}.c"
BIN="$DIR/tmp_3/${PROG}"
RESULTS="$DIR/results_3/results_${PROG}_$(date +%Y%m%d_%H%M%S).csv"
PLOT_SCRIPT="$DIR/plot_results.py"

# Параметры (процессы должны быть полными квадратами)
PROCS=(1 4)                 # можно менять (1,4,9,16,...)
SIZES=(500 1000 1500)         # размеры матриц n x n
REPEATS=3

mkdir -p "$DIR/tmp_3" "$DIR/results_3" "$DIR/graphics_3"

[ -f "$LAB_SRC" ] || { echo "Source $LAB_SRC not found" >&2; exit 1; }

# компиляция с -lm (sqrt используется)
mpicc "$LAB_SRC" -o "$BIN" -lm || { echo "mpicc failed" >&2; exit 2; }

echo "proc,samples,run,time_sec" > "$RESULTS"

for proc in "${PROCS[@]}"; do
  # quick sanity: proc must be a perfect square
  root=$(awk "BEGIN{r=sqrt($proc); printf (r==int(r)?int(r):-1)}")
  if [ "$root" -lt 0 ]; then
    echo "Skipping proc=$proc (not a perfect square)" >&2
    continue
  fi

  for n in "${SIZES[@]}"; do
    # ensure n divisible by sqrt(proc)
    if (( n % root != 0 )); then
      echo "Skipping size=$n for proc=$proc (n must be divisible by sqrt(proc)=$root)" >&2
      continue
    fi

    for ((r=1; r<=REPEATS; r++)); do
      echo "Running: proc=${proc} n=${n} run=${r}"
      output=$(mpiexec -x PMIX_MCA_gds=hash -n "$proc" "$BIN" "$n" 2>&1 || true)

      # Try to parse the Russian line:
      # "Максимальное время выполнения: %f секунд"
      time_sec=$(printf '%s\n' "$output" | grep -i "Максимальное время выполнения" | sed -E 's/.*: *([0-9]+(\.[0-9]+)?).*/\1/' | head -n1 || true)

      # fallback: try English/other patterns
      if [[ -z "$time_sec" ]]; then
        time_sec=$(printf '%s\n' "$output" | grep -i -E "Maximum time|Total time|max time" | sed -E 's/.*([0-9]+(\.[0-9]+)?).*/\1/' | head -n1 || true
)
      fi

      # If still missing, measure with /usr/bin/time
      if [[ -z "$time_sec" ]]; then
        tmpf=$(mktemp)
        /usr/bin/time -f "%e" -o "$tmpf" mpiexec -n "$proc" "$BIN" "$n" > /dev/null 2>&1 || true
        time_sec=$(cat "$tmpf" 2>/dev/null || echo "")
        rm -f "$tmpf"
      fi

      [[ -z "$time_sec" ]] && time_sec="-1"

      echo "${proc},${n},${r},${time_sec}" >> "$RESULTS"
    done
  done
done

# Построение графиков (скрипт plot_results.py должен принимать формат proc,samples,run,time_sec)
python3 "$PLOT_SCRIPT" --input "$RESULTS" --outdir "$DIR/graphics_3" || echo "plotting failed"
echo "plots: $DIR/graphics_3  results: $RESULTS"
