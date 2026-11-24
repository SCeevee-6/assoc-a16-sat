#!/usr/bin/env bash
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
git pull --ff-only origin main
command -v kaggle >/dev/null
command -v gcc >/dev/null
command -v g++ >/dev/null
command -v cadical >/dev/null || { echo "Install CaDiCaL first (sudo apt-get update && sudo apt-get install -y cadical)." >&2; exit 1; }

run_dir="batch3_codespaces/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$run_dir"
gcc -O3 -Wall -Wextra -Werror=int-conversion build_a16.c -o "$run_dir/build_a16"
"$run_dir/build_a16" 16 "$run_dir/A16_sym.cnf"
# One unrestricted CaDiCaL run with a positional DRAT proof output (CaDiCaL 1.7.4 syntax).
(cadical --seed=11 "$run_dir/A16_sym.cnf" "$run_dir/cadical.drat" >"$run_dir/cadical.log" 2>&1) &

# Five Kaggle kernels x four independent prefix shards = 20 extension workers.
# Together with the four CaDiCaL processes this occupies 24 CPU workers.
dirs=(kaggle_a16_w2 kaggle_a16_w3 kaggle_a16_w4 kaggle_a16_w5 kaggle_a16)
for index in ${!dirs[@]}; do
  dir=${dirs[$index]}
  id=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["id"])' "$dir/kernel-metadata.json")
  status=$(kaggle kernels status "$id")
  case "${status,,}" in *complete*|*error*|*cancel*) ;; *) echo "$id is not finished: $status" >&2; exit 1;; esac
  pkg=$(mktemp -d)
  python3 extension/package.py "$dir" "$pkg" "$index"
  kaggle kernels push -p "$pkg"
  rm -rf "$pkg"
done

echo "Launched 20 seed-extension workers and one proof-logging CaDiCaL job."
echo "Live SAT log: $run_dir/cadical.log"
echo "Proof file:    $run_dir/cadical.drat"
echo "Monitor Kaggle with: bash extension/rollover.sh status"

# Wait for Kaggle jobs to finish, then preserve their outputs before reuse.
while :; do
  busy=0
  for dir in "${dirs[@]}"; do
    id=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["id"])' "$dir/kernel-metadata.json")
    status=$(kaggle kernels status "$id") || { busy=1; continue; }
    case "${status,,}" in *complete*|*error*|*cancel*) ;; *) busy=1;; esac
  done
  (( busy == 0 )) && break
  sleep 120
done

stamp=$(date -u +%Y%m%dT%H%M%SZ)
archive="kaggle_rollover/seed_extension_$stamp"
for dir in "${dirs[@]}"; do
  id=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["id"])' "$dir/kernel-metadata.json")
  dest="$archive/${id##*/}"
  mkdir -p "$dest"
  kaggle kernels output "$id" -p "$dest"
  while IFS= read -r -d '' graph; do
    relative=${graph#"$dest"/}
    graph_copy="$archive/graphs/${id##*/}/$relative"
    mkdir -p "$(dirname "$graph_copy")"
    mv "$graph" "$graph_copy"
  done < <(find "$dest" -type f -name A16_graph.bin -print0)
  git add -f -- "$dest"
  git commit --only -m "Archive A16 seed-extension output for $id ($stamp)" -- "$dest"
  git push origin HEAD:main
done
