#!/usr/bin/env bash
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
mode=${1:-status}
case "$mode" in status|start|collect) ;; *) exit 2;; esac
command -v kaggle >/dev/null || { echo 'Run in the Kaggle-authenticated Codespace with the kaggle CLI installed.' >&2; exit 1; }
ps -eo pid,etimes,pcpu,rss,comm | awk 'NR==1 || $5 ~ /cadical|kissat|plingeling|a16_sls|ext-oracle/'
if [[ $mode != status ]]; then
  git fetch -q origin main
  [[ $(git rev-parse HEAD) == $(git rev-parse FETCH_HEAD) ]] || { echo "Use a clean worktree at origin/main."; exit 1; }
fi
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
stamp=$(date -u +%Y%m%dT%H%M%SZ)
archive="kaggle_rollover/extension_${stamp}_$$"
for dir in kaggle_a16 kaggle_a16_w{2,3,4,5}; do
  id=$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["id"])' "$dir/kernel-metadata.json")
  status=$(kaggle kernels status "$id")
  echo "$id $status"
  [[ $mode != status ]] || continue
  case "${status,,}" in *running*|*queued*|*pending*) continue;; *complete*|*error*|*cancel*) ;; *) echo "Unrecognized status; skipping $id"; continue;; esac
  dest="$archive/${id##*/}"
  mkdir -p "$dest"
  echo "$status" > "$dest/status.txt"
  # Never overwrite a version whose output could not be retrieved.
  kaggle kernels output "$id" -p "$dest" || { echo "Download failed; retaining $id"; continue; }
  # Keep the reproducible 143 MB graph locally, outside the Git-tracked output folder.\n  while IFS= read -r -d '' graph; do
    relative=${graph#"$dest"/}
    graph_copy="$archive/graphs/${id##*/}/$relative"
    mkdir -p "$(dirname "$graph_copy")"
    mv "$graph" "$graph_copy"
  done < <(find "$dest" -type f -name A16_graph.bin -print0)\n  [[ -n $(find "$dest" -type f ! -name status.txt -print -quit) ]] || { echo "No output files; retaining $id"; continue; }
  # Commit only this archive. Push must succeed before reusing its kernel.
  git add -f -- "$dest"
  git commit --only -m "Archive Kaggle outputs for $id ($stamp)" -- "$dest"
  git push origin HEAD:main
  [[ $mode == start && $dir != kaggle_a16 ]] || continue
  # A completed extension shard is collected, never automatically restarted.
  if find "$dest" -name extension_status.json -print -quit | grep -q .; then continue; fi
  pkg="$tmp/${id##*/}"
  mkdir "$pkg"
  python3 extension/package.py "$dir" "$pkg" "$(( ${dir##*w} - 2 ))"
  kaggle kernels push -p "$pkg"
done
if [[ $mode == start ]]; then
  echo "Waiting for jobs to finish before collecting outputs. Keep this shell or tmux session alive."
  while sleep 120; do
    busy=0
    for meta in kaggle_a16{,_w2,_w3,_w4,_w5}/kernel-metadata.json; do
      id=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["id"])' "$meta")
      status=$(kaggle kernels status "$id") || { busy=1; continue; }
      case "${status,,}" in *complete*|*error*|*cancel*) ;; *) busy=1;; esac
    done
    (( busy )) || break
  done
  bash extension/rollover.sh collect
fi
