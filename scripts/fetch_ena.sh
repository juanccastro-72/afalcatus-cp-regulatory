#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../data/genbank"
tail -n +2 ../genome_manifest_input.tsv | cut -f5 | while read acc; do
  [ "$acc" = "QB176266.1" ] && continue
  [ -s "$acc.embl" ] && continue
  curl -sS -m 120 -o "$acc.embl" "https://www.ebi.ac.uk/ena/browser/api/embl/$acc?download=true&expanded=true"
  echo "$acc $(head -1 $acc.embl | cut -c1-80)"
done
