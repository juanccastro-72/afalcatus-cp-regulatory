#!/usr/bin/env bash
# Parallel capped fetch of remaining mates. 60 MB/mate (~0.8M reads) each.
BYTES=60000000
B="https://ftp.sra.ebi.ac.uk/vol1/fastq"
declare -a U=(
  "SRR117/031/SRR11742231/SRR11742231_2"
  "SRR517/009/SRR5179119/SRR5179119_1"
  "SRR517/009/SRR5179119/SRR5179119_2"
  "SRR117/004/SRR1174704/SRR1174704_1"
  "SRR117/004/SRR1174704/SRR1174704_2"
  "SRR117/097/SRR11729897/SRR11729897_1"
  "SRR117/097/SRR11729897/SRR11729897_2"
  "SRR168/023/SRR16849223/SRR16849223_1"
  "SRR168/023/SRR16849223/SRR16849223_2"
)
ok(){ [ -s "$1" ] && [ "$(wc -c <"$1")" -gt 5000000 ] && [ "$(zcat "$1" 2>/dev/null | head -4 | wc -l)" -ge 4 ]; }
curl -sS -m 20 -o /dev/null "$B/" 2>/dev/null || true   # warm up host
fetch(){
  local path="$1" url out
  url="$B/${path}.fastq.gz"; out="reads/$(basename "$url")"
  for try in 1 2 3 4 5 6; do
    ok "$out" && { echo "ok $(basename "$out")"; return; }
    curl -sL --retry 3 --retry-delay 2 --connect-timeout 20 -r 0-$BYTES "$url" -o "$out" 2>/dev/null
    ok "$out" && { echo "got $(basename "$out") ($(du -h "$out"|cut -f1)) try $try"; return; }
    sleep 3
  done
  echo "FAIL $(basename "$out")"
}
for p in "${U[@]}"; do fetch "$p" & done
wait
echo "PAR_FETCH_COMPLETE"
ls -lh reads/
