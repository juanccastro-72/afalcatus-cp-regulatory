#!/usr/bin/env bash
# Robust capped downloader: hardcoded exact ENA URLs, warm-up, retry-until-valid per mate.
# 200 MB/mate (~2-3M reads) is ample for plastome coverage (plastome is a small % of RNA).
BYTES=200000000
B="https://ftp.sra.ebi.ac.uk/vol1/fastq"
declare -a M=(
  "SRR117/032/SRR11742232/SRR11742232"
  "SRR117/031/SRR11742231/SRR11742231"
  "SRR517/009/SRR5179119/SRR5179119"
  "SRR117/004/SRR1174704/SRR1174704"
  "SRR117/097/SRR11729897/SRR11729897"
  "SRR168/023/SRR16849223/SRR16849223"
)
ok(){ [ -s "$1" ] && [ "$(wc -c <"$1")" -gt 5000000 ] && [ "$(zcat "$1" 2>/dev/null | head -4 | wc -l)" -ge 4 ]; }
# warm up the connection to the SRA host once
curl -sS -m 20 -o /dev/null "$B/" 2>/dev/null || true
for path in "${M[@]}"; do
  for mate in 1 2; do
    url="$B/${path}_${mate}.fastq.gz"; out="reads/$(basename "$url")"
    if ok "$out"; then echo "ok   $(basename "$out") ($(du -h "$out"|cut -f1))"; continue; fi
    for try in 1 2 3 4 5 6; do
      curl -sL --retry 3 --retry-delay 2 --connect-timeout 20 -r 0-$BYTES "$url" -o "$out" 2>/dev/null
      if ok "$out"; then echo "got  $(basename "$out") ($(du -h "$out"|cut -f1)) try $try"; break; fi
      echo "retry $(basename "$out") try $try"; sleep 3
    done
  done
done
echo "GET_READS_COMPLETE"
ls -lh reads/
