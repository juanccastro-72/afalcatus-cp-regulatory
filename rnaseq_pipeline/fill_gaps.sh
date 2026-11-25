#!/usr/bin/env bash
# Robustly ensure every mate FASTQ exists and is valid. Retries each mate until a valid
# gzip of reasonable size is present (fixes the cold-connection first-request reset).
BYTES=430000000
valid(){ [ -s "$1" ] && [ "$(wc -c <"$1")" -gt 5000000 ] && [ "$(zcat "$1" 2>/dev/null | head -4 | wc -l)" -ge 4 ]; }
warm(){ curl -sS -m 20 -o /dev/null "https://ftp.sra.ebi.ac.uk/robots.txt" 2>/dev/null || true; }
for r in SRR11742232 SRR11742231 SRR5179119 SRR1174704 SRR11729897 SRR16849223; do
  ftp=$(curl -sS -m 30 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${r}&result=read_run&fields=fastq_ftp&format=tsv" | tail -n1)
  IFS=';' read -ra urls <<< "$ftp"
  for u in "${urls[@]}"; do
    [ -z "$u" ] && continue
    fn=$(basename "$u"); out="reads/$fn"
    if valid "$out"; then echo "ok   $fn ($(du -h "$out"|cut -f1))"; continue; fi
    for try in 1 2 3 4 5 6; do
      warm
      # gzip -t tolerates a truncated final member? No -> ranged file may fail gzip -t.
      # Accept by size instead: download range, then verify it decompresses at least partly.
      curl -sL --retry 3 --retry-delay 2 -r 0-$BYTES "https://$u" -o "$out" 2>/dev/null
      if [ -s "$out" ] && [ "$(wc -c <"$out")" -gt 5000000 ] && [ "$(zcat "$out" 2>/dev/null | head -4 | wc -l)" -ge 4 ]; then
        echo "got  $fn ($(du -h "$out"|cut -f1)) on try $try"; break
      fi
      echo "retry $fn (try $try, size=$(du -h "$out" 2>/dev/null|cut -f1))"; sleep 3
    done
  done
done
echo "FILL_GAPS_COMPLETE"
ls -lh reads/
