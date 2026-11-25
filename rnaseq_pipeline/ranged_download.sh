#!/usr/bin/env bash
# Ranged byte download (a normal completed GET, no SIGPIPE). ~430 MB/mate caps the big runs
# to a few million reads; the small (4M-read) runs are fully covered since their mate < 430 MB.
BYTES=430000000
for r in SRR11742232 SRR11742231 SRR5179119 SRR1174704 SRR11729897 SRR16849223; do
  echo "=== $r ==="
  ftp=$(curl -sS -m 30 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${r}&result=read_run&fields=fastq_ftp&format=tsv" | tail -n1)
  IFS=';' read -ra urls <<< "$ftp"
  for u in "${urls[@]}"; do
    [ -z "$u" ] && continue
    fn=$(basename "$u"); out="reads/$fn"
    if [ -s "$out" ] && [ "$(wc -c <"$out")" -gt 1000 ]; then echo "exists $out"; continue; fi
    curl -sL --retry 4 --retry-delay 3 -r 0-$BYTES "https://$u" -o "$out" 2>/dev/null
    echo "done $fn: size=$(du -h "$out" | cut -f1)"
  done
done
echo "ALL_DOWNLOADS_COMPLETE"
