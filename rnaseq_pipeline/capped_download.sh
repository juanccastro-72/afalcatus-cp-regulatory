#!/usr/bin/env bash
# Use EXACT fastq_ftp URLs from the ENA API; stream + cap at CAP reads/mate.
CAP=6000000
for r in SRR11742232 SRR11742231 SRR5179119 SRR1174704 SRR11729897 SRR16849223; do
  echo "=== $r ==="
  ftp=$(curl -sS -m 30 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${r}&result=read_run&fields=fastq_ftp&format=tsv" | tail -n1)
  IFS=';' read -ra urls <<< "$ftp"
  for u in "${urls[@]}"; do
    [ -z "$u" ] && continue
    fn=$(basename "$u"); out="reads/$fn"
    [ -s "$out" ] && { echo "exists $out"; continue; }
    (curl -sL --retry 3 "https://$u" 2>/dev/null | zcat 2>/dev/null | head -n $((CAP*4)) | gzip -1 > "$out") || true
    echo "done $fn: $(zcat "$out" 2>/dev/null | wc -l | awk '{print $1/4" reads"}')  size=$(du -h "$out"|cut -f1)"
  done
done
echo "ALL_DOWNLOADS_COMPLETE"
