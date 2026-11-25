#!/usr/bin/env bash
# Download the paired-end RNA-seq runs from ENA (needs normal internet access to
# ftp.sra.ebi.ac.uk). Run on any machine that can reach ENA.
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p reads
base="https://ftp.sra.ebi.ac.uk/vol1/fastq"
dl(){ # $1 = run accession
  r=$1; p6=${r:0:6}
  n=$(printf '%03d' $((10#${r: -2})) ) # last 2 digits -> 0xx subdir (ENA convention)
  for m in 1 2; do
    url="$base/$p6/0$n/$r/${r}_${m}.fastq.gz"
    echo "Downloading $url"
    curl -fSL --retry 3 -o reads/${r}_${m}.fastq.gz "$url" || echo "WARNING: $url failed (check ENA path)"
  done
}
for r in SRR11742232 SRR11742231 SRR5179119 SRR1174704 SRR11729897 SRR16849223; do dl $r; done
echo "If a path 404s, look up the exact fastq_ftp URL at https://www.ebi.ac.uk/ena/browser/view/<RUN>"
