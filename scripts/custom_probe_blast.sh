#!/usr/bin/env bash
#
# BLAST nominated LHS/RHS probe arms against a reference transcriptome and keep
# hits long enough to matter as off-target priming.
#
# Usage:
#   ./probe_blast.sh <query.fasta> <output.txt> [-d database] [-m min_length] [-t threads]
#
#
# Output table for read_blast_results() in probe_design.py.

set -euo pipefail

database="hs.refseq.rna"
min_length=20          # minimum alignment length to retain
threads=$(nproc 2>/dev/null || echo 1)
allow_empty=0

usage() {
    sed -n '3,12p' "$0" | sed 's/^# \?//'
    exit "${1:-1}"
}

if [[ $# -lt 2 ]]; then
    usage
fi

query=$1
output=$2
shift 2

while getopts ":d:m:t:eh" opt; do
    case $opt in
        d) database=$OPTARG ;;
        m) min_length=$OPTARG ;;
        t) threads=$OPTARG ;;
        e) allow_empty=1 ;;
        h) usage 0 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
        :)  echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# --- Preflight ---------
command -v blastn >/dev/null || { echo "ERROR: blastn not found on PATH." >&2; exit 1; }

[[ -s $query ]] || { echo "ERROR: query FASTA missing or empty: $query" >&2; exit 1; }

blastdbcmd -db "$database" -info >/dev/null 2>&1 || {
    echo "ERROR: BLAST database not found or unreadable: $database" >&2
    echo "       Set BLASTDB, or pass -d /path/to/db, or build it with make_hs_refseq_blastdb.sh." >&2
    exit 1
}

output_dir=$(dirname "$output")
mkdir -p "$output_dir"

n_queries=$(grep -c '^>' "$query")
echo "Querying $n_queries probe arms against $database using $threads thread(s)..."

# --- Run -----------------------------------------------------------------
# Written to a temp file and moved into place only on success, 
tmp=$(mktemp "${output}.XXXXXX")
trap 'rm -f "$tmp"' EXIT

blastn \
    -db "$database" \
    -query "$query" \
    -task blastn-short \
    -dust no \
    -num_threads "$threads" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    | awk -v min="$min_length" '$4 >= min' \
    > "$tmp"

n_hits=$(wc -l < "$tmp")

if [[ $n_hits -eq 0 && $allow_empty -eq 0 ]]; then
    echo "ERROR: zero hits passed the length filter. For a transcriptome-wide" >&2
    echo "       search this almost always means the wrong database, not a clean" >&2
    echo "       probe set. Re-run with -e if zero hits is genuinely expected." >&2
    exit 1
fi

mv "$tmp" "$output"
trap - EXIT

echo "Retained $n_hits hits of length >= ${min_length}: $output"