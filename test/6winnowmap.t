#!/bin/bash

set -o pipefail

REPO="$(dirname "$0")/.."
cd "$REPO"
export BASH_TAP_ROOT=test/bash-tap
source test/bash-tap/bash-tap-bootstrap
export LC_ALL=C

plan tests 4

aria2c -c -d /tmp -s 16 -x 16 --retry-wait 2 \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
is "$?" "0" "download GRCh38"

cargo build --release
gfabase="$(pwd)/target/release/gfabase"

export TMPDIR=$(mktemp -d --tmpdir gfabase_winnowmap_test_XXXXXX)

zstd -dc test/data/shasta-HG002-Guppy-3.6.0-run4-UL.chr21.gfa.zst | \
    $gfabase --verbose load -o "${TMPDIR}/shasta-HG002.chr21.gfab"
is "$?" "0" "gfabase load"

miniwdl run workflows/gfab_winnowmap.wdl \
    gfab="${TMPDIR}/shasta-HG002.chr21.gfab" \
    reference_fasta_gz=/tmp/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
    meryl_repetitive_kmers_gz=test/data/GRCh38_no_alt_analysis_set.repetitive.txt.gz \
    gfabase_exe="$gfabase" \
    add_mappings.quality=30 \
    --dir "$TMPDIR" --verbose
is "$?" "0" "miniwdl run"

gsql() {
    genomicsqlite "${TMPDIR}/_LAST/out/mapped_gfab/shasta-HG002.chr21.GRCh38_no_alt_analysis_set.gfab" "$1"| tail -n 1
}
export -f gsql
is "$(gsql 'select count(1) from gfa1_segment_mapping')" "173" "mappings count"

rm -rf "$TMPDIR"
