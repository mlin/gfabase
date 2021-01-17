#!/bin/bash

set -o pipefail

REPO="$(dirname "$0")/.."
cd "$REPO"
export BASH_TAP_ROOT=test/bash-tap
source test/bash-tap/bash-tap-bootstrap
export LC_ALL=C

plan tests 13

aria2c -c -d /tmp -o shasta-HG002-Guppy-3.6.0-run4-UL.gfa -s 16 -x 16 --retry-wait 2 \
    https://s3-us-west-2.amazonaws.com/czi.paolo-public/HG002-Guppy-3.6.0-run4-UL/Assembly.gfa
# SHA256: be83e6d77a848834811406fef0a53ad29eee035e3d2446b3e55fb03b33444cf0
is "$?" "0" "download gfa"

aria2c -c -d /tmp -o shasta-HG002-Guppy-3.6.0-run4-UL.paf -s 16 -x 16 --retry-wait 2 \
    https://s3-us-west-2.amazonaws.com/czi.paolo-public/HG002-Guppy-3.6.0-run4-UL/mapToHg38/Alignments-NoCigar.tsv
# a1153fd9f4de014cf4c892b1ff382a3247cbbfc29d593245d0f879a14e725449
is "$?" "0" "download paf"

gfabase="cargo run --release -- --verbose"

export TMPDIR=$(mktemp -d --tmpdir gfabase_shasta_test_XXXXXX)

time $gfabase load /tmp/shasta-HG002-Guppy-3.6.0-run4-UL.gfa -o "${TMPDIR}/Assembly.gfab" --compress 1
is "$?" "0" "gfabase load"

time $gfabase view "${TMPDIR}/Assembly.gfab" | wc -c
is "$?" "0" "gfabase view"

ls -lh "${TMPDIR}/Assembly.gfab"

gsql() {
    genomicsqlite "${TMPDIR}/Assembly.gfab" "$1"| tail -n 1
}
export -f gsql

is "$(gsql 'select count(1) from gfa1_segment_meta')"  "12520" "gfab segment count"
is "$(gsql 'select count(1) from gfa1_link')" "18456" "gfab link count"
is "$(gsql 'select sum(sequence_length) from gfa1_segment_meta')" "5100857444" "gfab base metacount"
is "$(gsql 'select sum(twobit_length(sequence_twobit)) from gfa1_segment_sequence')" "5100857444" "gfab base count"

#cp "${TMPDIR}/Assembly.gfab" "${TMPDIR}/Assembly.gfab.bak"
time $gfabase add-mappings "${TMPDIR}/Assembly.gfab" /tmp/shasta-HG002-Guppy-3.6.0-run4-UL.paf --quality 60
is "$?" "0" "add mappings"
is "$(gsql 'select count(1) from gfa1_segment_mapping')" "9180" "mapQ60 count"
time $gfabase add-mappings "${TMPDIR}/Assembly.gfab" /tmp/shasta-HG002-Guppy-3.6.0-run4-UL.paf --quality 30 --replace
is "$?" "0" "replace mappings"
is "$(gsql 'select count(1) from gfa1_segment_mapping')" "10238" "mapQ60 count"

is $($gfabase sub "${TMPDIR}/Assembly.gfab" --range --cutpoints 1 chr12:111766933-111817532 | grep "^S" | cut -f3 | LC_ALL=C sort | sha256sum | cut -f1 -d ' ') \
   "67cfe2746f311654a87cf11050e2711d52e26d433f146e6a527258286746af05" "ALDH2 segments"

rm -rf "$TMPDIR"
