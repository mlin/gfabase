#!/bin/bash

set -o pipefail

REPO="$(dirname "$0")/.."
cd "$REPO"
export BASH_TAP_ROOT=test/bash-tap
source test/bash-tap/bash-tap-bootstrap
export LC_ALL=C

plan tests 8

# aria2c -x 10 -j 10 -s 10 https://glennhickey.s3.amazonaws.com/share/GRCh38-f1-90-mc-mar13.gfa.gz
# pv GRCh38-f1-90-mc-mar13.gfa.gz | bgzip -dc | gfabase load  -o GRCh38-f1-90-mc-mar13.gfab --verbose --memory-gbytes 16
# gfabase sub GRCh38-f1-90-mc-mar13.gfab --view --connected --range GRCh38.chr21:1-1000000000 GRCh38.chr22:1-1000000000 | zstd -T0 -19 > GRCh38-f1-90-mc-mar13.chr21_chr22.zst
aria2c -c -d /tmp -s 16 -x 16 --retry-wait 2 https://github.com/mlin/gfabase/releases/download/v0.6.0/GRCh38-f1-90-mc-mar13.chr21_chr22.gfa.zst
is "$?" "0" "download gfa"

gfabase="cargo run --release -- --verbose"

export TMPDIR=$(mktemp -d --tmpdir gfabase_vg_walks_test_XXXXXX)

zstd -dc /tmp/GRCh38-f1-90-mc-mar13.chr21_chr22.gfa.zst | grep ^W | sort > "${TMPDIR}/original_walks" & pid=$!
zstd -dc /tmp/GRCh38-f1-90-mc-mar13.chr21_chr22.gfa.zst | time $gfabase load -o "${TMPDIR}/GRCh38-f1-90-mc-mar13.chr21_chr22.gfab" --compress 1
is "$?" "0" "gfabase load"
ls -lh "${TMPDIR}/GRCh38-f1-90-mc-mar13.chr21_chr22.gfab"

$gfabase view "${TMPDIR}/GRCh38-f1-90-mc-mar13.chr21_chr22.gfab" | grep ^W | sort | sha256sum > "${TMPDIR}/view_walks"
is "$?" "0" "gfabase view"

wait $pid
grep 'CHM13\|HG02148' "${TMPDIR}/original_walks" | sha256sum > "${TMPDIR}/original_sub_walks" & pid=$!
is "$(cat "${TMPDIR}/view_walks")" "$(cat "${TMPDIR}/original_walks" | sha256sum)" "roundtrip walks"

$gfabase sub "${TMPDIR}/GRCh38-f1-90-mc-mar13.chr21_chr22.gfab" --walk-samples CHM13,HG02148 \
    --view --connected --range GRCh38.chr21:1-1000000000 GRCh38.chr22:1-1000000000 \
    | grep ^W | sort | sha256sum > "${TMPDIR}/sub_walks"
is "$?" "0" "gfabase sub"

wait $pid
is "$(cat "${TMPDIR}/sub_walks")" "$(cat "${TMPDIR}/original_sub_walks")" "sub walks"

gsql() {
    genomicsqlite "${TMPDIR}/GRCh38-f1-90-mc-mar13.chr21_chr22.gfab" "$1" | tail -n 1
}
export -f gsql

is "$(gsql 'select count(1) from gfa1_walk_connectivity where component_id = 56712591')" "1052" "walk count 1"
is "$(gsql 'select count(1) from gfa1_walk_connectivity where component_id = 61550138')" "503" "walk count 2"

rm -rf "$TMPDIR"
