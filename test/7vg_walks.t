#!/bin/bash

set -o pipefail

REPO="$(dirname "$0")/.."
cd "$REPO"
export BASH_TAP_ROOT=test/bash-tap
source test/bash-tap/bash-tap-bootstrap
export LC_ALL=C

plan tests 4

aria2c -c -d /tmp -s 16 -x 16 --retry-wait 2 https://github.com/mlin/gfabase/releases/download/v0.5.0/GRCh38-f1-90-mc-mar13.chr21_chr22.zst
is "$?" "0" "download gfa"

gfabase="cargo run --release -- --verbose"

export TMPDIR=$(mktemp -d --tmpdir gfabase_vg_walks_test_XXXXXX)

zstd -dc /tmp/GRCh38-f1-90-mc-mar13.chr21_chr22.zst | grep ^W | sort | sha256sum > "${TMPDIR}/original_walks" & pid=$!
zstd -dc /tmp/GRCh38-f1-90-mc-mar13.chr21_chr22.zst | time $gfabase load -o "${TMPDIR}/GRCh38-f1-90-mc-mar13.chr21_chr22.gfab" --compress 1
is "$?" "0" "gfabase load"
ls -lh "${TMPDIR}/GRCh38-f1-90-mc-mar13.chr21_chr22.gfab"

$gfabase view "${TMPDIR}/GRCh38-f1-90-mc-mar13.chr21_chr22.gfab" | grep ^W | sort | sha256sum > "${TMPDIR}/view_walks"
is "$?" "0" "gfabase view"

wait $pid
is "$(cat "${TMPDIR}/view_walks")" "$(cat "${TMPDIR}/original_walks")" "roundtrip walks"

gsql() {
    genomicsqlite "$1" "$2"| tail -n 1
}
export -f gsql

# TODO: chr21/chr22 breakdown and sanity-check vs original whole-genome file

rm -rf "$TMPDIR"
