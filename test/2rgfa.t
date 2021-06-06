#!/bin/bash

# bash-tap tests for gfabase CLI ops
# requires cargo toolchain and:
#   apt install -y zstd && pip3 install --upgrade genomicsqlite

set -o pipefail

REPO="$(dirname "$0")/.."
cd "$REPO"
export BASH_TAP_ROOT=test/bash-tap
source test/bash-tap/bash-tap-bootstrap

plan tests 16

./cargo build --release
is "$?" "0" "cargo build"
PATH="$(pwd)/target/release:${PATH}"
gfabase version
is "$?" "0" "gfabase version"

if [[ -z $TMPDIR ]]; then
    TMPDIR=/tmp
fi
TMPDIR=$(mktemp -d "${TMPDIR}/gfabase_cli_test_XXXXXX")
export TMPDIR=$(realpath "$TMPDIR")

# extract test rGFA
zstd -dc test/data/GRCh38-20-0.10b.chr22_chrY.gfa.zst > "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfa"

# roundtrip it and check fidelity
gfabase load --compress 1 "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfa" -o "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab"
is "$?" "0" "gfabase load"
cat "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfa" \
    | gfabase --verbose load --compress 1 -o "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab"
is "$?" "0" "gfabase load stdin"

gfabase view "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab" > "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.roundtrip.gfa"
is "$?" "0" "gfabase view"

seq_digest=$(cat "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfa" | grep "^S" | cut -f3 | LC_ALL=C sort | sha256sum)
roundtrip_seq_digest=$(cat "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.roundtrip.gfa" | grep "^S" | cut -f3 | LC_ALL=C sort | sha256sum)
is "$roundtrip_seq_digest" "$seq_digest" "sequences roundtrip identical"

is `cat "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfa" | grep "^L" | wc -l | tr -d ' '` \
   `cat "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.roundtrip.gfa" | grep "^L" | wc -l | tr -d ' '` \
   "links count roundtrip identical"

# extract chr22 segments only
gfabase sub \
    "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab" -o "${TMPDIR}/GRCh38-20-0.10b.chr22only.gfab" \
    chr22:1-999999999 --range
chr22_digest=$(cat "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfa" | grep chr22 | cut -f3 | LC_ALL=C sort | sha256sum)
sub_chr22_digest=$(gfabase view "${TMPDIR}/GRCh38-20-0.10b.chr22only.gfab" | grep "^S" | cut -f3 | LC_ALL=C sort | sha256sum)
is "$sub_chr22_digest" "$chr22_digest" "chr22 reference segments"

# extract chr22 connected component
gfabase sub \
    "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab" -o "${TMPDIR}/GRCh38-20-0.10b.chr22.gfab" \
    chr22:1-999999999 --range --connected
is "$(gfabase view "${TMPDIR}/GRCh38-20-0.10b.chr22.gfab" | grep "^S" | wc -l | tr -d ' ')" "3319" "chr22-connected segments"
is "$(gfabase view "${TMPDIR}/GRCh38-20-0.10b.chr22.gfab" | grep "^L" | wc -l | tr -d ' ')" "4795" "chr22-connected links"

# sub --view to stream GFA directly
gfabase sub --view \
    "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab" \
    chr22:1-999999999 --range --connected \
    > "${TMPDIR}/GRCh38-20-0.10b.chr22.gfa"
is "$(grep "^S" "${TMPDIR}/GRCh38-20-0.10b.chr22.gfa" | wc -l | tr -d ' ')" "3319" "chr22-connected segments --view"
is "$(grep "^L" "${TMPDIR}/GRCh38-20-0.10b.chr22.gfa" | wc -l | tr -d ' ')" "4795" "chr22-connected links --view"

# sub --biconnected
gfabase sub "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab" \
    --range --biconnected --view chr22:11,000,000-12,000,000 \
    > "${TMPDIR}/megabase.gfa"
is $(cat "${TMPDIR}/megabase.gfa" | grep "^S" | cut -f3 | LC_ALL=C sort | sha256sum | cut -f1 -d ' ') \
   "4a7e30d1c601c809a47ed4dd430b21db7ad19aa0e37125b95d9a2deed950f51d" "sub --biconnected segments"
is $(cat "${TMPDIR}/megabase.gfa" | grep "^L" | wc -l | tr -d ' ') "21" "sub --biconnected links"

gfabase sub "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab" \
    --range --biconnected 2 --view chr22:11,000,000-12,000,000 \
    > "${TMPDIR}/megabase2.gfa"
is $(cat "${TMPDIR}/megabase2.gfa" | grep "^S" | cut -f3 | LC_ALL=C sort | sha256sum | cut -f1 -d ' ') \
   "40d135762d8991bca571e28615ad75c4d7bc76450922b9ceabfe4eeb2386240c" "sub --biconnected 2 segments"
is $(cat "${TMPDIR}/megabase2.gfa" | grep "^L" | wc -l | tr -d ' ') "57" "sub --biconnected 2 links"

rm -rf "$TMPDIR"
