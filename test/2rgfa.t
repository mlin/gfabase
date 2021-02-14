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

# sub --cutpoints
gfabase sub "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab" \
    --range --cutpoints 2 -o "${TMPDIR}/megabase.gfab" chr22:11,000,000-12,000,000 --verbose
gfabase view "${TMPDIR}/megabase.gfab" > "${TMPDIR}/megabase.gfa"
is $(cat "${TMPDIR}/megabase.gfa" | grep "^S" | cut -f3 | LC_ALL=C sort | sha256sum | cut -f1 -d ' ') \
   "17d49156acd0ccad3452fb938b932234132a5d31f25ce92e7c655bff0628c654" "sub --cutpoints segments"
is $(cat "${TMPDIR}/megabase.gfa" | grep "^L" | wc -l | tr -d ' ') "56" "sub --cutpoints links"

gfabase sub "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab" \
    --range --cutpoints 2 --cutpoints-nt 10000 --view chr22:11,000,000-12,000,000 \
    > "${TMPDIR}/megabase10k.gfa"
is $(cat "${TMPDIR}/megabase10k.gfa" | grep "^S" | cut -f3 | LC_ALL=C sort | sha256sum | cut -f1 -d ' ') \
   "5ebd4b41cce8f4a99ea2b42d090315d43162f87482678ef3f1c70c65bcf5ae51" "sub --cutpoints-nt segments"
is $(cat "${TMPDIR}/megabase10k.gfa" | grep "^L" | wc -l | tr -d ' ') "59" "sub --cutpoints-nt links"

rm -rf "$TMPDIR"
