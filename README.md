# gfabase

`gfabase` is a command-line tool for random-access storage of [Graphical Fragment Assembly (GFA)](https://github.com/GFA-spec/GFA-spec) data. It imports a .gfa file into a compressed **.gfab** file, from which it can later access subgraphs quickly (reading only the necessary parts), producing .gfa or .gfab. Beyond ID-based access, .gfab can index segments by their mappings to linear reference coordinates. This facilitates navigation within *de novo* assemblies or [rGFA reference graphs](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md), with quick access to subgraphs connected to linear coordinate ranges.

Effectively, .gfab is a new GFA-superset format with built-in compression and indexing. It is in fact a SQLite (+ [Genomics Extension](https://github.com/mlin/GenomicSQLite)) database populated with a [GFA-like schema](src/schema/GFA1.sql). Programmers have the option to access .gfab files directly using SQLite (+ Genomics Extension), without requiring `gfabase` nor even a low-level parser for .gfa/.gfab.

### Quick start

Each [Release](https://github.com/mlin/gfabase/releases) includes a prebuilt `gfabase` executable for modern Linux x86-64 hosts (with one caveat, shown there). The following examples also use the [`zstd` tool](https://github.com/facebook/zstd) for decompression.

```bash
./gfabase version

# stream two example GFA files into corresponding .gfab files:
# 1. metaSPAdes assembly of some simulated reads (from doi:10.1016/j.cell.2019.07.010)
curl -L "https://github.com/mlin/gfabase/blob/main/test/data/atcc_staggered.assembly_graph_with_scaffolds.gfa.zst?raw=true" \
    | zstd -dc \
    | ./gfabase load - atcc_staggered.metaspades.gfab
# 2. excerpt of a minigraph pangenome rGFA (from doi:10.1186/s13059-020-02168-z)
curl -L "https://github.com/mlin/gfabase/blob/main/test/data/GRCh38-20-0.10b.chr22_chrY.gfa.zst?raw=true" \
    | zstd -dc \
    | ./gfabase load - chr22_chrY.gfab

# extract a scaffold from the metagenome assembly (by GFA Path name)
./gfabase sub atcc_staggered.metaspades.gfab a_scaffold.gfab --path NODE_2_length_747618_cov_15.708553_3
# view GFA:
./gfabase view a_scaffold.gfab | less -S
# or in one command:
./gfabase sub atcc_staggered.metaspades.gfab - --view --path NODE_2_length_747618_cov_15.708553_3 | less -S

# from the rGFA graph, extract the subgraph starting from a certain megabase of chr22 and
# traversing up to 100,000 bp of linked segments
./gfabase sub chr22_chrY.gfab - --view --reference --radius-bp 100000 chr22:11,000,000-12,000,000 | less -S

# extract the whole connected component associated with chrY
./gfabase sub chr22_chrY.gfab chrY.gfab --reference --connected chrY:1-999,999,999
./gfabase view chrY.gfab | less -S
```

If we've also `pip3 install genomicsqlite`, then we can open a .gfab file in the SQLite interactive shell and poke around in SQL. (This does require up-to-date host SQLite and `sqlite3` shell.)

```
$ genomicsqlite chr22_chrY.gfab -readonly
> select count(1) from gfa1_segment_meta;
3833
> select sum(sequence_length) from gfa1_segment_meta where segment_id in
    (select distinct segment_id from gfa1_segment_mapping where _rowid_ in
       genomic_range_rowids('gfa1_segment_mapping', 'chr22', 11000000, 12000000)
     order by segment_id);
1441242
```

### Segment mappings

`gfabase load` currently understands two forms of linear mappings to make each segment discoverable by `gfabase sub --reference`,

1. The [rGFA tags](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) `SN:Z` and `SO:i` are present *and* the segment sequence length is known (from given sequence or `LN:i`)
2. Segment tag `rr:Z` giving a browser-style range like `rr:Z:chr1:2,345-6,789`

Notice that with `gfabase sub --reference --radius-bp N` you can usually extract the relevant subgraph given any segment with a "close enough" mapping.

Soon we plan to make it easy to source the ranges directly from a mapper run on the segment sequences. Please send ideas.

### Building from source

![CI](https://github.com/mlin/gfabase/workflows/CI/badge.svg?branch=main)

[Install rust toolchain](https://rustup.rs/) and:

```bash
git clone https://github.com/mlin/gfabase.git
cd gfabase
export RUSTFLAGS="-C link-args=-Wl,-rpath,\$ORIGIN"
cargo build --release
```

Then find the executable `target/release/gfabase`. <small>The `RUSTFLAGS` incantation makes it look for shared libraries alongside in the same directory before the usual system paths (useful for deploying a newer SQLite, as shown above).</small>
