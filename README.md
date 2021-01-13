# gfabase

`gfabase` is a command-line tool for random-access storage of [Graphical Fragment Assembly (GFA1)](https://github.com/GFA-spec/GFA-spec) data. It imports a .gfa file into a compressed **.gfab** file, from which it can later access subgraphs quickly (reading only the necessary parts), producing .gfa or .gfab. Beyond ID lookups, .gfab provides quick access to subgraphs by their mappings to linear reference coordinates, facilitating navigation within *de novo* assemblies or [rGFA reference graphs](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md).

Effectively, .gfab is a new GFA-superset format with built-in compression and indexing. It is in fact a SQLite (+ [Genomics Extension](https://github.com/mlin/GenomicSQLite)) database populated with a [GFA1-like schema](src/schema/GFA1.sql). Programmers have the option to access .gfab files directly using SQLite (+ Genomics Extension), without requiring `gfabase` nor even a low-level parser for .gfa/.gfab.

### Quick start

Each [Release](https://github.com/mlin/gfabase/releases) includes a prebuilt `gfabase` executable for modern Linux x86-64 hosts (with one caveat, shown there). The following examples also use the [`zstd` tool](https://github.com/facebook/zstd) for decompression.

**Example 1. metaSPAdes assembly of simulated metagenomic reads: access scaffold by Path name**

```bash
curl -L "https://github.com/mlin/gfabase/blob/main/test/data/atcc_staggered.assembly_graph_with_scaffolds.gfa.zst?raw=true" \
    | zstd -dc \
    | ./gfabase load - atcc_staggered.metaspades.gfab

# extract a scaffold from the metagenome assembly (by GFA Path name)
./gfabase sub atcc_staggered.metaspades.gfab a_scaffold.gfab --path NODE_2_length_747618_cov_15.708553_3
# view GFA:
./gfabase view a_scaffold.gfab | less -S
# or in one command:
./gfabase sub atcc_staggered.metaspades.gfab - --view --path NODE_2_length_747618_cov_15.708553_3 | less -S
```

<sup>Simulated reads from [Ye <em>et al.</em> (2019)](https://dx.doi.org/10.1016/j.cell.2019.07.010)</sup>

**Example 2. Human pangenome rGFA: access connected component & subgraph**

```bash
curl -L "https://github.com/mlin/gfabase/blob/main/test/data/GRCh38-20-0.10b.chr22_chrY.gfa.zst?raw=true" \
    | zstd -dc \
    | ./gfabase --verbose load - chr22_chrY.gfab

# from the rGFA graph, extract the whole connected component associated with chrY
./gfabase sub chr22_chrY.gfab chrY.gfab --range --connected chrY:1-999,999,999
./gfabase view chrY.gfab | less -S

# starting from the segments constituting a megabase of chr22, expand to the connected subgraph
# without crossing "cutpoints" -- defined below
./gfabase sub chr22_chrY.gfab - --view --range --cutpoints 1 chr22:11,000,000-12,000,000 | less -S
```
<sup>rGFA excerpted from from [Li, Feng, & Chu (2020)](https://dx.doi.org/10.1186/s13059-020-02168-z)</sup>


If we've also `pip3 install genomicsqlite`, then we can open a .gfab file in the SQLite interactive shell and poke around in SQL.

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

`gfabase load` currently understands two forms of linear mappings to make each segment discoverable by `gfabase sub --range`,

1. The [rGFA tags](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) `SN:Z` and `SO:i` are present *and* the segment sequence length is known (from given sequence or `LN:i`)
2. Segment tag `rr:Z` giving a browser-style range like `rr:Z:chr1:2,345-6,789`

Soon we plan to make it easy to source the ranges directly from a mapper run on the segment sequences. Please send ideas.

### Subgraph cutpoints

Using `gfabase sub` to extract a set or range of segments, we often want to get their "neighborhood" too, without loading in the whole connected component. The `.gfab` has an index of *cutpoints*, which are (informally) segments with no possible detours around them to walk from one end of the chromosome to the other. These are natural boundaries for local subgraph extraction:

* `--cutpoints 1` finds the subgraph connected to the command-line segments *without* crossing any cutpoint
* `--cutpoints N` finds the subgraph connected while crossing *at most N-1* cutpoints
* `--cutpoints N --cutpoints-nt L` *id.* but only cutpoint segments at least *L* nucleotides long count toward *N*

<sup>Fine print: you may get a bit more than you asked for; because the cutpoints are precomputed from an undirected segment graph, some repeat and inversion motifs are considered possible detours that really aren't.</sup>

### Building from source

![CI](https://github.com/mlin/gfabase/workflows/CI/badge.svg?branch=main)

[Install rust toolchain](https://rustup.rs/) and:

```bash
git clone https://github.com/mlin/gfabase.git
cd gfabase
export RUSTFLAGS="-C link-args=-Wl,-rpath,\$ORIGIN"
./cargo build --release
```

Then find the executable `target/release/gfabase`.

<sup>
1. ./cargo is a wrapper for cargo that generates Cargo.toml from Cargo.toml.in, filling in the crate version based on the git tag.
</sup>
<br/>
<sup>
2. The RUSTFLAGS setting makes it look for shared libraries alongside in the same directory before the usual system paths (useful for bundling SQLite when needed).
</sup>
