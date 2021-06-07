# gfabase

`gfabase` is a command-line tool for indexed storage of [Graphical Fragment Assembly (GFA1)](https://github.com/GFA-spec/GFA-spec) data. It imports a .gfa file into a compressed **.gfab** file, from which it can later access subgraphs quickly (reading only the necessary parts), producing .gfa or .gfab. Beyond ID lookups, .gfab indexes the graph by mappings onto reference genome coordinates, facilitating navigation within *de novo* assemblies and pangenome reference graphs.

Effectively, .gfab is a new GFA-superset format with built-in compression and indexing. It is in fact a SQLite (+ [Genomics Extension](https://github.com/mlin/GenomicSQLite)) database populated with a [GFA1-like schema](src/schema/GFA1.sql), which programmers have the option to access directly, without requiring `gfabase` nor even a low-level parser for .gfa/.gfab.

### Quick start

Each [Release](https://github.com/mlin/gfabase/releases) includes prebuilt `gfabase` executables for Linux and macOS x86-64 hosts. The executable provides subcommands:

* `gfabase load -o my.gfab [my.gfa]`: create .gfab from a .gfa file (or pipe decompression through standard input)
* `gfabase view my.gfab`: dump back to .gfa (if standard output is a terminal, automatically pipes to `less -S`)
* `gfabase sub my.gfab SEGMENT/PATH/RANGE... [--view]`: query for a subgraph, producing either .gfa or .gfab
* `gfabase add-mappings my.gfab mappings.paf`: add index of reference genome mappings for GFA segments

The following quick example accesses a scaffold by its Path name in a [metaSPAdes](https://cab.spbu.ru/software/meta-spades/) assembly of simulated metagenomic reads from [Ye <em>et al.</em> (2019)](https://dx.doi.org/10.1016/j.cell.2019.07.010); it also uses [`zstd`](https://github.com/facebook/zstd) for decompression.

```bash
curl -L "https://github.com/mlin/gfabase/blob/main/test/data/atcc_staggered.assembly_graph_with_scaffolds.gfa.zst?raw=true" \
    | zstd -dc \
    | ./gfabase load -o atcc_staggered.metaspades.gfab

# extract a scaffold from the metagenome assembly (by GFA Path name)
./gfabase sub atcc_staggered.metaspades.gfab -o a_scaffold.gfab --path NODE_2_length_747618_cov_15.708553_3
# view GFA:
./gfabase view a_scaffold.gfab
# or in one command:
./gfabase sub atcc_staggered.metaspades.gfab --view --path NODE_2_length_747618_cov_15.708553_3
```

The following in-depth notebooks demonstrate human genome uses, also integrating with [Bandage](https://rrwick.github.io/Bandage/) for visualization:

1. **[Navigating a human *de novo* assembly](https://nbviewer.jupyter.org/github/mlin/gfabase/blob/main/notebooks/gfabaseAssemblyNavigation.ipynb)**
2. **[Slicing a pangenome reference graph](https://nbviewer.jupyter.org/github/mlin/gfabase/blob/main/notebooks/gfabasePangenomeGraph.ipynb)**

<img width="500" alt="index" src="https://user-images.githubusercontent.com/356550/105319466-fd571080-5b68-11eb-9422-a0b3b01c7056.png">

### Segment mappings

Adding `--range` to `gfabase sub` means the other command-line arguments are linear sequence ranges (chrom:begin-end) to be resolved to overlapping segments. This relies on mappings of each segment to its own linear coordinates, which `gfabase load` understands in two forms:

1. The [rGFA tags](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) `SN:Z` and `SO:i` are present *and* the segment sequence length is known (from given sequence or `LN:i`)
2. Segment tag `rr:Z` giving a browser-style range like `rr:Z:chr1:2,345-6,789`

Furthermore, `gfabase add-mappings my.gfab mappings.paf` adds mappings of segment sequences generated by [minimap2](https://github.com/lh3/minimap2) or a similar tool producing [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md). The .gfab is updated in-place, so make a backup copy if needed.

### Connected subgraphs

Adding `--connected` to `gfabase sub` expands the subgraph to include the complete connected component(s) associated with the specified segments.

That may be overkill, if we're only interested in the segments' immediate neighborhood. In that case, set `--biconnected 1` instead to extract the associated *bi*connected component(s), which naturally stops the expansion at "cutpoints" (which any end-to-end walk of the chromosome must traverse). Setting `--biconnected 2` or higher also includes biconnected components *adjacent* to those including the query segments.

<sup>The `--connected` and `--biconnected` expansions treat the segment graph as undirected. Therefore the extracted subgraphs include, but are not limited to, directed "superbubbles."</sup>

### Web access

`gfabase view` and `gfabase sub` can read .gfab http/https URLs directly. The web server must support HTTP GET range requests, and the content must be immutable. This is mainly useful to query for a small subgraph, especially with `--no-sequences`. On the other hand, a series of queries expected to traverse a large fraction of the graph will be better-served by downloading the whole file upfront.

Here's an example invocation to inspect the subgraph surrounding the HLA locus in a Shasta ONT assembly, remotely accessing a .gfab served by GitHub. (See the above-linked notebooks for details about the flags given.)

```
./gfabase sub \
    https://github.com/mlin/gfabase/releases/download/v0.5.0/shasta-HG002-Guppy-3.6.0-run4-UL.gfab \
    --view --cutpoints 2 --no-sequences --guess-ranges --range \
    chr6:29,700,000-29,950,000
```

**To publish a .gfab on the web,** it's helpful to first "defragment" the file using the [`genomicsqlite` command-line tool](https://mlin.github.io/GenomicSQLite/guide_db/#genomicsqlite-interactive-shell) made available by `pip3 install genomicsqlite` or `conda install -c mlin genomicsqlite`:

```
genomicsqlite my.gfab --compact --inner-page-KiB 64 --outer-page-KiB 2
```

...generating `my.gfab.compact`, a defragmented version that'll be more efficient to access. (`mv my.gfab.compact my.gfab` if so desired.)

### Building from source

![CI](https://github.com/mlin/gfabase/workflows/CI/badge.svg?branch=main)

[Install rust toolchain](https://rustup.rs/) and:

```bash
git clone https://github.com/mlin/gfabase.git
cd gfabase
./cargo build --release
```

Then find the executable `target/release/gfabase`.

<sup>
./cargo is a wrapper for cargo that generates Cargo.toml from Cargo.toml.in, filling in the crate version based on the git tag.
</sup>
