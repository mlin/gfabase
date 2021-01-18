version development
# miniwdl workflow to map the segment sequences in a .gfab to a reference genome .fasta using
# winnowmap, generating a new mapped .gfab.
#
# Example invocation:
#
#   pip3 install miniwdl  # see https://github.com/chanzuckerberg/miniwdl
#
#   miniwdl run gfabase/workflows/gfab_winnowmap.wdl gfab=/path/to/myAssembly.gfab \
#       add_mappings.quality=30 --verbose --cfg gfabase/workflows/miniwdl.cfg
#
# Outputs the mapped .gfab as well as the paf.gz from winnowmap.
#
# On first use, takes extra time download the default GRCh38_no_alt_analysis_set reference, compile
# winnowmap, and generate the meryl/winnowmap repetitive k-mer list. The .cfg enables miniwdl's 
# cache so that subsequent invocations can reuse all of this work. The repetitive k-mer list is
# also output so that one can explicitly supply it to future invocations in case the cache isn't
# available. 
#
# * This WDL will need some work to use with runners other than miniwdl, as it currently uses some
#   convenience features specific to the latter.

workflow gfab_winnowmap {
    input {
        File gfab

        String reference_name = "GRCh38_no_alt_analysis_set"
        File reference_fasta_gz = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"

        # meryl/winnowmap repetitive k-mer list for the reference genome: if not supplied,
        # will be automatically computed (~20 minutes), and also output for future reuse.
        File? meryl_repetitive_kmers_gz

        String winnowmap_version = "2.0"
        File gfabase_exe = "https://github.com/mlin/gfabase/releases/download/v0.4.0/gfabase"
    }
 
    # Dockerfile recipe for building winnowmap executables (miniwdl-specific feature)
    Array[String] winnowmapDockerfile = [
        "FROM centos:7 AS builder",
        "RUN yum install -y -q which wget perl zlib-devel centos-release-scl",
        "RUN yum install -y -q devtoolset-8-gcc devtoolset-8-gcc-c++ devtoolset-8-make",
        "RUN mkdir /work",
        "WORKDIR /work",
        "RUN wget -nv https://github.com/marbl/Winnowmap/archive/v~{winnowmap_version}.tar.gz",
        "RUN tar zxf *.tar.gz",
        "RUN scl enable devtoolset-8 'make -j$(nproc) -C Winnowmap-*'",
        "RUN cp Winnowmap-*/bin/{meryl,winnowmap} .",
        "FROM ubuntu:20.04",
        "RUN apt-get -qq update && apt-get install -qq -y libgomp1 tabix",
        "COPY --from=builder /work/meryl /work/winnowmap /usr/local/bin/"
    ]

    # compute repetitive k-mers if needed
    if (!defined(meryl_repetitive_kmers_gz)) {
        call compute_meryl_repetitive_kmers {
            input:
            reference_name, reference_fasta_gz, winnowmapDockerfile
        }
    }
    File meryl_repetitive_kmers_gz_to_use = select_first([meryl_repetitive_kmers_gz, compute_meryl_repetitive_kmers.kmers_gz])

    # winnowmap segment sequences
    call winnowmap {
        input:
        meryl_repetitive_kmers_gz = meryl_repetitive_kmers_gz_to_use,
        gfab, reference_name, reference_fasta_gz, winnowmapDockerfile, gfabase_exe
    }

    # add mappings to .gfab
    call add_mappings {
        input:
        paf_gz = winnowmap.paf_gz,
        gfab, reference_name, gfabase_exe
    }

    output {
        File mapped_gfab = add_mappings.mapped_gfab
        File paf_gz = winnowmap.paf_gz
        File meryl_repetitive_kmers_gz_used = meryl_repetitive_kmers_gz_to_use
    }
}

task winnowmap {
    input {
        File gfab

        String reference_name
        File reference_fasta_gz
        File meryl_repetitive_kmers_gz

        String winnowmap_options = "-x asm5"

        Array[String] winnowmapDockerfile
        File gfabase_exe

        Int cpu = 16
        Int memGB = cpu
    }

    String output_name = basename(gfab, ".gfab") + ".~{reference_name}.paf.gz"

    command <<<
        set -euxo pipefail
        cp "~{gfabase_exe}" /usr/local/bin/gfabase
        chmod +x /usr/local/bin/gfabase

        bgzip -@ 4 -dc "~{reference_fasta_gz}" > /tmp/reference.fa & pid=$!
        bgzip -dc "~{meryl_repetitive_kmers_gz}" > /tmp/repetitive_kmers.txt
        # dump segments FASTA from .gfab
        gfabase view "~{gfab}" | awk '/^S/{print ">"$2"\n"$3}' | fold > /tmp/segments.fa
        wait $pid

        >&2 winnowmap --version
        winnowmap -W /tmp/repetitive_kmers.txt ~{winnowmap_options} -t ~{cpu} /tmp/reference.fa /tmp/segments.fa \
            | bgzip -@ 4 -c > "~{output_name}"
    >>>
    
    output {
        File paf_gz = output_name
    }

    runtime {
        cpu: cpu
        memory: "~{memGB} GB"
        inlineDockerfile: winnowmapDockerfile
    }
}

task compute_meryl_repetitive_kmers {
    input {
        String reference_name
        File reference_fasta_gz

        Int k = 19
        Float distinct = 0.9998

        Array[String] winnowmapDockerfile
    }

    String output_name = "~{reference_name}.repetitive.txt.gz"

    command <<<
        set -euxo pipefail
        bgzip -@ 4 -dc "~{reference_fasta_gz}" > /tmp/reference.fa
        >&2 meryl --version
        meryl count k=~{k} output merylDB /tmp/reference.fa
        meryl print greater-than distinct=~{distinct} merylDB | bgzip -@ 4 -c > "~{output_name}"
    >>>

    output {
        File kmers_gz = output_name
    }

    runtime {
        cpu: 16
        memory: "16 GB"
        inlineDockerfile: winnowmapDockerfile
    }
}

task add_mappings {
    input {
        File gfab

        String reference_name
        File paf_gz

        # mapQ & alignment block length filters
        Int quality = 0
        Int length = 0
        Boolean replace = true  # delete any existing mappings

        File gfabase_exe
    }

    String output_name = basename(gfab, ".gfab") + ".~{reference_name}.gfab"

    command <<<
        set -euxo pipefail
        cp "~{gfabase_exe}" /usr/local/bin/gfabase
        chmod +x /usr/local/bin/gfabase

        cp "~{gfab}" "~{output_name}"
        bgzip -@ 4 -dc "~{paf_gz}" | gfabase --verbose add-mappings "~{output_name}" \
            --quality ~{quality} \
            --length ~{length} \
            ~{if replace then "--replace" else ""}
    >>>

    output {
        File mapped_gfab = output_name
    }

    runtime {
        cpu: 4
        inlineDockerfile: [
            "FROM ubuntu:20.04",
            "RUN apt-get -qq update && apt-get install -qq -y libsqlite3-0 tabix"
        ]
    }
}
