version 1.0

workflow spades_atcc_staggered {
    input {
        # simulated FASTQs from Ye (2019) doi:10.1016/j.cell.2019.07.010, sha256sum:
        # dfbcfb16bc62a73dcb995e08718b748b48ca9b1455f08361746b482cbdeab307  atcc_staggered_1.fastq.bz2
        # 0717078a23754e3e8f5e54767a706f7651571af2a31e0a82b3fed1b9a4375573  atcc_staggered_2.fastq.bz2
        File fastq_bz2_1 = "gs://metax-bakeoff-2019/fastq/atcc_staggered_1.fastq.bz2"
        File fastq_bz2_2 = "gs://metax-bakeoff-2019/fastq/atcc_staggered_2.fastq.bz2"
        String output_name = "atcc_staggered"
    }

    scatter (bz2 in [fastq_bz2_1, fastq_bz2_2]) {
        call unbz2 {
            input:
            file_bz2 = bz2
        }
    }

    call metaspades {
        input:
        fastq_1 = unbz2.file[0],
        fastq_2 = unbz2.file[1],
        output_name = output_name
    }

    output {
        File assembly_graph_with_scaffolds_gfa_zst = metaspades.assembly_graph_with_scaffolds_gfa_zst
        File assembly_graph_after_simplification_gfa_zst = metaspades.assembly_graph_after_simplification_gfa_zst
        File strain_graph_gfa_zst = metaspades.strain_graph_gfa_zst
    }
}

task unbz2 {
    input {
        File file_bz2
    }

    String outname = basename(file_bz2, ".bz2")

    command <<<
        set -euxo pipefail
        apt-get -qq update && apt-get install -y bzip2
        bzip2 -dc "~{file_bz2}" > "~{outname}"
    >>>

    runtime {
        docker: "ubuntu:20.04"
    }

    output {
        File file = outname
    }
}

task metaspades {
    input {
        File fastq_1
        File fastq_2
        String output_name

        String version = "3.14.1"
        Int cpu = 16
        Int mem_GiB_per_cpu = 4
    }

    command <<<
        set -euxo pipefail
        apt-get -qq update && apt-get install -y wget zstd python3-pip time
        export TMPDIR=/tmp
        
        pushd "$TMPDIR"

        wget -nv "https://github.com/ablab/spades/releases/download/v~{version}/SPAdes-~{version}-Linux.tar.gz"
        tar zxf "SPAdes-~{version}-Linux.tar.gz"
        export PATH="$PATH:${TMPDIR}/SPAdes-~{version}-Linux/bin"

        >&2 time python3 "$(which metaspades.py)" \
            -1 "~{fastq_1}" -2 "~{fastq_2}" \
            -o spades_out --disable-gzip-output

        popd

        zstd -c -9 -T0 "${TMPDIR}/spades_out/assembly_graph_with_scaffolds.gfa" \
            > "~{output_name}.assembly_graph_with_scaffolds.gfa.zst"
        zstd -c -9 -T0 "${TMPDIR}/spades_out/assembly_graph_after_simplification.gfa" \
            > "~{output_name}.assembly_graph_after_simplification.gfa.zst"
        zstd -c -9 -T0 "${TMPDIR}/spades_out/strain_graph.gfa" \
            > "~{output_name}.strain_graph.gfa.zst"
    >>>

    runtime {
        docker: "ubuntu:20.04"
        cpu: cpu
        memory: "~{cpu*mem_GiB_per_cpu}GB"
    }

    output {
        File assembly_graph_with_scaffolds_gfa_zst = "~{output_name}.assembly_graph_with_scaffolds.gfa.zst"
        File assembly_graph_after_simplification_gfa_zst = "~{output_name}.assembly_graph_after_simplification.gfa.zst"
        File strain_graph_gfa_zst = "~{output_name}.strain_graph.gfa.zst"
    }
}
