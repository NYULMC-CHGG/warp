version 1.0

workflow GlimpseImputation {
    input {
        #list of reference files
        #File reference_chunks = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/glimpse_split_reference/glimpse_reference_chunks"
        File reference_chunks = "/gpfs/scratch/mccafj02/build_glimpse_reference_library/results/glimpse_ref_panel"
        File input_vcf
        File input_vcf_index

        String output_basename

        String docker = "us.gcr.io/broad-dsde-methods/glimpse:palantir-workflows_20c9de0"
        Int cpu_phase = 4
        Int mem_gb_phase = 5
        Int cpu_ligate = 4
        Int mem_gb_ligate = 4
    }

    scatter (reference_chunk in read_lines(reference_chunks)) {
        call GlimpsePhase {
            input:
                reference_chunk = reference_chunk,
                input_vcf = input_vcf,
                input_vcf_index = input_vcf_index,
                cpu = cpu_phase,
                mem_gb = mem_gb_phase,
                docker = docker

        }
    }
    call GlimpseLigate {
        input:
            imputed_chunks = GlimpsePhase.imputed_vcf,
            imputed_chunks_indices = GlimpsePhase.imputed_vcf_index,
            output_basename = output_basename,
            docker = docker,
            cpu = cpu_ligate,
            mem_gb = mem_gb_ligate
    }

    output {
        File imputed_vcf = GlimpseLigate.imputed_vcf
        File imputed_vcf_index = GlimpseLigate.imputed_vcf_index
    }

}

task GlimpsePhase {
    input {
        File input_vcf
        File input_vcf_index
        File reference_chunk
        Int mem_gb
        Int cpu
        String docker
        Int max_retries = 3
    }
    command <<<
        set -euo pipefail

        /bin/GLIMPSE2_phase \
        --input-gl ~{input_vcf} \
        --reference ~{reference_chunk} \
        --output phase_output.bcf \
        --threads ~{cpu} \
    >>>
    output {
        File imputed_vcf = "phase_output.bcf"
        File imputed_vcf_index = "phase_output.bcf.csi"
    }
    runtime {
        docker : docker
        memory : mem_gb +" GiB"
        cpu: cpu
        maxRetries: max_retries
        runtime_minutes: 10
    }
}

task GlimpseLigate {
    input {
        Array[File] imputed_chunks
        Array[File] imputed_chunks_indices
        String output_basename
        Int mem_gb
        Int cpu
        Int max_retries = 3
        String docker
        
        
    }
    command <<<
        set -xeuo pipefail
       
        /bin/GLIMPSE2_ligate \
        --input ~{write_lines(imputed_chunks)} \
        --output ~{output_basename}.imputed.vcf \
        --threads ~{cpu}

        bgzip ~{output_basename}.imputed.vcf
        tabix ~{output_basename}.imputed.vcf.gz
    >>>
    output {
        File imputed_vcf = "~{output_basename}.imputed.vcf.gz"
        File imputed_vcf_index = "~{output_basename}.imputed.vcf.gz.tbi"
        
    }
    runtime {
        docker : docker
        memory : mem_gb +" GiB"
        cpu: cpu
        maxRetries: max_retries
        runtime_minutes: 20
    }

}
