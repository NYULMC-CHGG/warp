version 1.0

workflow GlimpseImputation {
    input {
        #list of reference files
        File reference_chunks = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/glimpse_split_reference/snpsIndels/glimpse_ref_panel"
        File ref_dict= "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/dragen_reference/imputed.dict"
        File input_vcf
        File input_vcf_index

        String output_basename

        #String docker = "us.gcr.io/broad-dsde-methods/glimpse:palantir-workflows_20c9de0"
        String docker = "us.gcr.io/broad-dsde-methods/glimpse:odelaneau_f310862"
        Int cpu_phase = 1
        Int mem_gb_phase = 10
        Int cpu_ligate = 1
        Int mem_gb_ligate = 5
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
            ref_dict = ref_dict,
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
        --threads ~{cpu} 
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
        runtime_minutes: 30
    }
}

task GlimpseLigate {
    input {
        Array[File] imputed_chunks
        Array[File] imputed_chunks_indices
        File ref_dict
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
        --output ligated.vcf \
        --threads ~{cpu}
        # reset header
        ## need to add fasta.dict and reset header from reference fasta
        bcftools view -h --no-version ligated.vcf > old_header.vcf
        java -jar /picard.jar UpdateVcfSequenceDictionary -I old_header.vcf --SD ~{ref_dict} -O new_header.vcf 
        bcftools reheader -h new_header.vcf -o ~{output_basename}.imputed.vcf ligated.vcf
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
        runtime_minutes: 120
    }

}
