version 1.0

workflow Glimpse2Impute {
    input {
        File ref_binary = "/gpfs/scratch/mccafj02/build_glimpse_reference_library/results/glimpse_ref.tar.gz"
        File input_vcf
        File input_vcf_index
        String output_basename
        String docker = "us.gcr.io/broad-dsde-methods/glimpse:palantir-workflows_20c9de0"
        Int cpu = 5
        Int mem_gb = 100
        
    }

    call Impute {
        input:
            ref_binary = ref_binary,
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            output_basename = output_basename,
            docker = docker,
            cpu = cpu,
            mem_gb = mem_gb
    }
    output{
        File imputed_vcf = Impute.imputed_vcf
        File imputed_vcf_index = Impute.imputed_vcf_index
    }
}

task Impute{
    input{
        File ref_binary
        File input_vcf
        File input_vcf_index
        String output_basename
        String docker
        Int cpu
        Int mem_gb
    }
    command <<<
        set -euo pipefail

        tar -xzf ~{ref_binary}
        mkdir phase_output
        CHUNK=0
        for ref_chunk in ref_binary/*.bin
        do
        /bin/GLIMPSE2_phase \
        --input-gl ~{input_vcf} \
        --reference ${ref_chunk} \
        --output phase_output/${CHUNK}_phase_output.bcf \
        --threads ~{cpu} 
        (( CHUNK++ ))    
        done

        touch imputed_chunks
        touch temp_file
        CHUNK=0
        for chunk in phase_output/*
        do 
            echo $chunk >> temp_file
        done
        grep -v ".csi" temp_file > imputed_chunks

        /bin/GLIMPSE2_ligate \
        --input imputed_chunks \
        --output ~{output_basename}.imputed.vcf \
        --threads ~{cpu}

        bgzip ~{output_basename}.imputed.vcf
        tabix ~{output_basename}.imputed.vcf.gz
             
    >>>
    output{
        File imputed_vcf = "~{output_basename}.imputed.vcf.gz"
        File imputed_vcf_index = "~{output_basename}.imputed.vcf.gz.tbi"
    }
    runtime{
        docker : docker
        memory : mem_gb +" GiB"
        cpu: cpu
        runtime_minutes: 600
    }
}