version 1.0

# imports
import "../structs/dna_seq/DNASeqStructs.wdl"
import "../tasks/broad/GermlineVariantDiscovery.wdl" as Calling
import "../pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl" as ToWGSgermline

# workflow
workflow Deciphehr{

    #inputs
    input{
        File dnaSeqSSRef = "/gpfs/data/deciphEHRlab/pipeline/warp/structs/dna_seq/DNASeqSingleSampleReferences.json"
        File dragMapRef = "/gpfs/data/deciphEHRlab/pipeline/warp/structs/dna_seq/DragmapReference.json"
        File papiSettings = "/gpfs/data/deciphEHRlab/pipeline/warp/structs/dna_seq/PapiSettings.json"
        File variantCallScatter = "/gpfs/data/deciphEHRlab/pipeline/warp/structs/dna_seq/VariantCallingScatterSettings.json"
        SampleFastq sample_fastq
        File wgs_coverage_interval_list = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/v0/wgs_coverage_regions.hg38.interval_list"
        File allele_vcf = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/forceCall/variant_1000g_key.vcf.gz"
        File allele_index = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/forceCall/variant_1000g_key.vcf.gz.tbi"
        Boolean isXY = false

    }

    ## set up Structs
    ## read in data structures config files
    DNASeqSingleSampleReferences references = read_json(dnaSeqSSRef)
    DragmapReference dragmap_reference = read_json(dragMapRef)
    VariantCallingScatterSettings scatter_settings = read_json(variantCallScatter)
    PapiSettings papi_settings = read_json(papiSettings)

    ## convert Fastq to UBAM
    call ConvertFastqToUbam as FastqToUbam{
        input:
            sample_fastq = sample_fastq
    }
    ## get sample struct in json format
    call SampleStruct{
        input:
            sample_fastq = sample_fastq,
            ubam = FastqToUbam.ubam
    }

    # set up Struct for WGS workflow
    SampleAndUnmappedBams sample_and_unmapped_bams = read_json(SampleStruct.fjson)
    
    ## call WholeGermlineSingleSample workflow
    call ToWGSgermline.WholeGenomeGermlineSingleSample as wgs{
        input:
            sample_and_unmapped_bams = sample_and_unmapped_bams,
            references = references,
            wgs_coverage_interval_list = wgs_coverage_interval_list,
            dragmap_reference = dragmap_reference,
            scatter_settings = scatter_settings,
            papi_settings = papi_settings,
            dragen_functional_equivalence_mode = true,
            provide_bam_output = true,
            allow_empty_ref_alt = true,
            allele_vcf = allele_vcf,
            allele_index = allele_index

    }

    
    ## call XY typing
    call XYtyping {
        input: 
            samid = sample_fastq.sample_name,
            align_bam = wgs.output_bam,
            align_index = wgs.output_bam_index
    }
    #isXY = read_boolean(XYtyping.isXY)
    ## if XY
        if( read_boolean(XYtyping.isXY) ){
            ## call variants on X and Y par/non-par regions
            call XYregions{}
            Array[Region] regions = read_json(XYregions.json_regions)
            scatter( x in regions){
                call HaplotypeCallerXY{
                input: ## NEED TO FIGURE OUT THE REFERENCE AND INPUT FILES
                    fasta_ref = references.reference_fasta.ref_fasta,
                    fasta_index = references.reference_fasta.ref_fasta_index,
                    fasta_dict = references.reference_fasta.ref_dict,
                    input_bam = wgs.output_bam,
                    input_bam_index = wgs.output_bam_index,
                    region = x.region,
                    ploidy = x.ploidy,
                    mem_size_mb = 15000,
                    cpu = 10,
                    id = x.id,
                    allele_vcf = allele_vcf,
                    allele_index = allele_index

                }
            }
            #gather vcf outputs and merge into a single gvcf file
            Array[File] vcfs_to_merge = flatten(select_all([HaplotypeCallerXY.vcf]))
            Array[File] vcf_indices_to_merge = flatten(select_all([HaplotypeCallerXY.vcf_index]))
            call Calling.MergeVCFs as XYVCFs {
                input:
                    input_vcfs = vcfs_to_merge,
                    input_vcfs_indexes = vcf_indices_to_merge,
                    output_vcf_name = sample_fastq.sample_name+"_XY.g.vcf.gz",
            }
            # call FinalVCF {
            # input:
            #     fasta_ref = references.reference_fasta.ref_fasta,
            #     fasta_index = references.reference_fasta.ref_fasta_index,
            #     fasta_dict = references.reference_fasta.ref_dict,
            #     xy_vcfs = MergeVCFs.output_vcf,
            #     input_vcfs_indexes = MergeVCFs.output_vcf_index,
            #     autosome_vcf = wgs.output_vcf,
            #     autosome_vcf_index = wgs.output_vcf_index,
            #     output_name = sample_fastq.sample_name+"_XY_hard-filtered.g.vcf.gz"
            # }
        
        }
    # if( read_boolean(XYtyping.isXY) ){
    #     File output_vcf = FinalVCF.vcf
    #     File output_vcf_index = FinalVCF.vcf_index
    # }
    # if( !read_boolean(XYtyping.isXY)){
    #     File output_vcf = wgs.output_vcf
    #     File output_vcf_index = wgs.output_vcf_index
    # }
    

    ## output everything
    output{
       

        File? cross_check_fingerprints_metrics = wgs.cross_check_fingerprints_metrics

        

        File duplicate_metrics = wgs.duplicate_metrics
        File? output_bqsr_reports = wgs.output_bqsr_reports

        File gvcf_summary_metrics = wgs.gvcf_summary_metrics
        File gvcf_detail_metrics = wgs.gvcf_detail_metrics

        File? output_bam = wgs.output_bam
        File? output_bam_index = wgs.output_bam_index

        File output_vcf = wgs.output_vcf
        File output_vcf_index = wgs.output_vcf_index

        File xytyping = XYtyping.output_ratio
        File? xyvcf = XYVCFs.output_vcf
        File? xyvcf_index = XYVCFs.output_vcf_index
        
    }
    meta {
    allowNestedInputs: true
  }
    
}
task SampleStruct{
    input{
        SampleFastq sample_fastq
        File ubam
    }
    command <<<
        JSON_STRING='{"sample_name":"~{sample_fastq.sample_name}","base_file_name":"~{sample_fastq.base_file_name}","flowcell_unmapped_bams":["~{ubam}"],"final_gvcf_base_name":"~{sample_fastq.final_gvcf_base_name}","unmapped_bam_suffix":"~{sample_fastq.unmapped_bam_suffix}"}'
        echo $JSON_STRING > fout.json

    >>>
    output{
        File fjson = "fout.json"
    }
    runtime{
        cpu: 1
        memory: "1000 MiB"
        runtime_minutes: 2
    }
}
task ConvertFastqToUbam {
    input{
        SampleFastq sample_fastq
    }

    command <<<

       java -Xmx8G -jar /picard/picard.jar FastqToSam \
       FASTQ=~{sample_fastq.r1} \
       FASTQ2=~{sample_fastq.r2} \
       OUTPUT=~{sample_fastq.sample_name}_unmapped.bam \
       READ_GROUP_NAME="HLYLLDRXX" \
       SAMPLE_NAME=~{sample_fastq.sample_name} \
       LIBRARY_NAME="Solexa-272222" \
       PLATFORM_UNIT="HLYLLDRXX.2" \
       PLATFORM="illumina" \
       SEQUENCING_CENTER="NYU" \
       RUN_DATE="2014-08-20T00:00:00-0400"

    >>>

    output{
        File ubam = "~{sample_fastq.sample_name}_unmapped.bam"

    }
    runtime{
        docker: "us.gcr.io/broad-gotc-prod/dragmap:1.1.2-1.2.1-2.26.10-1.11-1643839530"
        cpu: 1
        memory: "9000 MiB"
        runtime_minutes: 150
    }
}
task XYtyping {
    input{
        String samid
        File? align_bam
        File? align_index
    }
    command <<<

        x_map=$(samtools idxstats ~{align_bam} | grep "chrX\s" | cut -f 3)
        x_len=$(samtools idxstats ~{align_bam} | grep "chrX\s" | cut -f 2)
        x_cov=$(echo "scale=10; ${x_map}/${x_len}" | bc)
        y_map=$(samtools idxstats ~{align_bam}| grep "chrY\s" | cut -f 3)
        y_len=$(samtools idxstats ~{align_bam} | grep "chrY\s" | cut -f 2)
        y_cov=$(echo "scale=10; ${y_map}/${y_len}" | bc)
        if(( $(echo "${y_cov} == 0" | bc -l) )); then
            ratio=${x_cov}
        else
            ratio=$(echo "scale=10; ${x_cov}/${y_cov}" | bc)
        fi

        if (( $(echo "${ratio} > 4.00" | bc -l) )); then
            sex="false"
        else
            sex="true"
        fi

        echo "sample xy_ratio isXY" | tr ' ' , > ~{samid}.XYratio.csv
        echo ~{samid} ${ratio} ${sex} | tr ' ' , >> ~{samid}.XYratio.csv
        echo ${sex} > isXY.txt

    >>>
    output{
        File output_ratio="~{samid}.XYratio.csv"
        File isXY="isXY.txt"
        #String xy=read_lines(stdout())
    }
    runtime{
        docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
        cpu: 1
        memory: "800 MiB"
        runtime_minutes: 30
    }
}
task XYregions{
    input{
       
    }
    command<<<
        JSON_STRING='[{"region":"chrX:1-10000","ploidy":"1","id":"r1"},{"region":"chrX:10001-2781479","ploidy":"2","id":"r2"},
                      {"region":"chrX:2781480-155701382","ploidy":"1","id":"r3"},{"region":"chrX:155701383-156030895","ploidy":"2","id":"r4"},
                      {"region":"chrX:156030896-156040895","ploidy":"1","id":"r5"},{"region":"chrY","ploidy":"1","id":"r6"}]'
        echo $JSON_STRING > fout.json
        
    >>>
    output{
        File json_regions = "fout.json"
    }
    runtime{
        cpu: 1
        memory: "800 MiB"
        runtime_minutes: 2
    }
}
task HaplotypeCallerXY{
    input{
        File fasta_ref
        File fasta_index
        File fasta_dict
        File? input_bam
        File? input_bam_index
        String region
        String ploidy
        String id
        Int mem_size_mb
        Int cpu
        File allele_vcf
        File allele_index
        
        
    }
    command <<<
        set -e 
        let java_memory_size_mb=$((~{mem_size_mb} - 1024))
        gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        -R ~{fasta_ref} \
        -I ~{input_bam} \
        -L ~{region} \
        -O ~{id}.vcf.gz \
        --alleles ~{allele_vcf} \
        --native-pair-hmm-threads ~{cpu} \
        -ploidy ~{ploidy} \
        -ERC GVCF \
        -G StandardAnnotation \
        -G StandardHCAnnotation \
        -G AS_StandardAnnotation \
        -RF OverclippedReadFilter \
        --smith-waterman FASTEST_AVAILABLE
    >>>
    output{
        File vcf = "~{id}.vcf.gz"
        File vcf_index = "~{id}.vcf.gz.tbi"
    }
    runtime{
        docker:"us.gcr.io/broad-gatk/gatk:4.3.0.0"
        cpu: cpu
        runtime_minutes: 900
        memory: "~{mem_size_mb} MiB"
    }
}
task FinalVCF {
    input{
        File fasta_ref
        File fasta_index
        File fasta_dict
        File xy_vcfs
        File autosome_vcf
        File input_vcfs_indexes
        File autosome_vcf_index
        String output_name
        String noXY = "noXY.vcf.gz"
    }
    
    command <<<

        gatk --java-options "-Xms5000m -Xmx5000m" \
        SelectVariants \
        -R ~{fasta_ref} \
        -V ~{autosome_vcf} \
        -XL chrY \
        -XL chrX \
        -O ~{noXY}

        gatk --java-options "-Xms5000m -Xmx5000m" \
        MergeVcfs \
        -I ~{xy_vcfs} \
        -I ~{noXY} \
        -O ~{output_name}
    >>>
    output{
        File vcf = "~{output_name}"
        File vcf_index = "~{output_name}.tbi"
    }
    runtime{
        docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
        cpu: 1
        memory: "6000 MiB"
        runtime_minutes: 90
    }
}