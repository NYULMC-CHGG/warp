version 1.0

#imports
import "../structs/dna_seq/DNASeqStructs.wdl"
import "../tasks/broad/UnmappedBamToAlignedBam.wdl" as ToBam
import "../tasks/broad/AggregatedBamQC.wdl" as AggregatedQC
import "../tasks/broad/Qc.wdl" as QC
import "../tasks/broad/BamToCram.wdl" as ToCram
import "../tasks/broad/Utilities.wdl" as Utilities
import "../pipelines/broad/dna_seq/germline/variant_calling/VariantCalling.wdl" as ToGvcf
import "../tasks/broad/GermlineVariantDiscovery.wdl" as Calling

#workflow
workflow Production{

    #inputs
    input{
        File dnaSeqSSRef = "/gpfs/data/deciphEHRlab/pipeline/warp/structs/dna_seq/DNASeqSingleSampleReferences.json"
        File dragMapRef = "/gpfs/data/deciphEHRlab/pipeline/warp/structs/dna_seq/DragmapReference.json"
        File papiSettings = "/gpfs/data/deciphEHRlab/pipeline/warp/structs/dna_seq/PapiSettings.json"
        File variantCallScatter = "/gpfs/data/deciphEHRlab/pipeline/warp/structs/dna_seq/VariantCallingScatterSettings.json"
        SampleFastq sample_fastq
        #File wgs_coverage_interval_list = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/v0/wgs_coverage_regions.hg38.interval_list"
        File allele_vcf = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/forceCall/variant_1000g_key.vcf.gz"
        File allele_index = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/forceCall/variant_1000g_key.vcf.gz.tbi"
        Boolean provide_bam_output = true
        Boolean use_gatk3_haplotype_caller = false

        Boolean dragen_functional_equivalence_mode = true
        Boolean dragen_maximum_quality_mode = false

        Boolean run_dragen_mode_variant_calling = false
        Boolean use_spanning_event_genotyping = true
        Boolean unmap_contaminant_reads = true
        Boolean perform_bqsr = true
        Boolean use_bwa_mem = true
        Boolean allow_empty_ref_alt = false
        Boolean use_dragen_hard_filtering = false
        Boolean skip_reblocking = true
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
            sample_fastq = sample_fastq,
    }
    ## get sample struct in json format
    call SampleStruct{
        input:
            sample_fastq = sample_fastq,
            ubam = FastqToUbam.ubam
    }

    # set up Struct for WGS workflow
    SampleAndUnmappedBams sample_and_unmapped_bams = read_json(SampleStruct.fjson)
    
    # Set DRAGEN-related arguments according to the preset arguments
  Boolean run_dragen_mode_variant_calling_ = if (dragen_functional_equivalence_mode || dragen_maximum_quality_mode) then true else run_dragen_mode_variant_calling
  Boolean use_spanning_event_genotyping_ = if dragen_functional_equivalence_mode then false else (if dragen_maximum_quality_mode then true else use_spanning_event_genotyping)
  Boolean unmap_contaminant_reads_ = if dragen_functional_equivalence_mode then false else (if dragen_maximum_quality_mode then true else unmap_contaminant_reads)
  Boolean perform_bqsr_ = if (dragen_functional_equivalence_mode || dragen_maximum_quality_mode) then false else perform_bqsr
  Boolean use_bwa_mem_ = if (dragen_functional_equivalence_mode || dragen_maximum_quality_mode) then false else use_bwa_mem
  Boolean use_gatk3_haplotype_caller_ = if (dragen_functional_equivalence_mode || dragen_maximum_quality_mode) then false else use_gatk3_haplotype_caller
  Boolean use_dragen_hard_filtering_ = if (dragen_functional_equivalence_mode || dragen_maximum_quality_mode) then true else use_dragen_hard_filtering

   # Not overridable:
  Float lod_threshold = -20.0
  String cross_check_fingerprints_by = "READGROUP"
  String recalibrated_bam_basename = sample_and_unmapped_bams.base_file_name + ".aligned.duplicates_marked.recalibrated"

  String final_gvcf_base_name = select_first([sample_and_unmapped_bams.final_gvcf_base_name, sample_and_unmapped_bams.base_file_name])

    #align sample
    call ToBam.UnmappedBamToAlignedBam{
        input:
            sample_and_unmapped_bams = sample_and_unmapped_bams,
            references = references,
            dragmap_reference = dragmap_reference,
            papi_settings = papi_settings,
            contamination_sites_ud = references.contamination_sites_ud,
            contamination_sites_bed = references.contamination_sites_bed,
             contamination_sites_mu = references.contamination_sites_mu,

             cross_check_fingerprints_by = cross_check_fingerprints_by,
             haplotype_database_file     = references.haplotype_database_file,
             lod_threshold               = lod_threshold,
             recalibrated_bam_basename   = recalibrated_bam_basename,
             perform_bqsr                = perform_bqsr_,
             use_bwa_mem                 = use_bwa_mem_,
             unmap_contaminant_reads     = unmap_contaminant_reads_,
             allow_empty_ref_alt         = allow_empty_ref_alt
    }
    ## call XY typing
    call XYtyping {
        input: 
            samid = sample_fastq.sample_name,
            align_bam = UnmappedBamToAlignedBam.output_bam,
            align_index = UnmappedBamToAlignedBam.output_bam_index
    }
    # XX sample
    if ( !read_boolean(XYtyping.isXY) ){
        call ToGvcf.VariantCalling as XXGvcf {
          input:
            run_dragen_mode_variant_calling = run_dragen_mode_variant_calling_,
            use_spanning_event_genotyping = use_spanning_event_genotyping_,
            calling_interval_list = references.calling_interval_list,
            evaluation_interval_list = references.evaluation_interval_list,
            haplotype_scatter_count = scatter_settings.haplotype_scatter_count,
            break_bands_at_multiples_of = scatter_settings.break_bands_at_multiples_of,
            input_bam = UnmappedBamToAlignedBam.output_bam,
            input_bam_index = UnmappedBamToAlignedBam.output_bam_index,
            ref_fasta = references.reference_fasta.ref_fasta,
            ref_fasta_index = references.reference_fasta.ref_fasta_index,
            ref_dict = references.reference_fasta.ref_dict,
            ref_str = references.reference_fasta.ref_str,
            dbsnp_vcf = references.dbsnp_vcf,
            dbsnp_vcf_index = references.dbsnp_vcf_index,
            base_file_name = sample_and_unmapped_bams.base_file_name,
            final_vcf_base_name = final_gvcf_base_name,
            agg_preemptible_tries = papi_settings.agg_preemptible_tries,
            use_gatk3_haplotype_caller = use_gatk3_haplotype_caller_,
            use_dragen_hard_filtering = use_dragen_hard_filtering_,
            skip_reblocking = skip_reblocking,
            allele_vcf = allele_vcf,
            allele_index = allele_index
      
        }


    }
    if ( read_boolean(XYtyping.isXY) ){
        call ToGvcf.VariantCalling as XYGvcf {
          input:
            run_dragen_mode_variant_calling = run_dragen_mode_variant_calling_,
            use_spanning_event_genotyping = use_spanning_event_genotyping_,
            calling_interval_list = references.calling_interval_list,
            evaluation_interval_list = references.evaluation_interval_list,
            haplotype_scatter_count = scatter_settings.haplotype_scatter_count,
            break_bands_at_multiples_of = scatter_settings.break_bands_at_multiples_of,
            input_bam = UnmappedBamToAlignedBam.output_bam,
            input_bam_index = UnmappedBamToAlignedBam.output_bam_index,
            ref_fasta = references.reference_fasta.ref_fasta,
            ref_fasta_index = references.reference_fasta.ref_fasta_index,
            ref_dict = references.reference_fasta.ref_dict,
            ref_str = references.reference_fasta.ref_str,
            dbsnp_vcf = references.dbsnp_vcf,
            dbsnp_vcf_index = references.dbsnp_vcf_index,
            base_file_name = sample_and_unmapped_bams.base_file_name,
            final_vcf_base_name = final_gvcf_base_name,
            agg_preemptible_tries = papi_settings.agg_preemptible_tries,
            use_gatk3_haplotype_caller = use_gatk3_haplotype_caller_,
            use_dragen_hard_filtering = use_dragen_hard_filtering_,
            skip_reblocking = skip_reblocking,
            allele_vcf = allele_vcf,
            allele_index = allele_index
      
        }
        ## call variants on X and Y par/non-par regions
        call XYregions{}
        Array[Region] regions = read_json(XYregions.json_regions)
        scatter( x in regions){
                call HaplotypeCallerXY{
                input: ## NEED TO FIGURE OUT THE REFERENCE AND INPUT FILES
                    fasta_ref = references.reference_fasta.ref_fasta,
                    fasta_index = references.reference_fasta.ref_fasta_index,
                    fasta_dict = references.reference_fasta.ref_dict,
                    input_bam = UnmappedBamToAlignedBam.output_bam,
                    input_bam_index = UnmappedBamToAlignedBam.output_bam_index,
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
            call CorrectedVCF {
            input:
                fasta_ref = references.reference_fasta.ref_fasta,
                fasta_index = references.reference_fasta.ref_fasta_index,
                fasta_dict = references.reference_fasta.ref_dict,
                xy_vcfs = XYVCFs.output_vcf,
                input_vcfs_indexes = XYVCFs.output_vcf_index,
                autosome_vcf = XYGvcf.output_vcf,
                autosome_vcf_index = XYGvcf.output_vcf_index,
                output_name = sample_fastq.sample_name+".hard-filtered.g.vcf.gz"
            }

    }

    if (provide_bam_output) {
    File provided_output_bam = UnmappedBamToAlignedBam.output_bam
    File provided_output_bam_index = UnmappedBamToAlignedBam.output_bam_index
    }

    #set up outputs
    File vcf_output = select_first([XXGvcf.output_vcf, CorrectedVCF.output_vcf])
    File vcf_output_index = select_first([XXGvcf.output_vcf_index, CorrectedVCF.output_vcf_index])
    File summary_metrics = select_first([XXGvcf.vcf_summary_metrics, XYGvcf.vcf_summary_metrics])
    File detail_metrics = select_first([XXGvcf.vcf_detail_metrics, XYGvcf.vcf_detail_metrics])
  

## output everything
    output{
       

        File? cross_check_fingerprints_metrics = UnmappedBamToAlignedBam.cross_check_fingerprints_metrics
        File duplicate_metrics = UnmappedBamToAlignedBam.duplicate_metrics
        File? output_bqsr_reports = UnmappedBamToAlignedBam.output_bqsr_reports

        File gvcf_summary_metrics = summary_metrics
        File gvcf_detail_metrics = detail_metrics

        File? output_bam = provided_output_bam
        File? output_bam_index = provided_output_bam_index

        File output_vcf = vcf_output
        File output_vcf_index = vcf_output_index

        File xytyping = XYtyping.output_ratio
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
        String LIBRARY_NAME = "Solexa-272222"
        String PLATFORM_UNIT = "HLYLLDRXX.2"
        String PLATFORM = "illumina"
        String SEQUENCING_CENTER = "NYU"
        String READ_GROUP_NAME = "HLYLLDRXX"
    }

    command <<<

       current_date_time=$(date "+%Y-%m-%dT%H:%M:%S%z")
       #java -Xmx8G -jar /picard/picard.jar FastqToSam \
       gatk --java-options "-Xmx8G" \
       FastqToSam \
       FASTQ=~{sample_fastq.r1} \
       FASTQ2=~{sample_fastq.r2} \
       OUTPUT=~{sample_fastq.sample_name}_unmapped.bam \
       READ_GROUP_NAME=~{READ_GROUP_NAME} \
       SAMPLE_NAME=~{sample_fastq.sample_name} \
       LIBRARY_NAME=~{LIBRARY_NAME} \
       PLATFORM_UNIT=~{PLATFORM_UNIT} \
       PLATFORM=~{PLATFORM} \
       SEQUENCING_CENTER=~{SEQUENCING_CENTER} \
       RUN_DATE=$current_date_time

    >>>

    output{
        File ubam = "~{sample_fastq.sample_name}_unmapped.bam"

    }
    runtime{
        #docker: "us.gcr.io/broad-gotc-prod/dragmap:1.1.2-1.2.1-2.26.10-1.11-1643839530"
        #cpu: 1
        docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
        cpu: 1
        memory: "9000 MiB"
        runtime_minutes: 240
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
task CorrectedVCF {
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
        File output_vcf = "~{output_name}"
        File output_vcf_index = "~{output_name}.tbi"
    }
    runtime{
        docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
        cpu: 1
        memory: "6000 MiB"
        runtime_minutes: 180
    }
}
