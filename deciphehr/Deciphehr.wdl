version 1.0

# imports
import "../structs/dna_seq/DNASeqStructs.wdl"
import "../pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl" as ToWGSgermline

# workflow
workflow Deciphehr{

    #inputs
    input{
        File dnaSeqSSRef = "/gpfs/data/chaklab/home/mccafj02/home/repos/warp/structs/dna_seq/DNASeqSingleSampleReferences.json"
        File dragMapRef = "/gpfs/data/chaklab/home/mccafj02/home/repos/warp/structs/dna_seq/DragmapReference.json"
        File papiSettings = "/gpfs/data/chaklab/home/mccafj02/home/repos/warp/structs/dna_seq/PapiSettings.json"
        File variantCallScatter = "/gpfs/data/chaklab/home/mccafj02/home/repos/warp/structs/dna_seq/VariantCallingScatterSettings.json"
        SampleFastq sample_fastq
        File wgs_coverage_interval_list = "/gpfs/data/chaklab/home/mccafj02/pipelines/reference/hg38/v0/wgs_coverage_regions.hg38.interval_list"

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
            allow_empty_ref_alt = true

    }

    
    ## call sex typing
    
    ## if Male
        ## call variants on X and Y par/non-par regions

        ## merge vcfs

        ## filter

    ## merge with dragen vcf

    ## output everything
    output{
        Array[File] quality_yield_metrics = wgs.quality_yield_metrics

        Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = wgs.unsorted_read_group_base_distribution_by_cycle_pdf
        Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = wgs.unsorted_read_group_base_distribution_by_cycle_metrics
        Array[File] unsorted_read_group_insert_size_histogram_pdf = wgs.unsorted_read_group_insert_size_histogram_pdf
        Array[File] unsorted_read_group_insert_size_metrics = wgs.unsorted_read_group_insert_size_metrics
        Array[File] unsorted_read_group_quality_by_cycle_pdf = wgs.unsorted_read_group_quality_by_cycle_pdf
        Array[File] unsorted_read_group_quality_by_cycle_metrics = wgs.unsorted_read_group_quality_by_cycle_metrics
        Array[File] unsorted_read_group_quality_distribution_pdf = wgs.unsorted_read_group_quality_distribution_pdf
        Array[File] unsorted_read_group_quality_distribution_metrics = wgs.unsorted_read_group_quality_distribution_metrics

        File read_group_alignment_summary_metrics = wgs.read_group_alignment_summary_metrics
        File read_group_gc_bias_detail_metrics = wgs.read_group_gc_bias_detail_metrics
        File read_group_gc_bias_pdf = wgs.read_group_gc_bias_pdf
        File read_group_gc_bias_summary_metrics = wgs.read_group_gc_bias_summary_metrics

        File? cross_check_fingerprints_metrics = wgs.cross_check_fingerprints_metrics

        File selfSM = wgs.selfSM
        Float contamination = wgs.contamination

        File calculate_read_group_checksum_md5 = wgs.calculate_read_group_checksum_md5

        File agg_alignment_summary_metrics = wgs.agg_alignment_summary_metrics
        File agg_bait_bias_detail_metrics = wgs.agg_bait_bias_detail_metrics
        File agg_bait_bias_summary_metrics = wgs.agg_bait_bias_summary_metrics
        File agg_gc_bias_detail_metrics = wgs.agg_gc_bias_detail_metrics
        File agg_gc_bias_pdf = wgs.agg_gc_bias_pdf
        File agg_gc_bias_summary_metrics = wgs.agg_gc_bias_summary_metrics
        File agg_insert_size_histogram_pdf = wgs.agg_insert_size_histogram_pdf
        File agg_insert_size_metrics = wgs.agg_insert_size_metrics
        File agg_pre_adapter_detail_metrics = wgs.agg_pre_adapter_detail_metrics
        File agg_pre_adapter_summary_metrics = wgs.agg_pre_adapter_summary_metrics
        File agg_quality_distribution_pdf = wgs.agg_quality_distribution_pdf
        File agg_quality_distribution_metrics = wgs.agg_quality_distribution_metrics
        File agg_error_summary_metrics = wgs.agg_error_summary_metrics

        File? fingerprint_summary_metrics = wgs.fingerprint_summary_metrics
        File? fingerprint_detail_metrics = wgs.fingerprint_detail_metrics

        File wgs_metrics = wgs.wgs_metrics
        File raw_wgs_metrics = wgs.raw_wgs_metrics

        File duplicate_metrics = wgs.duplicate_metrics
        File? output_bqsr_reports = wgs.output_bqsr_reports

        File gvcf_summary_metrics = wgs.gvcf_summary_metrics
        File gvcf_detail_metrics = wgs.gvcf_detail_metrics

        File? output_bam = wgs.output_bam
        File? output_bam_index = wgs.output_bam_index

        File output_cram = wgs.output_cram
        File output_cram_index = wgs.output_cram_index
        File output_cram_md5 = wgs.output_cram_md5

        File validate_cram_file_report = wgs.validate_cram_file_report

        File output_vcf = wgs.output_vcf
        File output_vcf_index = wgs.output_vcf_index
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
        cpu: 2
        memory: "2000 MiB"
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
        cpu: 2
        memory: "10000 MiB"
    }
}