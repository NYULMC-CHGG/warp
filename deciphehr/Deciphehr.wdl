version 1.0

# imports
import "../structs/dna_seq/DNASeqStructs.wdl"
import "../tasks/broad/GermlineVariantDiscovery.wdl" as Calling
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
            allow_empty_ref_alt = true

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
                    mem_size_mb = 20000,
                    id = x.id

                }
            }
            #gather vcf outputs and merge into a single gvcf file
            Array[File] vcfs_to_merge = flatten(select_all([HaplotypeCallerXY.vcf]))
            Array[File] vcf_indices_to_merge = flatten(select_all([HaplotypeCallerXY.vcf_index]))
            call Calling.MergeVCFs as MergeVCFs {
                input:
                    input_vcfs = vcfs_to_merge,
                    input_vcfs_indexes = vcf_indices_to_merge,
                    output_vcf_name = sample_fastq.sample_name+"_XY_autosomes_merged.vcf.gz",
            }
            call FinalVCF {
            input:
                fasta_ref = references.reference_fasta.ref_fasta,
                fasta_index = references.reference_fasta.ref_fasta_index,
                fasta_dict = references.reference_fasta.ref_dict,
                xy_vcfs = MergeVCFs.output_vcf,
                input_vcfs_indexes = MergeVCFs.output_vcf_index,
                autosome_vcf = wgs.output_vcf,
                autosome_vcf_index = wgs.output_vcf_index,
                output_name = sample_fastq.sample_name+"_Final_merged.g.vcf.gz"
            }
        
        }

    

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

        File xytyping = XYtyping.output_ratio
        File? finaVCF = FinalVCF.vcf
        File? finalVCF_index = FinalVCF.vcf_index
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

        echo "sample x/y sex" | tr ' ' , > ~{samid}.XYratio.csv
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
        cpu: 2
        memory: "5000 MiB"
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
        -ploidy ~{ploidy} \
        -ERC GVCF \
        -G StandardAnnotation \
        -G StandardHCAnnotation \
        -G AS_StandardAnnotation \
        -RF OverclippedReadFilter
    >>>
    output{
        File vcf = "~{id}.vcf.gz"
        File vcf_index = "~{id}.vcf.gz.tbi"
    }
    runtime{
        docker:"us.gcr.io/broad-gatk/gatk:4.3.0.0"
        cpu: 5
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
        cpu: 5
        memory: "10000 MiB"
    }
}