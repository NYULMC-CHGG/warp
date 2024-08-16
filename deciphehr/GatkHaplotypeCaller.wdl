version 1.0

#imports
import "../pipelines/broad/dna_seq/germline/variant_calling/VariantCalling.wdl" as VarCall


workflow GatkHaplotypeCaller{
    
    input{

        Boolean run_dragen_mode_variant_calling = true
        Boolean use_spanning_event_genotyping = false
        File calling_interval_list = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/v0/wgs_calling_regions.hg38.interval_list"
        File evaluation_interval_list = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/v0/wgs_evaluation_regions.hg38.interval_list"
        File ref_fasta = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/dragen_reference/Homo_sapiens_assembly38_masked.fasta"
	    File ref_fasta_index = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/dragen_reference/Homo_sapiens_assembly38_masked.fasta.fai"
	    File ref_dict = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/dragen_reference/Homo_sapiens_assembly38_masked.dict"
	    File ref_str = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/hg38.str.zip"
        File dbsnp_vcf = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
        File dbsnp_vcf_index =  "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
        
        String sample_name
        Int contamination = 0
        Int haplotype_scatter_count = 10
        Int break_bands_at_multiples_of = 100000
        File input_bam
        File input_bam_index
        Int agg_preemptible_tries = 0
        Boolean make_gvcf = true
        Boolean make_bamout = false
        Boolean use_gatk3_haplotype_caller = false
        Boolean use_dragen_hard_filtering = true
        File allele_vcf = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/forceCall/variant_1000g_key.vcf.gz"
        File allele_index = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/forceCall/variant_1000g_key.vcf.gz.tbi"
        Int cpu = 10
        String partition = "deciphehr"

    }

    call VarCall.VariantCalling as vc {
        input:
            run_dragen_mode_variant_calling = run_dragen_mode_variant_calling,
            use_spanning_event_genotyping = use_spanning_event_genotyping,
            calling_interval_list = calling_interval_list,
            evaluation_interval_list = evaluation_interval_list,
            haplotype_scatter_count = haplotype_scatter_count,
            break_bands_at_multiples_of = break_bands_at_multiples_of,
            contamination = contamination,
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            ref_str = ref_str,
            dbsnp_vcf = dbsnp_vcf,
            dbsnp_vcf_index = dbsnp_vcf_index,
            base_file_name = sample_name,
            final_vcf_base_name = sample_name,
            agg_preemptible_tries = agg_preemptible_tries,
            make_gvcf = make_gvcf,
            make_bamout = make_bamout,
            use_gatk3_haplotype_caller = use_gatk3_haplotype_caller,
            skip_reblocking = true,
            use_dragen_hard_filtering = use_dragen_hard_filtering,
            allele_vcf = allele_vcf,
            allele_index = allele_index
    }

    output{
        File summary = vc.vcf_summary_metrics
        File detail = vc.vcf_detail_metrics
        File vcf = vc.output_vcf
        File vcf_index = vc.output_vcf_index
    }

    






}

