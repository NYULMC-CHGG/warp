version 1.0

# imports
import "../structs/dna_seq/DNASeqStructs.wdl"
import "Deciphehr.wdl" as RunWorkFlow

workflow DragenGATK {
    input{
        Array[SampleFastq] samples
    }

    scatter(sample in samples)
    {
        call RunWorkFlow.Deciphehr as dec
        {
           input:
               sample_fastq = sample
        }
        


    }
    
    output{
        Array[File?] cross_check_fingerprints_metrics = dec.cross_check_fingerprints_metrics
        Array[File] duplicate_metrics = dec.duplicate_metrics
        Array[File?] output_bqsr_reports = dec.output_bqsr_reports

        Array[File] gvcf_summary_metrics = dec.gvcf_summary_metrics
        Array[File] gvcf_detail_metrics = dec.gvcf_detail_metrics

        Array[File?] output_bam = dec.output_bam
        Array[File?] output_bam_index = dec.output_bam_index

        Array[File] output_vcf = dec.output_vcf
        Array[File] output_vcf_index = dec.output_vcf_index

        Array[File] xytyping = dec.xytyping
        Array[File?] xyvcf = dec.xyvcf
        Array[File?] xyvcf_index = dec.xyvcf_index

    }


}

