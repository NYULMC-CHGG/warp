version 1.0


# imports
import "../structs/dna_seq/DNASeqStructs.wdl"

workflow DeciphehrBatch{

    #inputs
    input{
        Array[SampleFastq] samples
    }


  scatter (sample in samples)
  {
    call Atask 
    {
        input:
            sample_fastq = sample
    }
  }

  output{
        Array[File] outfile1 = Atask.fjson
        Array[File] outfile2 = Atask.num2
        
  }

}

task Atask{
    input{
        SampleFastq sample_fastq
    }
    command <<<
        JSON_STRING='{"sample_name":"~{sample_fastq.sample_name}","base_file_name":"~{sample_fastq.base_file_name}","final_gvcf_base_name":"~{sample_fastq.final_gvcf_base_name}","unmapped_bam_suffix":"~{sample_fastq.unmapped_bam_suffix}"}'
        echo $JSON_STRING > ~{sample_fastq.sample_name}.json
        echo "HELLO" > ~{sample_fastq.sample_name}_2.json

    >>>
    output{
        File fjson = "~{sample_fastq.sample_name}.json"
        File num2 = "~{sample_fastq.sample_name}_2.json"
    }
    runtime{
        cpu: 1
        memory: "1000 MiB"
        runtime_minutes: 2
    }
}


