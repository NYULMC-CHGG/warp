version 1.0

# imports

#workflow 
workflow ImputeFormat{
    input{
        File vcf
        File vcf_index
        File force_call = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/forceCall/sort_vv.txt"
        File ref_fasta = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/dragen_reference/Homo_sapiens_assembly38_masked.fasta"
        File ref_index = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/dragen_reference/Homo_sapiens_assembly38_masked.fasta.fai"
        File ref_dict = "/gpfs/data/deciphEHRlab/pipeline/reference/hg38/dragen_reference/Homo_sapiens_assembly38_masked.dict"
        Boolean isParabricks = false
        String sampleID
        String partition = 'deciphehr'

    }

    call Format {
        input:
            vcf = vcf,
            vcf_index = vcf_index,
            force_call = force_call,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            ref_dict = ref_dict,
            sampleID = sampleID,
            isParabricks = isParabricks,
            cpu = 1,
            mem_gb = 4,
            partition = partition,
            runtime_minutes = 180
    }

    output{
        File vcf_ouput = Format.formated_vcf
        File vcf_index_output = Format.formated_vcf_index
    }
}

task Format{
    input{
        File vcf
        File vcf_index
        File force_call
        File ref_fasta
        File ref_index
        File ref_dict
        Boolean isParabricks
        String sampleID
        Int max_indel_length = 50
        Int min_meanDP = 1
        String partition
        Int cpu
        Int runtime_minutes
        Int mem_gb

    }
    command <<<
        set -euxo pipefail

        module load bcftools
        module load tabix
        module load gatk/4.2.1.0
        module load vcftools

        # RM blocks. Not NON REF, Can be NON_REF or * depending on caller
        rvcf=~{sampleID}_called.vcf
        if [ ~{isParabricks} == "false" ]; then
            zcat ~{vcf} | awk '$5 == "<NON_REF>" { next } { print }' > $rvcf
        else
            zcat ~{vcf} | awk '$5 == "<*>" { next } { print }' > $rvcf
        fi
        #SPLIT ALLELES & RM Long INDELS
        echo "Splitting alleles..."
        svcf=~{sampleID}_called_split_gatk.vcf
        gatk LeftAlignAndTrimVariants -R ~{ref_fasta} -V $rvcf -O $svcf --split-multi-allelics --max-indel-length ~{max_indel_length}

        #RM NON-REF empth alleles. Not empty allele can be NON_REF or * depending on caller
        nrvcf=~{sampleID}_called_split_nr.vcf
        if [ ~{isParabricks} == "false" ]; then
            awk '$5 == "<NON_REF>" { next } { print }' $svcf > $nrvcf
        else
            awk '$5 == "<*>" { next } { print }' $svcf > $nrvcf
        fi

        #RM variants w/no GL for downstream analysis
        vcfm=~{sampleID}_called_split_rd1
        echo "remove zero reads..."
        vcftools --vcf $nrvcf --out $vcfm --min-meanDP ~{min_meanDP} --recode --recode-INFO-all

        #Add in missing variants in force call list as missing
        cv=~{sampleID}_test_called_vars.txt
        vcfmv=${vcfm}'.recode.vcf'
        ln=$(grep -n '#' $vcfmv | tail -1 |cut -d':' -f1)
        a=$(($ln+1))
        tail -n+$a $vcfmv | awk 'BEGIN{OFS=":"} {print $1,$2,$4,$5}' > $cv

        #Compare to 1000G
        cvs=~{sampleID}_sort_cv.txt
        echo "Comparing..."
        sort $cv > $cvs
        mfile=~{sampleID}_compare.txt
        comm -23 ~{force_call} $cvs > $mfile
        #Add back missing vars
        awk -F':' 'BEGIN {OFS="\t"}{print $1,$2,".",$3,$4,".",".",".","GT","./."}' $mfile >> $vcfmv

        #Re-sort and Zip
        echo "Re-sorting and zipping..."
        bcftools sort -O z -o ~{sampleID}_called_all_sort.vcf.gz $vcfmv
        #bgzip -c ~{sampleID}_called_all_sort.vcf.gz
        bcftools index -t ~{sampleID}_called_all_sort.vcf.gz

        echo "Finito!"

    >>>
    output {
        File formated_vcf = "~{sampleID}_called_all_sort.vcf.gz"
        File formated_vcf_index = "~{sampleID}_called_all_sort.vcf.gz.tbi"
    }
    runtime{
        cpu: cpu
        memory : mem_gb+" GiB"
        runtime_minutes: runtime_minutes
        partition : partition

    }
}