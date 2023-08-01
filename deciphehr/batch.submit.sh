#!/bin/bash
####################################
# Submit jobs using csv file of samples
# sampleID,/path/to/R1.fastq,/path/to/R2.fastq
# requirements: 
# jq (https://jqlang.github.io/jq/download/)
# run: 
# $ batch.submit.sh sample.manifest
# admin: jonathan.mccafferty@nyulangone.org
####################################

BATCH=$1
# check for input directory of json files
if [ ! -d 'inputs' ]
then
    mkdir inputs
fi
if [ ! -d 'jobLog' ]
then
    mkdir jobLog
fi

DATE=$(date +%Y%m%d)
echo "sample slurmID runID date R1 R2" | tr ' ' , >> jobLog/${DATE}_run.log
# read in samples and submit jobs to slurm
while IFS=',' read s r1 r2
do
    # this will generate a random number and attach it to the sample id for tracking and submit jobs to slurm
    # all .json input files will go in to inputs/ directory
    # job run log will output to jobLog/ directory
    RID='rid'$RANDOM
    SAM=${s}_${RID}${DATE}
    F1=${r1}
    F2=${r2}
    BAM=".bam"
    JSON='{"Deciphehr.sample_fastq":{"sample_name":"'${SAM}'","r1":"'${F1}'","r2":"'${F2}'","base_file_name":"'${SAM}'","final_gvcf_base_name":"'${SAM}'","unmapped_bam_suffix":"'${BAM}'"}}'
    echo ${JSON} | jq . > inputs/${SAM}.json
    #cat inputs/${SAM}.json | jq .
    jobid=$(sbatch --parsable deciphehr.batch.sh inputs/${SAM}.json)
    echo ${s} ${jobid} ${RID} ${DATE} ${F1} ${F2} | tr ' ' , >> jobLog/${DATE}_run.log
    
    
done < $BATCH
