#!/bin/bash
#SBATCH --job-name=dcipher
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=96:00:00
#SBATCH --partition=cpu_medium
#SBATCH -o ./jobLog/JOB%j.out
#SBATCH -e ./jobLog/JOB%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jonathan.mccafferty@nyulangone.org

# Load any necessary modules
module load jdk

JAR=/gpfs/scratch/mccafj02/DECIPHEHR/pipeline/cromwell/cromwell-85.jar
CONFIG=/gpfs/scratch/mccafj02/DECIPHEHR/pipeline/cromwell/cromwell.config
WDL=/gpfs/scratch/mccafj02/DECIPHEHR/pipeline/warp/deciphehr/Deciphehr.wdl
INPUT=$1
OPTIONS=/gpfs/scratch/mccafj02/DECIPHEHR/pipeline/driver/wdl_options.json
CE=$PWD/cromwell-executions

java -jar -Xms2g -Xmx4g -Dconfig.file=${CONFIG} ${JAR} run ${WDL} -i ${INPUT} -o ${OPTIONS}
