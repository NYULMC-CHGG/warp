#!/bin/bash
#SBATCH --job-name=dcipher
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --partition=cpu_medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@nyulangone.org

# Load any necessary modules
module load jdk

JAR=/gpfs/data/chaklab/data/DECIPHEHR/pipeline/cromwell/cromwell-85.jar
CONFIG=/gpfs/data/chaklab/data/DECIPHEHR/pipeline/cromwell/cromwell.config
WDL=/gpfs/data/chaklab/data/DECIPHEHR/pipeline/warp/deciphehr/Deciphehr.wdl
INPUT=/gpfs/data/chaklab/data/DECIPHEHR/pipeline/warp/deciphehr/ex_inputs.json
OPTIONS=/gpfs/data/chaklab/data/DECIPHEHR/pipeline/warp/deciphehr/wdl_options.json
CE=$PWD/cromwell-executions

java -jar -Xms10g -Xmx18g -Dconfig.file=${CONFIG} ${JAR} run ${WDL} -i ${INPUT} -o ${OPTIONS}
