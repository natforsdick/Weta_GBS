#!/bin/bash -e
#SBATCH -J faststr
#SBATCH --time 02:00:00 #
#SBATCH -c 2
#SBATCH --mem=2G
#SBATCH --output faststr.%j.%a.out # CHANGE number for new run
#SBATCH --error faststr.%j.%a.err #  CHANGE number for new run
#SBATCH --array=0-9

############
# 08-faststructure.sl
# Nat Forsdick, 2021-04-06
# Script to run Weta GBS data in faststructure.
################

###########
# MODULES #
module purge; module load fastStructure/1.0-gimkl-2020a-Python-2.7.18

############
# PARAMS   #
INDIR=/path/to/03_ref_map_PE/all2/
INBED=weta.snpmiss60.mac3.thin.bed
filename=$(basename $INBED .bed)

echo "Beginning faststructure run for $k at "$(date)
mkdir -p /path/to/04-faststructure/faststr_${SLURM_ARRAY_TASK_ID}
cd /path/to/04-faststructure/faststr_${SLURM_ARRAY_TASK_ID}

for K in {1..5};
do
    structure.py -K $K --input=${INDIR}${filename} \
        --output=${filename}_${K} \
        --cv=100 --full --prior=simple
done
