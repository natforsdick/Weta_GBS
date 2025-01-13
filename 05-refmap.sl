#!/bin/bash -e
#SBATCH -J ref_map
#SBATCH --time=01:30:00
#SBATCH --mem=3G
#SBATCH --cpus-per-task=12
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err

###########
# 03_ref_map.sl
# Nat Forsdick, 2021-01-21
# This script is to run ref_map for the GBS data.
# In this case, these are weta D. mahoenui GBS reads mapping to \
# the draft D. fallai reference genome.
###########

###########
# MODULES
module purge
module load Stacks/2.65-GCC-11.3.0
module list
############

############
# PARAMS
INDIR=/path/to/03_bowtie/BAM/
OUTDIR=/path/to/03_ref_map_PE/
list=/path/to/ref/
poplist="Weta_GBS_Batch2_all"
############
for pop in $poplist
do
    if [ ! -e ${OUTDIR}${pop} ]; then
        mkdir -p ${OUTDIR}${pop}
    fi

    cd ${OUTDIR}${pop}

    echo "Running ref_map for ${pop}"
    srun ref_map.pl -T 24 --samples $INDIR --popmap ${list}${pop}.txt -o ${OUTDIR}${pop}
    echo "Finished ref_map for ${pop}"
done
