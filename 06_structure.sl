#!/bin/bash -e
#SBATCH -A uoo02327
#SBATCH -J faststr_bial
#SBATCH --time 00:05:00 #
#SBATCH -c 12
#SBATCH --mem=2G
#SBATCH --partition=large
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --mail-type=FAIL,END
#SBATCH --output structure.%j.out # CHANGE number for new run
#SBATCH --error structure.%j.err #  CHANGE number for new run

############
# 06_structure.sl
# Nat Forsdick, 2021-04-06
# Script to run Weta GBS data in structure.
################

###########
# MODULES #
module load Structure/2.3.4-GCCcore-7.4.0

############
# PARAMS   #
poplist="Fallai Het Mahoenui_all"

INDIR=/nesi/nobackup/ga03186/Weta_GBS_bowtie_B1_B2/04_populations/
OUTDIR=/nesi/nobackup/ga03186/Weta_GBS_bowtie_B1_B2/05_structure/

Klist="1 2 3 4 5 6 7 8 9 10"
############

for K in $Klist
    do

        structure -K $K -o $OUTDIR/${K} -m ./mainparams_90.txt -e ./extraparams_all.txt

        done
