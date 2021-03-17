#!/bin/bash -e
#SBATCH -J ref_map
#SBATCH -A ga03186
#SBATCH --time=01:00:00
#SBATCH --mem=3G
#SBATCH --cpus-per-task=6
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz

###########
# 03_ref_map.sl
# Nat Forsdick, 2021-01-21
# This script is to run ref_map for the GBS data.
# In this case, these are weta D. mahoenui GBS reads mapping to \
# the draft D. fallai reference genome.
# This is based on the workflow of Victoria Twort. 
###########

###########
# MODULES
module purge
module load Stacks/2.41-gimkl-2018b
module list
############

############
# PARAMS
INDIR=/nesi/nobackup/ga03186/Weta_GBS_Batchcombo_adap/03_stacks_sorted/ 
OUTDIR=/nesi/nobackup/ga03186/Weta_GBS_Batchcombo_adap/04_ref_map/
REF=/nesi/nobackup/ga03186/reference/
list=/nesi/project/ga03186/ref/
#refstack=Weta_GBS_Batch2_POP.txt
#poplist="Het Mah Fallai 
poplist="weta_all_popmap"
############
for pop in $poplist
do
if [ ! -e ${OUTDIR}Weta_GBS_adap_${pop} ]; then
mkdir -p ${OUTDIR}Weta_GBS_adap_${pop}
fi

cd ${OUTDIR}Weta_GBS_adap_${pop}

echo "Running ref_map for ${pop}"
srun ref_map.pl -T 10 --samples $INDIR --popmap ${list}${pop}.txt -o ${OUTDIR}Weta_GBS_adap_${pop}
echo "Finished ref_map for ${pop}"
done
