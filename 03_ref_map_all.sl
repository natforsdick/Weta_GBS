#!/bin/bash -e
#SBATCH -J ref_map_bwa
#SBATCH -A ga03186
#SBATCH --time=02:00:00
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
INDIR1=/nesi/nobackup/ga03186/Weta_GBS_bwa_B1_B2/02_bwa_sorted/ 
OUTDIR1=/nesi/nobackup/ga03186/Weta_GBS_bwa_B1_B2/03_ref_map/
REF=/nesi/nobackup/ga03186/reference/
list=/nesi/project/ga03186/ref/
#refstack=Weta_GBS_Batch2_POP.txt
#poplist="Het Mah Fallai 
poplist1="weta_all_popmap Het Mahoenui_all Fallai"

INDIR2=/nesi/nobackup/ga03186/Weta_GBS_bowtie_B1_B2/02_bowtie_sorted/
OUTDIR2=//nesi/nobackup/ga03186/Weta_GBS_bowtie_B1_B2/03_ref_map/
############
for pop in $poplist1
do
if [ ! -e ${OUTDIR1}${pop} ]; then
mkdir -p ${OUTDIR1}${pop}
mkdir -p ${OUTDIR2}${pop}
fi

cd ${OUTDIR1}${pop}

echo "Running ref_map for ${pop}, ${OUTDIR1}"
srun ref_map.pl -T 10 --samples $INDIR1 --popmap ${list}${pop}.txt -o ${OUTDIR1}${pop}
echo "Finished ref_map for ${pop}, ${OUTDIR1}"

cd ${OUTDIR2}${pop}
 
echo "Running ref_map for ${pop}, ${OUTDIR2}"
srun ref_map.pl -T 10 --samples $INDIR2 --popmap ${list}${pop}.txt -o ${OUTDIR2}${pop}
echo "Finished ref_map for ${pop}, ${OUTDIR2}"

done



