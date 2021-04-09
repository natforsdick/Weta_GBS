#!/bin/bash -e
#SBATCH -J ref_map
#SBATCH -A ga03186
#SBATCH --time=00:30:00
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
#INDIR=/nesi/nobackup/ga03186/Weta_GBS_Batchcombo_adap/03_stacks_sorted/ 
INDIR=/nesi/nobackup/ga03186/Weta_GBS_Batch2/02_bowtie_PE_B2/BAM/
OUTDIR=/nesi/nobackup/ga03186/Weta_GBS_Batch2/03_ref_map_PE/
REF=/nesi/nobackup/ga03186/reference/
list=/nesi/project/ga03186/ref/
#poplist="Het Mah Fallai Mahoenui_all Weta_GBS_Batch2_POP_MI Weta_GBS_Batch2_POP_SR"
#poplist="Weta_GBS_Batch2_POP Weta_GBS_Batch2_POP_MI Weta_GBS_Batch2_POP_SR"
poplist="Weta_GBS_Batch2_POP_noU"
############
for pop in $poplist
do
if [ ! -e ${OUTDIR}${pop} ]; then
mkdir -p ${OUTDIR}${pop}
fi

cd ${OUTDIR}${pop}

echo "Running ref_map for ${pop}"
srun ref_map.pl -T 10 --samples $INDIR --popmap ${list}${pop}.txt -o ${OUTDIR}${pop}
echo "Finished ref_map for ${pop}"
done
