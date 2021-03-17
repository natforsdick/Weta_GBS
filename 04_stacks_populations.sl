#!/bin/bash -e
#SBATCH -J stacks_pop
#SBATCH -A ga03186
#SBATCH --time=01:45:00
#SBATCH --mem=100M
#SBATCH --cpus-per-task=4
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz

############
# 04_stacks_populations.sl
# Nat Forsdick, 2021-01-25
# For running the populations tool in stacks to call and filter SNPs.
############

############
# MODULES
module purge 
module load Stacks/2.41-gimkl-2018b
############

############
# PARAMS
INDIR=/nesi/nobackup/ga03186/Weta_GBS_Batchcombo_adap/04_ref_map/
OUTDIR=/nesi/nobackup/ga03186/Weta_GBS_Batchcombo_adap/05_populations/
poplist="Het Mah Fallai Mahoenui_all Weta_GBS_Batch2_POP_MI Weta_GBS_Batch2_POP_SR" 
#POPMAP=/nesi/project/ga03186/ref/Weta_GBS_Batch2_POP_blankrem.txt
REFDIR=/nesi/project/ga03186/ref/
############
which populations

for pop in $poplist;
do
if [ ! -e ${OUTDIR}${pop}_a/ ]; then
mkdir -p ${OUTDIR}${pop}_a/
mkdir -p ${OUTDIR}${pop}_b/
fi

echo "Running stacks populations for ${pop}, no missing data"
populations -P ${INDIR}Weta_GBS_adap_${pop} -O ${OUTDIR}${pop}_a/ -M ${REFDIR}${pop}.txt -t 8 --min-maf 0.05 --hwe --fstats --smooth-popstats --smooth --bootstrap --vcf --structure --genepop -r 1 --write-single-snp

echo "Running stacks populations for ${pop}, 30% missing data"
populations -P ${INDIR}Weta_GBS_adap_${pop} -O ${OUTDIR}${pop}_b/ -M ${REFDIR}${pop}.txt -t 8 --min-maf 0.05 --hwe --fstats --smooth-popstats --smooth --bootstrap --vcf --structure --genepop -r 0.7 --write-single-snp
echo "Completed stacks processing for ${pop}"
done
