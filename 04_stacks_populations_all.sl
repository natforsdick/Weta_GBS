#!/bin/bash -e
#SBATCH -J stacks_pop
#SBATCH -A ga03186
#SBATCH --time=02:40:00
#SBATCH --mem=300M
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
INDIR1=/nesi/nobackup/ga03186/Weta_GBS_bwa_B1_B2/03_ref_map/
INDIR2=/nesi/nobackup/ga03186/Weta_GBS_bowtie_B1_B2/03_ref_map/
OUTDIR1=/nesi/nobackup/ga03186/Weta_GBS_bwa_B1_B2/04_populations/
OUTDIR2=/nesi/nobackup/ga03186/Weta_GBS_bowtie_B1_B2/04_populations/
#poplist="Mah Fallai Batch2_pop weta_all_popmap Mahoenui_all"
#poplist="Mahoenui_all"
#poplist="weta_all_popmap"
poplist="weta_all_popmap_batch"
#POPMAP=/nesi/project/ga03186/ref/Weta_GBS_Batch2_POP_blankrem.txt
REFDIR=/nesi/project/ga03186/ref/
############
which populations

for pop in $poplist;
do

if [ ! -e ${OUTDIR1}${pop}_a/ ]; then
mkdir -p ${OUTDIR1}${pop}_a/
mkdir -p ${OUTDIR1}${pop}_b/
fi

if [ ! -e ${OUTDIR2}${pop}_a/  ]; then
mkdir -p ${OUTDIR2}${pop}_a/
mkdir -p ${OUTDIR2}${pop}_b/
fi


#echo "Running stacks populations for ${pop}, no missing data"
#populations -P ${INDIR1}${pop} -O ${OUTDIR1}${pop}_a/ -M ${REFDIR}${pop}.txt -t 8 --min-maf 0.05 --hwe --fstats --smooth-popstats --smooth --bootstrap --vcf --structure --genepop -r 1 --write-single-snp

#echo "Running stacks populations for ${pop}, 30% missing data"
#populations -P ${INDIR1}${pop} -O ${OUTDIR1}${pop}_b/ -M ${REFDIR}${pop}.txt -t 8 --min-maf 0.05 --hwe --fstats --smooth-popstats --smooth --bootstrap --vcf --structure --genepop -r 0.7 --write-single-snp
#echo "Completed stacks processing for ${pop} bwa"

#echo "Running stacks populations for ${pop} bowtie, no missing data"
#populations -P ${INDIR2}${pop} -O ${OUTDIR2}${pop}_a/ -M ${REFDIR}${pop}.txt -t 8 --min-maf 0.05 --hwe --fstats --smooth-popstats --smooth --bootstrap --vcf --structure --genepop -r 1 --write-single-snp
   

echo "Running stacks populations for ${pop}, 30% missing data"
populations -P ${INDIR2}${pop} -O ${OUTDIR2}${pop}_b/ -M ${REFDIR}${pop}.txt -t 8 --min-maf 0.05 --hwe --fstats --smooth-popstats --smooth --bootstrap --vcf --structure --genepop -r 0.7 --write-single-snp
echo "Completed stacks processing for ${pop} bowtie"

done
