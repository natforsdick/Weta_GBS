#!/bin/bash -e
#SBATCH -J vcf2adegenet
#SBATCH -A ga03186
#SBATCH --time=00:03:00
#SBATCH --mem=300M
#SBATCH --cpus-per-task=4
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz


############
# 05_vcf2adegenet.sl
# Nat Forsdick, 2021-03-31
# For converting vcfs to RAW and MAP formats for Adegenet
############

############
# MODULES
module purge
module load Stacks/2.41-gimkl-2018b
############

############
# PARAMS
INDIR1=/nesi/nobackup/ga03186/Weta_GBS_bwa_B1_B2/04_populations/
INDIR2=/nesi/nobackup/ga03186/Weta_GBS_bowtie_B1_B2/04_populations/
#poplist="Fallai Batch2_pop weta_all_popmap Mahoenui_all weta_all_popmap_batch Mahoenui_all_pop" #Mah
poplist="Mahoenui_all_pop" #
#poplist="weta_all_popmap_batch"
REFDIR=/nesi/project/ga03186/ref/
OUTDIR1=/nesi/nobackup/ga03186/Weta_GBS_bwa_B1_B2/05_plink/
OUTDIR2=/nesi/nobackup/ga03186/Weta_GBS_bowtie_B1_B2/05_plink/
############

for pop in $poplist;
do

#mkdir ${OUTDIR1}${pop}_a
#mkdir ${OUTDIR1}${pop}_b
#mkdir ${OUTDIR2}${pop}_a
#mkdir ${OUTDIR2}${pop}_b

    populations -V ${INDIR1}${pop}_a/${pop}_a.vcf -O ${OUTDIR1}${pop}_a/ -M ${REFDIR}${pop}.txt -t 8 --plink
    populations -V ${INDIR1}${pop}_b/${pop}_b.vcf -O ${OUTDIR1}${pop}_b/ -M ${REFDIR}${pop}.txt -t 8 --plink

    populations -V ${INDIR2}${pop}_a/${pop}_a.vcf -O ${OUTDIR2}${pop}_a/ -M ${REFDIR}${pop}.txt -t 8 --plink
    populations -V ${INDIR2}${pop}_b/${pop}_b.vcf -O ${OUTDIR2}${pop}_b/ -M ${REFDIR}${pop}.txt -t 8 --plink

    done
