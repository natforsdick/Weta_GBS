#!/bin/bash -e
#SBATCH -J stacks_pop
#SBATCH --time=00:10:00
#SBATCH --mem=200M
#SBATCH --cpus-per-task=2
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err

############
# 04_stacks_populations.sl
# Nat Forsdick, 2021-01-25
# For running the populations tool in stacks to call and filter SNPs.
############

############
# MODULES
module purge
module load Stacks/2.65-GCC-11.3.0
############

############
# PARAMS
INDIR=/path/to/03_ref_map_PE/Weta_GBS_Batch2_all/
OUTDIR=/path/to/03_ref_map_PE/all2/

POPMAP=/path/to/ref/Weta_GBS_Batch2_all.txt
############

which populations

mkdir -p ${OUTDIR}

echo "Running stacks"
populations -P ${INDIR} -O ${OUTDIR} -M ${POPMAP} -t 8 --vcf --write-random-snp --ordered-export

# We can then use vcftools to check that the output only contains biallelic SNPs and no indels
ml purge; ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
vcftools --vcf populations.snps.vcf --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out weta-bialsnps-check
