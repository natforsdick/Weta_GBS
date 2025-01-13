#!/bin/bash -e
#SBATCH -A ga03048
#SBATCH -J bowtie_index
#SBATCH -c 12
#SBATCH --mem=20G
#SBATCH --partition=large
#SBATCH --time=02:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err

############
# run_bowtie2_index
# Nat Forsdick, 2021-01-13
############

# This script is for indexing a reference genome using Bowtie2, prior to read alignment.

############
# MODULES
module purge
module load Bowtie2/2.3.5-GCC-7.4.0
############

############
# PARAMS
REFFILE=Deinacrida-fallai-asm-10kb
REFDIR=/nesi/project/ga03186/ref/
REF=${REFDIR}${REFFILE}
############

cd $REFDIR

if [ ! -e ${REFFILE}.fasta ]; then
gunzip ${REFFILE}.fasta.gz
fi

echo "Indexing $REFFILE"
bowtie2-build -f --t 8 --large-index ${REF}.fasta ${REF}
echo "completed indexing"

