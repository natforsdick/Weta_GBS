#!/bin/bash -e
#SBATCH -J stacks_demux
#SBATCH --time=10:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=2
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH --mail-type=FAIL,END

###########
# 01-stacks-process-radtags.sl
# Nat Forsdick, 2021-01-21
###########

# This script is to demultiplex, clean, and filter paired-end GBS sequence reads with combinatorial barcodes.
# Based on the script used by Victoria Twort for single-end GBS sequencing data for Weta Batch 1.

###########
#  MODULES
module purge
module load Stacks/2.41-gimkl-2018b
###########

###########
# PARAMS
INDIR=/path/to/data/weta-GBS/
OUTDIR=/path/to/01_stacks_demux_PE2/
samp="HFJKLCCX2_6_201027_FD09254671_Other__R_200824_ROBELS1_LIBX10_M003" # sequence fastq file
###########

if [ ! -e ${OUTDIR}  ]; then
    mkdir -p ${OUTDIR}
    fi

# process_radtags -1 pair_1 -2 pair_2 [-b barcode_file] -o out_dir -e enzyme [-c] [-q] [-r] [-t len]
#  -P = paired-end data, -p = in_dir, -i = input filetype, -b = barcode file, -c = clean data, \
# -q = quality filter
srun process_radtags -P -1 ${INDIR}${samp}_R1.fastq.gz \
    -2 ${INDIR}${samp}_R2.fastq.gz \
    -i gzfastq -b ${INDIR}barcodes_stacks_batch2.csv -o ${OUTDIR} \
    -e pstI -c -q --inline_inline \
    --adapter_1 AGATCGGAAGAGC --adapter_mm 3
