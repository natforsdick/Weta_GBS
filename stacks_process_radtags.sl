#!/bin/bash -e
#SBATCH -J stacks_demux
#SBATCH -A ga03186
#SBATCH --time=10:00:00
#SBATCH --mem=800M
#SBATCH --cpus-per-task=2
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz

###########
# stacks_process_radtags.sl
# Nat Forsdick, 2021-01-21
###########

# This script is to demultiplex, clean, and filter paired-end GBS sequence reads with combinatorial barcodes.
# Based on the script used by Victoria Twort for single-end GBS sequencing data for Weta Batch 1.
# This is intended for processing Weta GBS Batch 2. 

###########
#  MODULES
module purge
module load Stacks/2.41-gimkl-2018b
###########
mkdir /nesi/nobackup/ga03186/Weta_GBS_Batch2/01_stacks_demux_PE

# process_radtags -1 pair_1 -2 pair_2 [-b barcode_file] -o out_dir -e enz [-c] [-q] [-r] [-t len]
#  -P = paired-end data, -p = in_dir, -i = input filetype, -b = barcode file, -c = clean data, \
# -q = quality filter, --inline_null: barcode is inline with sequence, occurs only on single-end read
srun process_radtags -P -1 /nesi/project/ga03186/data/Weta_GBS_Batch2/raw/HFJKLCCX2_6_201027_FD09254671_Other__R_200824_ROBELS1_LIBX10_M003_R1.fastq.gz -2 /nesi/project/ga03186/data/Weta_GBS_Batch2/raw/HFJKLCCX2_6_201027_FD09254671_Other__R_200824_ROBELS1_LIBX10_M003_R2.fastq.gz -i gzfastq -b /nesi/project/ga03186/data/Weta_GBS_Batch2/ref/barcodes_stacks_batch2.txt -o /nesi/nobackup/ga03186/Weta_GBS_Batch2/01_stacks_demux_PE/ -e
pstI -c -q --inline_inline --adapter_1 AGATCGGAAGAGC --adapter_2 AGATCGGAAGAGC --adapter_mm 3 

#srun process_radtags -p /nesi/project/ga03186/data/Weta_GBS_Batch2/raw/ -i gzfastq -f HFJKLCCX2_6_201027_FD09254671_Other__R_200824_ROBELS1_LIBX10_M003_R1.fastq.gz -b /nesi/project/ga03186/data/Weta_GBS_Batch2/ref/barcodes_stacks_batch2.txt -o /nesi/nobackup/ga03186/Weta_GBS_Batch2/01_stacks_demux/ -e pstI -c -q --inline_inline -y gzfastq
