#!/bin/bash -e
#SBATCH -J B1_stacks_demux
#SBATCH -A ga03186
#SBATCH --time=05:00:00
#SBATCH --mem=500M
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
mkdir -p /nesi/nobackup/ga03186/Weta_GBS_Batch1/01_stacks_demux_adap
# adapter sequences based on those in the Stacks manual for Illumina HiSeq data.
srun process_radtags -p /nesi/project/ga03186/data/Weta_GBS_Batch1/raw/ -i gzfastq -b /nesi/project/ga03186/ref/B1_barcodes_stacks.txt -o /nesi/nobackup/ga03186/Weta_GBS_Batch1/01_stacks_demux_adap/ -e pstI -c -q --inline_null --adapter_1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --adapter_2 GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 2 
#Originally run with:
#--adapter_1 AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTT 

#srun process_radtags -p /nesi/project/ga03186/data/Weta_GBS_Batch1/raw/ -f H7YLCBCXX-NZGL01305_L001_R1.fastq.gz -i gzfastq /nesi/project/ga03186/data/Weta_GBS_Batch1/ref/barcodes_stacks.txt -o /nesi/nobackup/ga03186/Weta_GBS_Batch1/01_stacks_demux -e pstI -c -q --inline_null --adapter_1 AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTT --adapter_mm 3


