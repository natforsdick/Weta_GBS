#!/bin/bash -e
#SBATCH --account=ga03186
#SBATCH --job-name=gbs_trim
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --partition=large
#SBATCH --time=05:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out 
#SBATCH --error %x.%j.err 

cd /nesi/nobackup/ga03186/Weta_GBS_Batch2/01_stacks_demux_PE/01a_stacks_demux_raw/

module purge
module load TrimGalore/0.6.4-gimkl-2018b

perl /nesi/project/ga03186/scripts/GBS-PreProcess/batch_trim.pl /nesi/project/ga03186/ref/trimgalore_barcodes4.txt
