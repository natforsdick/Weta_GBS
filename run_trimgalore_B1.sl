#!/bin/bash -e
#SBATCH --job-name=mappingB1_adap
#SBATCH -A ga03186
#SBATCH --time=00:03:00 # 04:00:00
#SBATCH --mem=500M
#SBATCH --cpus-per-task=10 # 6
#SBATCH --array=1-96%8 # Tailor to samp # 
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz

###########
#run_trimgalore_B1.sl
# Nat Forsdick, 2021-03-12
# Script to quality and adapter trim GBS data
###########

###########
# PARAMS  #
###########
cutadapt=/opt/nesi/mahuika/cutadapt/2.3-gimkl-2018b-Python-3.7.3/bin/cutadapt
adapter="CTGCAAGATCGGAAGAGC"
INDIR1=/nesi/nobackup/ga03186/Weta_GBS_Batch1/01_stacks_demux_adap/
samplist=/nesi/project/ga03186/ref/B1_barcodes_stacks.txt

OUTDIR1=/nesi/nobackup/ga03186/Weta_GBS_Batch1/01_stacks_demux_adap/01b_B1_trimmed/
QC1=/nesi/nobackup/ga03186/Weta_GBS_Batch1/01_stacks_demux_adap/01b_B1_trimmed/trimmed_QC/
###########

###########
# MODULES #
###########
module purge
module load TrimGalore/0.6.4-gimkl-2018b FastQC/0.11.9 cutadapt/2.3-gimkl-2018b-Python-3.7.3 
###########

cd ${INDIR1}
QUERY=`cat ${samplist} | awk -v line=$SLURM_ARRAY_TASK_ID '{if(NR == line) print $2}'`

echo "Trimming ${QUERY}" 
trim_galore -q 20 -a $adapter --fastqc --fastqc_args "--outdir ${QC1}" --length 20 -o ${OUTDIR1} --path_to_cutadapt $cutadapt -j 4 ${QUERY}.fq.gz
echo "Completed ${QUERY}"
