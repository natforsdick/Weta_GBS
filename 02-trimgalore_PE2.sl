#!/bin/bash -e
#SBATCH --job-name=trim_B2
#SBATCH --time=00:06:00
#SBATCH --mem=400M
#SBATCH --cpus-per-task=4
#SBATCH --array=1-96%8
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err

###########
#run_trimgalore_B2.sl
# Nat Forsdick, 2021-03-12
# Script to quality and adapter trim paired-end GBS data
###########

###########
# PARAMS  #
###########
cutadapt=/path/to/cutadapt/2.3-gimkl-2018b-Python-3.7.3/bin/cutadapt
INDIR1=/path/to/01_stacks_demux_PE2/
samplist=/path/to/weta-GBS/Weta_GBS_Batch2_filelist.txt

OUTDIR1=/path/to/01_stacks_demux_PE2/02_trimmed_PE2/
QC1=/path/to/01_stacks_demux_PE2/02_trimmed_PE2/QC/
###########

###########
# MODULES #
###########
module purge
module load TrimGalore/0.6.4-gimkl-2018b FastQC/0.11.9 cutadapt/2.3-gimkl-2018b-Python-3.7.3 
###########

if [ ! -e $OUTDIR1 ]; then
    mkdir -p $OUTDIR1
    fi
if [ ! -e $QC1 ]; then
    mkdir -p ${QC1}
    fi

cd ${INDIR1}
QUERY=`cat ${samplist} | awk -v line=$SLURM_ARRAY_TASK_ID '{if(NR == line) print $1}'`

echo "Trimming ${QUERY}" 
trim_galore --paired --illumina -q 25 --fastqc --fastqc_args "--outdir ${QC1}" --length 50 -o ${OUTDIR1} --path_to_cutadapt $cutadapt -j 4 ${QUERY}.1.fq.gz ${QUERY}.2.fq.gz
echo "Completed ${QUERY}"
