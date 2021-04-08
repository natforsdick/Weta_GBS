#!/bin/bash -e
#SBATCH --job-name=trim_B2
#SBATCH -A ga03186
#SBATCH --time=00:10:00 # 04:00:00
#SBATCH --mem=300M
#SBATCH --cpus-per-task=12 
#SBATCH --array=1-96%8 #-96%8# # 97-192%16 #1-100%20 #101-192%20 # Tailor to samp # 
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz

###########
#run_trimgalore_B2.sl
# Nat Forsdick, 2021-03-12
# Script to quality and adapter trim paired-end GBS data
###########

###########
# PARAMS  #
###########
cutadapt=/opt/nesi/mahuika/cutadapt/2.3-gimkl-2018b-Python-3.7.3/bin/cutadapt
adapter1="AGATCGGAAGAGC" # SE adapter: "CTGCAAGATCGGAAGAGC"
adapter2="AGATCGGAAGAGC" # Standard Illumina, as determined by automatic detection.
#adapter="CTGCAAGATCGGAAGAGC"
INDIR1=/nesi/nobackup/ga03186/Weta_GBS_Batch2/01_stacks_demux_PE/01a_stacks_demux_PE_raw/
samplist=/nesi/project/ga03186/ref/Weta_GBS_Batch2_filelist.txt

OUTDIR1=/nesi/nobackup/ga03186/Weta_GBS_Batch2/01_stacks_demux_PE/01b_PE_B2_trimmed/
QC1=/nesi/nobackup/ga03186/Weta_GBS_Batch2/01_stacks_demux_PE/01b_PE_B2_trimmed/PE_B2_trimmed_QC/
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
trim_galore --paired -a ${adapter1} -a2 ${adapter2} -q 20 --fastqc --fastqc_args "--outdir ${QC1}" --length 20 -o ${OUTDIR1} --path_to_cutadapt $cutadapt -j 4 ${QUERY}.1.fq.gz ${QUERY}.2.fq.gz
echo "Completed ${QUERY}"
