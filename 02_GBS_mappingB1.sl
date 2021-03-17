#!/bin/bash -e
#SBATCH --job-name=mappingB1_adap
#SBATCH -A ga03186
#SBATCH --time=00:08:00 # 04:00:00
#SBATCH --mem=1000M # 22G
#SBATCH --cpus-per-task=12 # 6
#SBATCH --array=1-2#96%12 # Tailor to samp # - limited to 20 simultaneously.
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz

###########
# 02_GBS_mapping.sl
# Nat Forsdick, 2021-01-21
# This script is to map demultiplexed GBS reads to a reference genome.
# In this case, these are weta D. mahoenui GBS reads mapping to \
# the draft D. fallai reference genome.
# This is based on the workflow of Victoria Twort. 
###########

###########
# MODULES
module purge
module load Bowtie2/2.3.5-GCC-7.4.0
module list
# Newer versions of Bowtie2 throw an error when using --no-unal. Victoria hasn't specified which version she used.
###########

###########
# PARAMS
refdir=/nesi/nobackup/ga03186/reference/
reffile=GBS_scaffolds
ref=$refdir$reffile
listdir=/nesi/project/ga03186/ref/
samplist=B1_barcodes_stacks.txt
INDIR=/nesi/nobackup/ga03186/Weta_GBS_Batch1/01_stacks_demux/
OUTSAM=/nesi/nobackup/ga03186/Weta_GBS_Batch1/02_stacks_adap_trim_map/SAM/
OUTBAM=/nesi/nobackup/ga03186/Weta_GBS_Batch1/02_stacks_adap_trim_map/BAM/
OUTSORT=/nesi/nobackup/ga03186/Weta_GBS_Batch1/03_stacks_adap_trim__map_sorted/
###########

if [ ! -e ${OUTSORT} ]; then
mkdir -p ${OUTSAM}
mkdir -p ${OUTBAM}
mkdir -p ${OUTSORT}
fi

###########

QUERY=`cat ${listdir}${samplist} | awk -v line=$SLURM_ARRAY_TASK_ID '{if(NR == line) print $0}'`
echo $SLURM_ARRAY_TASK_ID
#`ls -1 *.fq | sed 's/.fq//'`
cd $INDIR
# Prior to running you must index the reference file with Bowtie2


# -x=index dir + filename prefix, 
echo 'Aligning '${QUERY}
#### When using Stacks demux inputs:
# Trimmed output format: MI31.1_val_1.fq.gz
bowtie2 --very-sensitive-local --no-unal --threads 6 \
-x $refdir$reffile \
-1 $QUERY.1.fq.gz -2 $QUERY.2.fq.gz -S $OUTSAM$QUERY.sam

### When using axe demux inputs:
#bowtie2 --very-sensitive-local --no-unal --threads 6 \
#-x $refdir$reffile \
#-1 Weta_GBS_Batch2_R1_${QUERY}_R1.fastq.gz -2 Weta_GBS_Batch2_R2_${QUERY}_R2.fastq.gz -S $OUTSAM$QUERY.sam

module load SAMtools/1.9-GCC-7.4.0 # CHECK VERSION
module list

echo 'Sorting '${QUERY}
samtools view -S -b $OUTSAM$QUERY.sam > $OUTBAM$QUERY.bam

samtools sort $OUTBAM$QUERY.bam -o $OUTSORT$QUERY.sorted.bam
echo 'Finished processing '${QUERY}

