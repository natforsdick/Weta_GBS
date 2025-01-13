#!/bin/bash -e
#SBATCH --job-name=mapping
#SBATCH --time=00:30:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=12
#SBATCH --array=1-96%8
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err

###########
# 04-bowtie-align.sl
# Nat Forsdick, 2021-01-21
# This script is to map demultiplexed GBS reads to a reference genome.
# In this case, these are weta Deinacrida mahoenui GBS reads mapping to \
# the draft D. fallai reference genome.
# This is based on the workflow of Victoria Twort. 
# Prior to running, run 03-index-ref.sl to index the reference assembly.
###########

###########
# MODULES
module purge
module load Bowtie2/2.3.5-GCC-7.4.0 SAMtools/1.9-GCC-7.4.0
module list
###########

###########
# PARAMS
refdir=/path/to/ref/
reffile=Deinacrida-fallai-asm-10kb
ref=$refdir$reffile
samplist=/path/to/ref/Weta_GBS_Batch2_filelist.txt
INDIR=/path/to/01_stacks_demux_PE2/02_trimmed_PE2/
OUTSAM=/path/to/03_bowtie/SAM/
OUTBAM=/path/to/03_bowtie/BAM/
###########

if [ ! -e ${OUTBAM} ]; then
mkdir -p ${OUTSAM}
mkdir -p ${OUTBAM}
fi

###########

QUERY=`cat ${samplist} | awk -v line=$SLURM_ARRAY_TASK_ID '{if(NR == line) print $0}'`
echo $SLURM_ARRAY_TASK_ID
#`ls -1 *.fq | sed 's/.fq//'`
cd $OUTBAM
# Prior to running you must index the reference file with Bowtie2

echo 'Aligning '${QUERY}
#### When using Stacks demux inputs:
# Trimmed output format: MI31.1_val_1.fq.gz
# For paired-end reads:
bowtie2 --very-sensitive-local --no-unal --threads 12 \
-x $refdir$reffile \
-1 ${INDIR}$QUERY.1_val_1.fq.gz -2 ${INDIR}$QUERY.2_val_2.fq.gz -S $OUTSAM$QUERY.sam

echo 'Sorting '${QUERY}
samtools view -bSh $OUTSAM$QUERY.sam | samtools sort - -o $OUTBAM$QUERY.sorted.bam
samtools index ${OUTBAM}${QUERY}.sorted.bam

# Now let's grab some mapping stats:

map=$(samtools view -F4 -c ${OUTBAM}$QUERY.sorted.bam)
unmap=$(samtools view -f4 -c ${OUTBAM}${QUERY}.sorted.bam)
total=$(($map + $unmap))
perc_mapped=`echo "scale=4;($map/$total)*100" | bc`

echo "$QUERY.bam" >> ${OUTBAM}B2_bwa_mapping_stats.txt
echo "mapped $map" >> ${OUTBAM}B2_bwa_mapping_stats.txt
echo "perc_mapped $perc_mapped" >> ${OUTBAM}B2_bwa_mapping_stats.txt
echo "unmapped $unmap" >> ${OUTBAM}B2_bwa_mapping_stats.txt

echo "completed $QUERY"

echo 'Finished processing '${QUERY}

