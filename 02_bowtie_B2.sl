#!/bin/bash -e
#SBATCH --job-name=bowtie_B2
#SBATCH -A ga03186
#SBATCH --time=00:40:00 # 04:00:00
#SBATCH --mem=18G # 22G
#SBATCH --cpus-per-task=16 # 6
#SBATCH --array=1-2%2 #-96%8# # 97-192%16 #1-100%20 #101-192%20 # Tailor to samp # - limited to 20 simultaneously.
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
module load Bowtie2/2.3.5-GCC-7.4.0 SAMtools/1.9-GCC-7.4.0
module list
# Newer versions of Bowtie2 throw an error when using --no-unal. Victoria hasn't specified which version she used.
###########

###########
# PARAMS
refdir=/nesi/nobackup/ga03186/reference/
reffile=GBS_scaffolds
ref=$refdir$reffile
samplist=/nesi/project/ga03186/ref/Weta_GBS_Batch2_filelist2.txt
#INDIR=/nesi/nobackup/ga03186/Weta_GBS_Batch1/01_stacks_demux_adap/01b_B1_trimmed/
#OUTSAM=/nesi/nobackup/ga03186/Weta_GBS_Batch1/02_bowtie_B1/SAM/
#OUTBAM=/nesi/nobackup/ga03186/Weta_GBS_Batch1/02_bowtie_B1/BAM/
#OUTSORT=/nesi/nobackup/ga03186/Weta_GBS_Batch1/02_bowtie_B1/02_bowtie_B1_sorted/
INDIR=/nesi/nobackup/ga03186/Weta_GBS_Batch2/01_stacks_demux_PE/01b_PE_B2_trimmed/
OUTSAM=/nesi/nobackup/ga03186/Weta_GBS_Batch2/02_bowtie_PE_B2/SAM/
OUTBAM=/nesi/nobackup/ga03186/Weta_GBS_Batch2/02_bowtie_PE_B2/BAM/
###########

if [ ! -e ${OUTSORT} ]; then
mkdir -p ${OUTSAM}
mkdir -p ${OUTBAM}
mkdir -p ${OUTSORT}
fi

###########

QUERY=`cat ${samplist} | awk -v line=$SLURM_ARRAY_TASK_ID '{if(NR == line) print $0}'`
echo $SLURM_ARRAY_TASK_ID
#`ls -1 *.fq | sed 's/.fq//'`
cd $OUTBAM
# Prior to running you must index the reference file with Bowtie2


# -x=index dir + filename prefix, 
echo 'Aligning '${QUERY}
#### When using Stacks demux inputs:
# Trimmed output format: MI31.1_val_1.fq.gz
# For paired-end reads:
bowtie2 --very-sensitive-local --no-unal --threads 8 \
-x $refdir$reffile \
-1 ${INDIR}$QUERY.1_val_1.fq.gz -2 ${INDIR}$QUERY.2_val_2.fq.gz -S $OUTSAM$QUERY.sam

# For single-end reads:
#bowtie2 --very-sensitive-local --no-unal --threads 8 -x $refdir$reffile -U ${INDIR}${QUERY}.1_val_1.fq.gz -S ${OUTSAM}${QUERY}.sam

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

