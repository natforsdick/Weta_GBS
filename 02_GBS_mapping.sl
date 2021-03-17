#!/bin/bash -e
#SBATCH --job-name=mapping_stats_bowtie_B1
#SBATCH -A ga03186
#SBATCH --time=00:03:00 # 04:00:00
#SBATCH --mem=2G # 22G
#SBATCH --cpus-per-task=4 # 6
#SBATCH --array=1-1 #96%10 #-96%8# # 97-192%16 #1-100%20 #101-192%20 # Tailor to samp # - limited to 20 simultaneously.
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
samplist1=/nesi/project/ga03186/ref/B1_filelist.txt
samplist2=/nesi/project/ga03186/ref/Weta_GBS_Batch2_filelist.txt
INDIR=/nesi/nobackup/ga03186/Weta_GBS_Batch1/01_stacks_demux_adap/01b_B1_trimmed/
OUTSAM=/nesi/nobackup/ga03186/Weta_GBS_Batch1/02_bowtie_B1/SAM/
OUTBAM=/nesi/nobackup/ga03186/Weta_GBS_Batch1/02_bowtie_B1/BAM/
OUTSORT=/nesi/nobackup/ga03186/Weta_GBS_Batch1/03_bowtie_B1_sorted/
#INDIR2=/nesi/nobackup/ga03186/Weta_GBS_Batchcombo/
#OUTSAM2=/nesi/nobackup/ga03186/Weta_GBS_Batchcombo/02_stacks_trim_map/SAM/
#OUTBAM2=/nesi/nobackup/ga03186/Weta_GBS_Batchcombo/02_stacks_trim_map/BAM/
#OUTSORT2=/nesi/nobackup/ga03186/Weta_GBS_Batchcombo/03_stacks_sorted/
#samplist2=/nesi/project/ga03186/ref/Weta_GBS_Batch2_filelist1.txt
###########
if [ ! -e ${OUTSAM}   ]; then
    mkdir -p $OUTSAM
    mkdir -p $OUTBAM
    mkdir -p $OUTSORT
        fi

QUERY1=`cat ${samplist1} | awk -v line=$SLURM_ARRAY_TASK_ID '{if(NR == line) print $0}'`
echo $SLURM_ARRAY_TASK_ID
cd $OUTBAM
# Prior to running you must index the reference file with Bowtie2

module load SAMtools/1.10-GCC-9.2.0
module list

#samtools sort $OUTBAM$QUERY.bam -o $OUTSORT$QUERY.sorted.bam
samtools index ${OUTBAM}${QUERY1}.sorted.bam
  
   # Now let's grab some mapping stats:
    
map=$(samtools view -F4 -c ${OUTSORT}$QUERY1.sorted.bam)
unmap=$(samtools view -f4 -c ${OUTSORT}${QUERY1}.sorted.bam)
total=$(($map + $unmap))
perc_mapped=`echo "scale=4;($map/$total)*100" | bc`
         
echo "$QUERY1.bam" >> ${OUTSORT}B1_bowtie_mapping_stats.txt
echo "mapped $map" >> ${OUTSORT}B1_bowtie_mapping_stats.txt
echo "% mapped $perc_mapped" >> ${OUTSORT}B1_bowtie_mapping_stats.txt
echo "unmapped $unmap" >> ${OUTSORT}B1_bowtie_mapping_stats.txt
         
echo "completed $QUERY1"

echo 'Finished processing '${QUERY1}
