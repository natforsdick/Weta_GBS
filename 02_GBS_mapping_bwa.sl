#!/bin/bash -e
#SBATCH -J bwa_mappingB1
#SBATCH -A ga03186
#SBATCH --time=00:07:00 
#SBATCH --mem=15G 
#SBATCH --cpus-per-task=12
#SBATCH --array=1-96%4 # Tailor to samp number - limited to 12 simultaneously.
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz

############################
# GBS_mapping_bwa.sl
# Weta GBS mapping pipeline
# Nat Forsdick, 2021-03-04
# This script is to map demultiplexed GBS reads to a reference genome.
# In this case, these are wētā GBS reads mapping to 
# the draft D. fallai reference genome.
############################

###########
# MODULES
module purge
module load  BWA/0.7.17-GCC-9.2.0 SAMtools/1.10-GCC-9.2.0
module list
###########

###########
# PARAMS
fq_list1=/nesi/project/ga03186/ref/B1_filelist.txt
fq_list2=/nesi/project/ga03186/ref/Weta_GBS_Batch2_filelist.txt

reffile=GBS_scaffolds
refdir=/nesi/nobackup/ga03186/reference/
REF=$refdir$reffile

INDIR1=/nesi/nobackup/ga03186/Weta_GBS_Batch1/01_stacks_demux_adap/01b_B1_trimmed/
OUTDIR1=/nesi/nobackup/ga03186/Weta_GBS_Batch1/02_bwa_B1/
INDIR2=/nesi/nobackup/ga03186/Weta_GBS_Batch2/01_stacks_demux_PE/01b_B2_trimmed/
OUTDIR2=/nesi/nobackup/ga03186/Weta_GBS_Batch2/02_bwa_B2/
###########

###########
# MKDIR
#if [ ! -e ${OUTDIR1}  ]; then
#    mkdir -p $OUTDIR1
#    mkdir -p $OUTDIR2
#    fi
###########

###########
# Index reference
#if [ ! -f ${REF}.amb  ]; then
#    cd $refdir
#    echo "The reference has not been indexed. Indexing now"
# Added -b to speed up the indexing - tells it how large the block size should be when processing, and is apparently \
#    optimal when set to GenomeSize / 8. 
#    bwa index -a bwtsw $REF.fasta -b 500000000
#    else
#        echo "BWA index file found" 
#        fi

###########
# Mapping
cd $OUTDIR1
QUERY1=`cat ${fq_list1} | awk -v line=$SLURM_ARRAY_TASK_ID '{if(NR == line) print $0}'`
echo "Task $SLURM_ARRAY_TASK_ID"

echo "Processing $QUERY1"

echo "bwa mem -t 8 ${REF}.fasta ${INDIR1}${QUERY1}.fq.gz"
bwa mem -t 8 ${REF} ${INDIR1}${QUERY1}_trimmed.fq.gz | samtools view -bSh - | samtools sort - -o ${OUTDIR1}${QUERY1}.sorted.bam
samtools index ${OUTDIR1}${QUERY1}.sorted.bam

# Now let's grab some mapping stats:
echo "Getting stats"
map=$(samtools view -F4 -c ${OUTDIR1}$QUERY1.sorted.bam)
unmap=$(samtools view -f4 -c ${OUTDIR1}${QUERY1}.sorted.bam)
total=$(($map + $unmap))
perc_mapped=`echo "scale=4;($map/$total)*100" | bc`

echo "$QUERY1.bam" >> ${OUTDIR1}B1_bwa_mapping_stats.txt
echo "mapped $map" >> ${OUTDIR1}B1_bwa_mapping_stats.txt
echo "% mapped $perc_mapped" >> ${OUTDIR1}B1_bwa_mapping_stats.txt
echo "unmapped $unmap" >> ${OUTDIR1}B1_bwa_mapping_stats.txt

echo "completed $QUERY1"
