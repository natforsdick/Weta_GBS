#!/bin/bash -e
#SBATCH --account=ga03186
#SBATCH --job-name=cutadapt
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --partition=large
#SBATCH --time=00:30:00
#SBATCH --array=1 #-96%10 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out 
#SBATCH --error %x.%j.err 

############
# gbs_trim_cutadapt.sl
# Nat Forsdick, 2021-01-29
# Run adapter trimming of raw demux GBS data with cutadapt
###########

###########
# MODULES
module purge
module load cutadapt/2.10-gimkl-2020a-Python-3.8.2
###########

###########
# PARAMS
datadir=/nesi/nobackup/ga03186/Weta_GBS_Batch2/01_stacks_demux/
refdir=/nesi/project/ga03186/scripts/
filtered=/nesi/nobackup/ga03186/Weta_GBS_Batch2/01_stacks_demux/01b_filtered/
trimmed=/nesi/nobackup/ga03186/Weta_GBS_Batch2/01_stacks_demux/01c_trimmed/

start=`date`
echo "Logfile for GBS pipeline run on $start" > logfile.txt
logfile=${datadir}/logfile.txt
###########
cd $datadir
#samplist=`ls -1 *.1.fq.gz | sed 's/\.1\.fq\.gz//'`
samplist=`ls -1 MI20.1.fq.gz | sed 's/\.1\.fq\.gz//'`

echo "Filtering and trimming with cutadapt"
echo "Filtering and trimming with cutadapt" >> $logfile
if [ ! -e ${filtered} ]; then
mkdir -p ${filtered}
mkdir -p ${trimmed}
fi

echo "Filtering summary" > ${filtered}filtering_summary.txt
echo "Trimming summary" > ${trimmed}trimming_summary.txt

#####################
for samp in $samplist

do
echo "Processing $samp" >> $logfile
echo "Filtering $samp"

gunzip ${samp}.1.fq.gz
gunzip ${samp}.2.fq.gz

#####################
cd ${filtered}

echo "$samp" >> filtering_summary.txt
grep '^@' ${datadir}${samp}.1.fq | wc -l >> ${datadir}filtering summary.txt
grep '^@' ${datadir}${samp}.2.fq | wc -l >> ${datadir}filtering summary.txt

grep -B1 -A2 '^TGCAG' ${datadir}${samp}.1.fq | sed '/^--$/d' > ${samp}.1.fq
grep -B1 -A2 '^TGCAG' ${datadir}${samp}.2.fq | sed '/^--$/d' > ${samp}.2.fq

grep '^@' ${samp}.1.fq >> ${datadir}filtering_summary.txt
grep '^@' ${samp}.2.fq >> ${datadir}filtering_summary.txt

cd ${datadir}
######################
echo "Trimming ${samp}" >>${datadir}${logfile}

grep '^@' ${filtered}${samp}.1.fq | wc -l >> ${datadir}trimming_summary.txt
grep '^@' ${filtered}${samp}.2.fq | wc -l >> ${datadir}trimming_summary.txt
cd $datadir$trimmed

echo "$samp" >> trimming summary.txt

cutadapt -a file:${refdir}adapter_complete.fa -m 30 -o ${samp}_trimmed.1.fastq -p ${samp}_trimmed.2.fastq ${filtered}${samp}.1.fq ${filtered}${samp}.2.fq

grep '^@' ${samp}_trimmed.1.fq | wc -l >> trimming_summary.txt
grep '^@' ${samp}_trimmed.2.fq | wc -l >> trimming_summary.txt

cd ${datadir}

echo "$samp processed" 
echo "$samp processed" >> $logfile

gzip ${datadir}${samp}.1.fq
gzip ${datadir}${samp}.2.fq

end=`date`
echo "Filtering and trimming done $end" >> $logfile
