#!/bin/bash -e

ml PLINK/1.09b6.16
INDIR=/nesi/nobackup/ga03048/Weta_GBS_Batch2/03_ref_map_PE/all2/
VCF=weta.snpmiss60.mac3.thin
plink --vcf ${INDIR}${VCF}.vcf --make-bed --aec --out ${INDIR}${VCF}
head ${INDIR}${VCF}.bim

plink --vcf ${INDIR}${VCF}.vcf --aec --recode A --out ${INDIR}${VCF}

plink --vcf ${INDIR}${VCF}.vcf --aec --recode rlist --out ${INDIR}${VCF}
