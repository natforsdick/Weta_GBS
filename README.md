# Wētā genotyping by sequencing analysis

Developed by Nat Forsdick, 2021. This project is led and funded by [Manaaki Whenua - Landcare Research](https://landcareresearch.co.nz/).

This repo contains scripts used to analyse single- and paired-end GBS data from weta species, including _Deinacrida heteracantha_, _D. fallai_, and _D. mahoenui_.

This work is associated with Forsdick et al., _Population genomic analysis of Mahoenui giant wētā (Deinacrida mahoenui) reveals no reduction in genomic diversity following translocation_, _in progress_, focussing on _D. mahoenui_.

Scripts were originally run on the [NeSI](https://www.nesi.org.nz/) platform via SLURM workload manager, except for R scripts which were run locally. 

The workflow moves through demultiplexing, quality control, and mapping, before processing through Stacks _ref_map_ and _populations_ pipelines after which data are output in formats for analysis via genetics packages such as Adegenet and SNPRelate in R, and STRUCTURE. 

## Software

* Stacks v2.41
* FastQC v0.11.9
* TrimGalore v0.6.4
* cutadapt v2.3
* SAMtools v1.9
* VCFtools v
* PLINK v
* Structure v2.3.4
* R v4.3.1
  * adegenet v2.1.3
  * SNPRelate v1.34.1 
