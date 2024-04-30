# Wētā genotyping by sequencing analysis

Developed by Nat Forsdick, 2021. This project is led and funded by [Manaaki Whenua - Landcare Research](https://landcareresearch.co.nz/).

This repo contains scripts used to analyse single- and paired-end genotyping-by-sequencing (GBS) data from giant wētā species, including _Deinacrida heteracantha_, _D. fallai_, and _D. mahoenui_.

This work is associated with Forsdick et al., _Population genomic analysis of Mahoenui giant wētā (Deinacrida mahoenui) reveals no reduction in genomic diversity following translocation_, (in progress), focussing on _D. mahoenui_, using a reference genome from _D. fallai_.

Scripts were originally run on the [NeSI](https://www.nesi.org.nz/) platform via SLURM workload manager, except for R scripts which were run locally. 

The workflow moves through demultiplexing, quality control, and mapping, before processing through Stacks _ref_map_ and _populations_ pipelines after which data are output in formats for analysis via genetics packages such as adegenet and SNPRelate in R, and STRUCTURE. 

## Software

* [Stacks](https://catchenlab.life.illinois.edu/stacks/) v2.41

* [TrimGalore](https://github.com/FelixKrueger/TrimGalore) v0.6.4
  * [FastQC](https://github.com/s-andrews/FastQC) v0.11.9
  * [cutadapt](https://cutadapt.readthedocs.io/en/v2.3/) v2.3
* SAMtools v1.9
* VCFtools v
* PLINK v
* [Structure](https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/structure_doc.pdf) v2.3.4
* R v4.3.1
  * [adegenet](https://cran.r-project.org/web/packages/adegenet/index.html) v2.1.3
  * [SNPRelate](https://github.com/zhengxwen/SNPRelate) v1.34.1
  * [pophelper](https://github.com/royfrancis/pophelper) v2.3.1

## Pipeline

1. [stacks_process_radtags.sl](stacks_process_radtags.sl) - Demultiplex raw paired-end GBS with Stacks _process_radtags_.
2. [run_trimgalore_B2.sl](run_trimgalore_B2.sl) - Trim and adapter removal
3. [run_bowtie2_index.sl](run_bowtie2_index.sl) - Index reference genome
4. [02_bowtie_B2.sl](02_bowtie_B2.sl) - Map individual data, collect mapping statistics
5. [03_ref_map.sl](03_ref_map.sl) - Run Stacks _ref_map.pl_
6. [04_stacks_populations_B2.sl](04_stacks_populations_B2.sl) - Call and filter variants, allowing either 30% or 0% missing data, collect preliminary statistics, and output as VCF and PLINK 
7. [05_vcf2adegenet.sl](05_vcf2adegenet.sl) - Convert VCF to PLINK format for conversion to other formats for downstream processing
8. Analysis of final SNP sets in R
   * []() - Discriminant analysis of principal components and more with adegenet
   * []() - Principal component analysis, Fst, and more with SNPRelate
9. [06_structure.sl](06_structure.sl) - Analysis of final SNP sets with STRUCTURE
10. Visualisation of combined STRUCTURE outputs in R with 
