# Mahoenui giant wētā population genomics

<img align="right" src="https://github.com/user-attachments/assets/439d45e4-cb5e-4f52-816f-7cff3b93dfbb">

Developed by Natalie Forsdick, 2021. This project is led by [Manaaki Whenua - Landcare Research](https://landcareresearch.co.nz/), and funded by the New Zealand Ministry of Business, Innovation and Employment through New Zealand's Biological Heritage National Science Challenge and the Strategic Science Investment Fund for Crown Research Institutes, with support from [Genomics Aotearoa](https://www.genomics-aotearoa.org.nz/).

This repo contains scripts used to analyse paired-end genotyping-by-sequencing (GBS) data from Mahoenui giant wētā, _Deinacrida mahoenui_.

This work is associated with the publication: Forsdick et al., 2025. 'Population genomic analysis of Mahoenui giant wētā (_Deinacrida mahoenui_) reveals minimal reduction in genomic diversity following translocation'. Insect Conservation and Diversity. DOI: [10.1111/icad.12810](https://doi.org/10.1111/icad.12810)

Scripts were run on the [NeSI](https://www.nesi.org.nz/) platform via SLURM workload manager, except for `R` scripts which were run locally. 

The workflow moves through demultiplexing, quality control, and mapping, before processing through `Stacks` _ref_map_ and _populations_ pipelines after which data are output in formats for analysis via genetics packages such as `adegenet` and `SNPRelate` in `R`, and `FastSTRUCTURE`. 

## Software

* [Stacks](https://catchenlab.life.illinois.edu/stacks/) v2.65
* [TrimGalore](https://github.com/FelixKrueger/TrimGalore) v0.6.4
  * [FastQC](https://github.com/s-andrews/FastQC) v0.11.9
  * [cutadapt](https://cutadapt.readthedocs.io/en/v2.3/) v2.3
* Bowtie2 v2.3.5
* SAMtools v1.9
* VCFtools v0.1.15
* PLINK v1.09b6.16
* BayeScan v2.1
* [fastStructure]() v1.0
* NeEstimator v2.1
* CLUMPP v1.1.2
* R v4.4.0
  * [SNPfiltR]() v1.0.1
  * [adegenet](https://cran.r-project.org/web/packages/adegenet/index.html) v2.1.10
  * [SNPRelate](https://github.com/zhengxwen/SNPRelate) v1.38.0
  * [pophelper](https://github.com/royfrancis/pophelper) v2.3.1
  * PopGenReport v3.1

## Pipeline

1. [01-stacks_process_radtags.sl](01-stacks_process_radtags.sl) - Demultiplex raw paired-end GBS with Stacks _process_radtags_.
2. [02-trimgalore_PE2.sl](02-trimgalore_PE2.sl) - Trim and adapter removal
3. [03-bowtie-index-ref.sl](03-bowtie-index-ref.sl) - Index the Poor Knights giant wētā reference genome assembly
4. [04-bowtie-align.sl](04-bowtie-align.sl) - Align individual data to the reference genome assembly, collect mapping statistics
5. [05-refmap.sl](05-refmap.sl) - Run Stacks _ref_map.pl_
6. [06-stacks-popns.sl](06-stacks-popns.sl) - Call and filter variants and collect preliminary statistics 
7. [07-export-format.sh](07-export-format.sh) - Convert VCF to various formats for downstream processing
8. [SNP-filtering.Rmd](SNP-filtering.Rmd) - Additional SNP filtering analysis
9. Population genomic analyses:
   * [PCA-popgen.Rmd](PCA-popgen.Rmd) - PCA, Fst, and more
   * [08-faststructure.sl](08-faststructure.sl) - population structure analysis
   * [fastStructure-viz.Rmd](fastStructure-viz.Rmd) - visualisation of the results of population structure analysis

