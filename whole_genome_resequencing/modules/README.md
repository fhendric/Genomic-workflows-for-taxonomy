# Whole genome resequencing data

## Intro
This section describes the workflow to infer patterns of species divergence using whole-genome sequencing data of multiple individuals. The analysis assumes the availability of Illumina or DNB sequencing reads covering the entire genome of multiple individuals of the same or related species in fastq format, as well as a reference genome of the focal or a closely related species. The workflow describes how to preprocess (clean) the raw sequencing reads, map the reads to the reference genome, identify single nucleotide polymorphisms (SNPs), perform different types analysis on the genetic relationship between the individuals and visualize the results. 

## Set-up the directory structure

Before we start, it is most convenient to set up an organized directory structure to store the different datatypes used throughout this workflow. Using the bash command `mkdir <foldername>`, you can create the following directories in your project working directory:


| Directory      | Description |
|----------------|-------------|
| `astral`       | Stores the consensus trees produced by ASTRAL. |
| `bam`          | Stores the alignment files in BAM format. |
| `bed`          | Stores BED files containing genomic intervals such as predefined windows, proteins, Ultra-Conserved Elements (UCEs), etc., for which phylogenetic trees will be constructed. |
| `cleaned`      | Stores cleaned reads (after adapter trimming and quality filtering), if preprocessing is performed. |
| `consensus`    | Stores the diploid consensus genome sequence for each individual. |
| `genome`       | Stores the reference genome sequence in FASTA format. |
| `reads`        | Stores the raw sequencing reads in `fastq.gz` format. |
| `samples`      | Contains a file `samples.txt` listing all sample names (individual IDs). These names should match the basenames of the files in both the `reads` and `bam` directories. |
| `scripts`      | Stores your custom scripts used in the workflow. |
| `trees`        | Stores phylogenetic trees inferred for each genomic window. |
| `vcf`          | Stores Variant Call Format (VCF) files containing SNP and other variant information for all individuals. |
| `multifasta`   | Stores FASTA sequence alignments of all individuals for each genomic interval defined in the BED file. |

## Required programs and tools

| Tool      | Description |
|-----------|-------------|
| ASTRAL    | Phylogenomic inference tool to estimate consensus tree based on a series of gene trees. See also this updated version: (https://github.com/chaoszhang/ASTER?tab=readme-ov-file). |
| BCFTOOLS  | SNP calling and filtering vcf files (https://www.htslib.org/) |
| BEDTOOLS  | Manipulating and extracting sequences using genomic intervals specified in BED files (https://bedtools.readthedocs.io/en/latest/)|
| BWA       | Mapping low-divergent short reads to a reference genome (https://github.com/lh3/bwa) |
| IQTREE    | Maximum likelihood tree analysis (https://iqtree.github.io/) |
| PICARD    | Manipulating BAM/SAM files (https://github.com/broadinstitute/picard) |
| SAMTOOLS  | Suite of programs for interacting with high-throughput sequencing data (https://www.htslib.org/) |
| R packages| ape (Analysis of Phylogenetics and Evolution): functions for analysis of DNA sequences and phylogenetic trees (https://cran.r-project.org/web/packages/ape/index.html)|
|           | phangorn (Phylogenetic Reconstruction and Analysis): Estimation of phylogenetic trees and networks (https://cran.r-project.org/web/packages/phangorn/index.html) |
