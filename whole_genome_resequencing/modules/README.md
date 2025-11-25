# Whole genome resequencing data

## Intro
This section describes the workflow to infer patterns of species divergence using whole-genome sequencing data of multiple individuals. The analysis assumes the availability of Illumina or DNB sequencing reads covering the entire genome of multiple individuals of the same or related species in fastq format and a single reference genome of the focal or a closely related species. The workflow describes how to preprocess (clean) the raw sequencing data, map the reads to the reference genome, identify single nucleotide polymorphisms (SNPs), perform different types analysis on the genetic relationship between the individuals and visualize the results. 

## Setting the directory structure

Before we start, it is most convenient to set up an organized directory structure to store the different datatypes used throughout this workflow.  
Using the bash command `mkdir <foldername>`, you can create the following directories in your project working directory:


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
