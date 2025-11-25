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

## Input data formats

The standard data format of Illumina or DNB sequencing conducted at a sequencing facility is a pair of fastq files for each sequenced sample, often denoted with a “1” and “2” in the filenames. These two files contain the forward and reverse sequencing read of the sequenced DNA fragments. The order of the reads is identical in the two files, wherein the n’th sequence in the first and second fastq file correspond to the forward and reverse sequencing read of the same DNA fragment. In the fastq format, each sequencing read is represented by four lines. The first line gives the name (@ID) of the sequencing read, the second line shows the actual nucleotide sequence, the third line gives a ‘+’ and the fourth line a quality score of each nucleotide. These quality scores are ASCII coded Phred scores and represent the probability that a base was incorrectly called (se  https://en.wikipedia.org/wiki/Phred_quality_score)
Forward (sampleID_1.fq.gz) file:

### Resequencing data (fastq) 

`@E250063195L1C001R00300000730/1
ATGTCAGATAAATTACTGGTTCCTAAGTTACGAATTAGCTAACCTACTTTTTTCACGTGTTAAAATACAACAATAACATTCATGTACTGCCATTTGCGTCGACCGGCAACGCTAATGTCC
+
FE@EFFFE8FFBE;EFFFFD:FFEF:FEEFFFFEEBFBCFFCEFCFCFCF4EFFFFFFFEFFF>FAFFFFFFBDFFFFFCEFEFDFFFFFFFF:BFFFEFEEFFF=F@FFFF>EFCFCFF`

Reverse (sampleID_2.fq.gz) file:

`@E250063195L1C001R00300000730/2
ATACCAGAAGACGGTCAGTGTCGATATTAAAACCTCTTCTGCCCATCCTATCACGTGATTTGACAGTATAACGGACATTAGCGTTGCCGGTCGACGCAAATGGCAGTACATGAATGTTAT
+
FFFFFFEFCEEGEFEDFF@FDFFFFF?FEFEFFFFFF<BFFFDFFFFF8GFFFFFFDF6;EFFFFF=F;GFEFAFFF9FFFFFFBGFGEFFFEFFFFFF2AGFFFF9GFF%GFF6E(;F?`

For subsequent analyses, it is recommended to use the sampleID as basename of the fastq files and should match the sampleID’s in the “samples.txt” file.  

### Genome (fasta) 

The genome sequence generally comes in regular fasta format with each genomic sequence (contig, scaffold or chromosome) coming in two lines. The first line, starting with “>”, gives the name of the sequence, while the second line gives the actual nucleotide sequence.  
