# Whole genome resequencing data

## Intro
This section describes the workflow to infer patterns of species divergence using whole-genome sequencing data of multiple individuals. The analysis assumes the availability of Illumina or DNB sequencing reads covering the entire genome of multiple individuals of the same or related species in fastq format and a single reference genome of the focal or a closely related species. The workflow describes how to preprocess (clean) the raw sequencing data, map the reads to the reference genome, identify single nucleotide polymorphisms (SNPs), perform different types analysis on the genetic relationship between the individuals and visualize the results. 

## Setting the directory structure

Before we start, it’s most convenient to set up an organized directory structure to store different datatypes. Using the bash command 
'mkdir <foldername>' the following directories could be created in your project working directory:
astral		stores the consensus trees produced by astral
bam 		stores the alignment files in BAM format
bed	stores the BED files, being files that list genomic intervals such as predefined genomic windows, proteins, Ultra-Conserved Elements [UCE], etc.. for which trees will be constructed
cleaned 	if reads were first cleaned for adapter contamination and low quality bases, they can be stored here
consensus		stores the diploid consensus genomic sequence for each individual
genome		stores the reference genome sequence in fasta format
reads 		stores the raw sequencing reads in fastq.gz format
samples	contains a text file “samples.txt” listing the sample names (individuals). These sample names should match the basename of the fastq.gz files in the “reads” folder and bam files in the “bam” folder..
scripts		stores your scripts
trees		stores the phylogenetic trees for the different genomic windows
vcf	stores the Variant Call Format [VCF] files, which contain the info on the called SNPs and other types of variation between individuals
multifasta	stores fasta sequence alignments of all individuals for each genomic interval defined in a BED file  
