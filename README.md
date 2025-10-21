# Genomic Workflows for Taxonomy

Authors: Frederik Hendrickx and Carl Vangestel, 2025

License: CC-BY 4.0

Welcome to the **"Genomic Workflows for Taxonomy"** repository. This repository provides modular workflows to infer patterns of species divergence using multilocus sequencing data of multiple individuals. The analysis assumes the availability of Illumina or DNB sequencing reads covering the entire genome of multiple individuals of the same or related species in fastq format and a single reference genome of the focal or a closely related species. The workflow describes how to preprocess (clean) the raw sequencing data, map the reads to the reference genome, identify single nucleotide polymorphisms (SNPs), perform different types analysis on the genetic relationship between the individuals and visualize the results. .

## Workflows
- [Whole Genome Resequencing Data](./whole_genome_resequencing/)
- [RAD/GBS without Reference Genome](./rad_gbs_without_reference/)
- [RAD/GBS with Reference Genome](./rad_gbs_with_reference/)
