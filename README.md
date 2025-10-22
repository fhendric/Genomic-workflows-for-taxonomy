# Genomic Workflows for Taxonomy

Authors: Frederik Hendrickx and Carl Vangestel, 2025

License: CC-BY 4.0

This repository is designed to help taxonomists seamlessly integrate genomic data into their research on species identification, divergence, and speciation. Our goal is to provide clear, step-by-step workflows that take you from raw sequencing reads all the way to advanced analyses of genome-wide patterns of species differentiation.
The workflows are especially tailored for complexes of closely related species, including cryptic species and recent radiations, where subtle genetic differences are key.
Here, you’ll find guidance on:
**- Preprocessing raw sequencing data** to ensure high-quality inputs
**- Mapping reads to reference genomes**
**- Identifying single nucleotide polymorphisms (SNPs)**
**- Analyzing genetic relationships among individuals**
**- Visualizing results to uncover patterns of divergence and relatedness**
Depending on the type of sequencing data that you have available, three different workflows are proposed. The first one assumes that you have sequencing data available from the analysis that cover the entire genome (whole genome resequencing data) and the reference genome of the focal or a related species. The second workflow assumes the availability of Restriction-site Associated DNA (RADseq) or Genotype-by-sequencing (GBS) data available as well as a reference genome. The last option assumes that no reference genome is available.
Depending on the type of sequencing data you have, we offer three tailored workflows:
**1.	Whole-Genome Resequencing** – For studies with sequencing data that cover the entire genome and a reference genome is available for the focal species or a close relative.
**2.	RAD/GBS with Reference Genome** – For datasets generated via Restriction-site Associated DNA (RADseq) or Genotyping-by-Sequencing (GBS), where a reference genome is available.
**3.	RAD/GBS without Reference Genome** – For situations where RAD or GBS data are available but no reference genome exists, providing a de novo approach to analyze genetic variation.
Each workflow is modular and reproducible, guiding you from raw data to meaningful biological insights.


## Workflows
- [Whole Genome Resequencing Data](./whole_genome_resequencing/)
- [RAD/GBS without Reference Genome](./rad_gbs_without_reference/)
- [RAD/GBS with Reference Genome](./rad_gbs_with_reference/)
