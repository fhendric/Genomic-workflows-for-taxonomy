# Whole genome resequencing data

## Intro
This section describes the workflow to infer patterns of species divergence using whole-genome sequencing data of multiple individuals. The analysis assumes the availability of Illumina or DNB sequencing reads covering the entire genome of multiple individuals of the same or related species in fastq format, as well as a reference genome of the focal or a closely related species. The workflow describes how to preprocess (clean) the raw sequencing reads, map the reads to the reference genome, identify single nucleotide polymorphisms (SNPs), perform different types analysis on the genetic relationship between the individuals and visualize the results. The workflow is provided in different components that should be run sequentially.

### Workflow components
- [1. Directory structure, data format and tools](./01_dataformat_and_tools/)
- [2. Clean reads](./02_data_cleaning/)
- [3. Map reads](./03_read_mapping)
- [4. Calling Single Nucleotide Variants (SNPs)](./04_SNP_calling)
- [5. Filter SNPs](./05_SNP_filtering)


> [!TIP]
> Steps 2–5 are the most **computationally demanding** components of the workflow and typically need to be executed on a high-performance computing (HPC) cluster. Running these steps therefore requires not only access to such infrastructure, but also practical experience with cluster environments (e.g., job submission, resource allocation, queue systems). Since this expertise is not always available, these initial steps have been bundled into an integrated workflow called **GenTax**, which is available on the **European Galaxy server** (https://usegalaxy.eu). The Galaxy platform is a widely used, open-source, web-based environment that enables users to run bioinformatics analyses without command-line knowledge. With GenTax, you can simply upload your raw FASTQ files and the reference genome; the workflow then automatically performs Steps 2–5 and produces all required output files, including the variant (SNP) data needed for downstream analyses.
>
> More information and a user manual for running the GenTax workflow can be found here: *[link to be inserted]*.
