## Consensus sequence for each individual

Once we obtained our VCF file, we can start to generate a genomic consensus fasta sequence for each individual. There are several ways to generate a consensus sequence, each with their pro’s and con’s. For example, one can substitute nucleotides in the reference sequence by SNPs that are called with high quality. Although this ensures that only high quality SNPs are incorporated, positions without any genotypic information will be represented by the nucleotide in the reference sequence and not reflect the true genotype of that individual on this position. We therefore only want to incorporate positions with genotypic information, with positions without reliable information to be represented by “N”, and heterozygous positions represented as ambiguous IUPAC codes. To ensure the alignment between the sequences is retained, we do not include indels or mnps and use the `project.noindel.minDP5.vcf.gz` VCF generated in previous step. Using bcftools, the consensus sequence for `sample01` can be generated using the following command: 

```bash
# Define genome file
GENOME="./genome/refgenome.fasta"

# Define VCF
VCF_NOINDEL="./vcf/project.noindel.minDP5.vcf.gz "

# Generate consensus sequence
bcftools consensus –-fasta-ref "$GENOME" --missing N --samples sample01 -o ./consensus/sample01.fa "$VCF_NOINDEL"
```

If you want to run the individual analyses in parallel on a cluster, you can use the following script:

```bash
#!/bin/bash

cd ~/project

# Load modules
module load BCFtools
module load SAMtools

# Read sample name from file based on array task ID
SAMPLE=$(sed -n "${PBS_ARRAYID}p" ./samples/samples.txt)

# Input files
VCF_NOINDEL="./vcf/project.noindel.minDP5.vcf.gz" 
GENOME="./genome/refgenome.fasta"
OUT_DIR="./consensus"

# Generate consensus sequence
bcftools consensus --fasta-ref "$GENOME" --missing N --samples "$SAMPLE"  -o "$OUT_DIR/${SAMPLE}.fa" "$VCF_NOINDEL"

# Index consensus FASTA
samtools faidx "$OUT_DIR/${SAMPLE}.fa"
```
