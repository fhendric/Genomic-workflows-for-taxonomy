## Mapping Reads to the Reference Genome

To gain meaningful insights from sequencing data, we first need to determine where each read originates from in the genome. This process is called **mapping** or **alignment**, where sequencing reads are matched to a reference genome. There are many mapping tools available, often tailored for specific data types:

- STAR or HISAT for RNA-seq  
- minimap2 for long reads  
- miniprot for protein sequences  

In this guide, we’ll use **BWA-MEM**, a widely used algorithm within the BWA software suite, to align paired-end reads to a reference genome.

---

### 1. Index the Reference Genome

Before mapping, the reference genome must be indexed with BWA. Assuming `$GENOME` points to the reference FASTA file (e.g., `GENOME="./genome/refgenome.fasta"`), indexing is performed using:

```bash
bwa index $GENOME
```

For downstream applications, it is often useful to also generate a .fai index using SAMtools:

```bash
samtools faidx $GENOME
```

