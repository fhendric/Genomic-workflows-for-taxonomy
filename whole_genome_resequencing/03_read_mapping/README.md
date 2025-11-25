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

### 2. Map Paired-End Reads

Assuming that `$READ1` and `$READ2` point to the the forward and reverse fastq files (optionally gzipped) of an individual (e.g. `READ1=”./reads/sample01_1.fq.gz”` and `READ2=”./reads/sample01_2.fq.gz”`), `$SAMPLE` the sampleID of an individual (`SAMPLE=”sample01”`) and `$BAM_OUT` the output file (`BAM_OUT=”./bam/sample01.bam”`) the full command to map paired reads to a reference genome is:

```bash
bwa mem -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" “$GENOME” “$READ1” “$READ2” | samtools view -bS | samtools sort -o “$BAM_OUT”
```

Let's break down the command. The first part, `bwa mem “$GENOME” “$READ1” “$READ2”` does the actual alignment of the reads to the reference genome and outputs a SAM file—a text-based format containing alignment information. Since SAM files are large, they're typically converted to the more compact, binary BAM format using the program **samtools**. This is done on-the-fly using a pipe (|), which passes the SAM output directly into **samtools view**. The `-bS` flags tell **samtools** view to convert from SAM (`-S`) to BAM (`-b`). The `-R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA"`  option in **bwa mem** adds read group information to the header of the resulting BAM file. This includes:


