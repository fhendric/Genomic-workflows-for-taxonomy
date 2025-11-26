## Generating multifasta files for genomic windows

The consensus sequences generated in the previous step can now be combined into **multifasta** files. A multifasta file is simply a FASTA file that contains the consensus sequences of multiple individuals for the same genomic region, and are generally used as input files for phylogenetic anayses. These regions can be consecutive genomic windows, coding regions, etc.... Creating these multifasta files is a two-step process wherein we first make a **BED** file that lists the genomic windows, and then generate the multifasta for each window.

### Step 1: Create a BED file
To define which regions you want to extract, you first need a BED file. A BED file is a plain text file with three columns:
1.	Chromosome or scaffold name
2.	Start position (starting from 0!)
3.	End position

> [!TIP]
> 
> If you want to generate regions of fixed size (e.g., 10 kb windows) at fixed intervals (e.g., every 100 kb), you can use the **makewindows** tool from **bedtools**.
> With **makewindows**, you provide:
> - The reference genome index file (`*.fai`) created earlier with SAMtools (`-g`)
> - The desired window size (`-w`)
> - The step size (`-s`)
> The BED file can then be generated with the following command:
> ```bash
>
> module load BEDtools
> 
> GENOME="./genome/refgenome.fasta"
> BED = “./bed/windows.bed”
> bedtools makewindows -g "$GENOME".fai -w 10000 -s 100000 > "$BED"
> ```

### Step 2: Extract sequences into multifasta files
Once the BED file is ready, we can extract the sequences. The bash script below will create multifasta files for each defined window. Note that the script also includes the generation of the BED file (Step1) and this step should be silenced if you already have a BED file available. 

```bash
# Set input files

GENOME ="./genome/refgenome.fasta"
CONSENSUS="./consensus"
BED = “./bed/windows.bed”
OUTPUT_DIR="./multifasta "
WINDOW_SIZE=10000
STEP_SIZE=100000

# Generate windows [optional if not generated before]

bedtools makewindows -g $GENOME.fai -w $WINDOW_SIZE -s $STEP_SIZE >  $BED
echo "Processing windows..."

# Loop through each window and append sequences to multifasta

while read chrom start end; do
   	win_id="${chrom}_${start}_${end}"
    	out="${OUTPUT_DIR}/${win_id}.fa"

    	# Initialize empty FASTA
    	> "$out"

    	for fasta in $CONSENSUS/*.fa; do
        	sample_name=$(basename "$fasta" .fa)

        	start1=$((start + 1))
        	seq=$(samtools faidx "$fasta" "${chrom}:${start1}-${end}")  # +1 to convert BED to 1-based
        	# Strip header and collapse sequence into one line
        	seq_only=$(echo "$seq" | tail -n +2 | tr -d '\n')

        	echo ">$sample_name" >> "$out"
        	echo "$seq_only" >> "$out"
	done
done < $BED

echo "Multi-FASTA generation completed.”
```
