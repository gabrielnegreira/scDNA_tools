#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=64
#SBATCH --time=12:00:00
#SBATCH --job-name=run_sinto
#SBATCH --mail-type=BEGIN,END,FAIL

#this is a simple bash script to run sinto in the Vlaams supercomputer.

#get the modules
module --force purge
#module load  calcua/2023a calcua/2024a calcua/2025a calcua/all calcua/system 
module load calcua/2024a calcua/all
#module load BEDTools/2.31.1-GCC-13.3.0 SAMtools/1.21-GCC-13.3.0

#load sinto container
unset PYTHONPATH
export PATH="/scratch/antwerpen/205/vsc20542/containers/scDNAseq/bin:$PATH"

#parameters
FILE="/scratch/antwerpen/205/vsc20542/cellranger-dna/outputs/GC1001386_D12/GC1001386_D12/outs/possorted_bam.bam"
REF_PATH="/scratch/antwerpen/grp/aitg/jcdujardin/reference_Ldon_10X/refdata-LdPBQ7G3I2I8X/fasta/genome.fa.fai"
OUTPUT_DIR="/scratch/antwerpen/205/vsc20542/cellranger-dna/outputs/sinto_outs/"

mkdir -p "$OUTPUT_DIR"
cd $OUTPUT_DIR

#create the bed file where read counts will be stored
out_file="frags.bed"
#bedtools makewindows -g "$REF_PATH" -w $BIN_SIZE > "$out_file"

sinto fragments \
  -b "$FILE" \
  -f "$out_file" \
  -t CB \
  --use_chrom "(?i)^Ld"
  -m 30 \
  -p 64 \
  --shift_plus 0 --shift_minus 0  

#sort it and bgzip it
sort -k1,1 -k2,2n "$out_file" > "sort_$out_file"
mv "sort_$out_file" $out_file
