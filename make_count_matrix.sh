#!/bin/bash
#SBATCH --ntasks=1 --cpus-per-task=64
#SBATCH --time=50:00:00
#SBATCH --job-name=split_bam
#SBATCH --mail-type=BEGIN,END,FAIL

# --- Parameters ---
REF_FILE="/scratch/antwerpen/grp/aitg/jcdujardin/reference_Ldon_10X/refdata-LdPBQ7G3I2I8X/fasta/genome.fa"
BAM_FILES_DIR="/scratch/antwerpen/205/vsc20542/cellranger-dna/outputs/bam_split_with_samtools/per_cell_bams/"
OUT_DIR="/scratch/antwerpen/205/vsc20542/cellranger-dna/outputs/count_matrix/"
BIN_SIZE=5000

# --- other variables ---
JOBS=${SLURM_CPUS_PER_TASK:-1}
THREADS_PER_BAM=1

# --- modules ---
module --force purge
module load calcua/2024a calcua/all
module load BEDTools/2.31.1-GCC-13.3.0 SAMtools/1.21-GCC-13.3.0 parallel/20240722-GCCcore-13.3.0

mkdir -p "${OUT_DIR}" "${OUT_DIR}/tmp"
cd "${OUT_DIR}"

#make sure reference is indexed
if [ ! -f "${REF_FILE}.fai" ]; then
  echo "[INFO] Indexing BAM..."
  samtools faidx "$REF_FILE"
else
  echo "[INFO] reference already indexed; skipping."
fi

#create the bed file for the genomic bins
bedtools makewindows -g "${REF_FILE}.fai" -w "$BIN_SIZE" > "tmp/genome_${BIN_SIZE}bp_bins.bed"
BINS_BED="tmp/genome_${BIN_SIZE}bp_bins.bed"

#create a file which lists all bam files 
find "$BAM_FILES_DIR" -type f -name "*.bam" | sort > tmp/bam_list.tpm

process_one() {
  local BAM="$1" REF_FILE="$2" BINS_BED="$3" OUT_DIR="$4" THREADS="$5"
  local sample tmp
  sample=$(basename "${BAM}" .bam)
  [[ -s "${OUT_DIR}/${sample}.counts.bed" ]] && { echo "[${sample}] skip (done)"; return 0; }

  tmp="${OUT_DIR}/tmp/${sample}_$$"
  mkdir -p "${tmp}"

  echo "[${sample}] start"
  # Minimal threads per small BAM: 1 (reduces overhead)
    samtools sort -n -@ "${THREADS}" -O BAM -T "${OUT_DIR}/tmp.${sample}.nsort" "${BAM}" \
  | samtools fixmate -m -@ "${THREADS}" - - \
  | samtools sort    -@ "${THREADS}" -O BAM -T "${OUT_DIR}/tmp.${sample}.csort" - \
  | samtools markdup -r -@ "${THREADS}" - - \
  | samtools view -b -F 4 -q 30 -@ "${THREADS}" - \
  | bedtools coverage -a "${BINS_BED}" -b stdin -counts -sorted -g "${REF_FILE}.fai" \
    > "${OUT_DIR}/${sample}.counts.bed"

  rm -rf "${tmp}"
  echo "[${sample}] done"
}
export -f process_one
export REF_FILE BINS_BED OUT_DIR THREADS_PER_BAM

# Up to 64 cells in parallel, each single-threaded
parallel --jobs "${JOBS}" --line-buffer process_one {} "${REF_FILE}" "${BINS_BED}" "${OUT_DIR}/tmp/" "${THREADS_PER_BAM}" :::: tmp/bam_list.tpm

echo "Counting done, now merging counts into a single file"

# get the files
files=( tmp/*.counts.bed )

# Strip extension for names
names=$(for f in "${files[@]}"; do basename "$f" .counts.bed; done)

#merge
bedtools unionbedg -i ${files[@]} -header -names $names -empty -g "${REF_FILE}.fai" -filler "NA" \
  > "count_matrix_${BIN_SIZE}bp.tsv"

#remove left overs
rm -rf tmp/
