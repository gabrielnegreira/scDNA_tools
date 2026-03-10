#!/bin/bash
#SBATCH --ntasks=1 --cpus-per-task=64
#SBATCH --time=50:00:00
#SBATCH --job-name=split_bam
#SBATCH --mail-type=BEGIN,END,FAIL

# --- modules ---
module --force purge
module load calcua/2024a calcua/all
module load BEDTools/2.31.1-GCC-13.3.0 SAMtools/1.21-GCC-13.3.0

#set usage helper
usage() {
  cat <<EOF
Usage: sbatch split_bam.sh [options]

Options:
  -i <bam file>           (REQUIRED) path to the input bam file.
  -o <output path>        (REQUIRED) path to the output directory.
  -b <barcodes>           (OPTIONAL) A list of barcodes (one per line) to get reads from. If not set, all barcodes will be split.
  -z <compress>           (OPTIONAL) Decide if the output directory should be compressed (true) into a gzipped file or not (false). 
  -h <help>               Show this help.
EOF
  exit 1
}

# parse options (note new -l)
while getopts ":i:o:b" opt; do
  case "$opt" in
    i) FILE="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    b) BARCODES="$OPTARG" ;;
    z) COMPRESS="$OPTARG" ;;
    h) usage ;;
    \?) echo "Invalid option -$OPTARG" >&2; usage ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Basic validation
if [[ ! ${FILE+x} || ! ${OUT_DIR+x} ]]; then
  echo "ERROR: -i <bam file> and -o <output_path> are required." >&2
  usage
fi

if [[ -n "$BARCODES" ]]; then
  BARCODES_PRINT="$BARCODES"
else
  BARCODES_PRINT="not set"
fi  
#set compress to false if not set.
COMPRESS=${COMPRESS:-false}

if [[ "$COMPRESS" != true ]]; then
    COMPRESS=false
fi

echo ">>> Parameters:"
echo "    Bam file: $FILE"
echo "    Output directory: $OUT_DIR"
echo "    Barcodes: $BARCODES_PRINT"
echo "    Compress: $COMPRESS"

# --- thread planning ---
# Total CPUs allocated (fallback 64 when running outside SLURM)
N_TOTAL=${SLURM_CPUS_PER_TASK:-64}

# Helper: ensure >=1
one_or_more() { v=$1; [ "$v" -lt 1 ] && echo 1 || echo "$v"; }

# Indexing: samtools uses main+workers; use N_TOTAL-1 workers (cap a bit to be kind to FS)
INDEX_WORKERS=$(one_or_more $(( N_TOTAL - 1 )))
# Practical cap (optional): uncomment to limit compression threads during index
# INDEX_WORKERS=$(( INDEX_WORKERS > 32 ? 32 : INDEX_WORKERS ))

# Pipeline case (view | split): avoid oversubscription.
# Two procs run concurrently, each has 1 main thread; budget workers = N_TOTAL-2.
if [ -n "$BARCODES" ]; then
  WORKER_BUDGET=$(( N_TOTAL - 2 ))
  [ "$WORKER_BUDGET" -lt 2 ] && WORKER_BUDGET=2  # guarantee at least 1+1
  # Bias ~75% to split (heavier I/O/compression), rest to view
  SPLIT_WORKERS=$(( (3 * WORKER_BUDGET) / 4 ))
  VIEW_WORKERS=$(( WORKER_BUDGET - SPLIT_WORKERS ))
  SPLIT_WORKERS=$(one_or_more "$SPLIT_WORKERS")
  VIEW_WORKERS=$(one_or_more "$VIEW_WORKERS")
else
  # Single-process case: only split runs; give it N_TOTAL-1 workers
  SPLIT_WORKERS=$(one_or_more $(( N_TOTAL - 1 )))
fi

echo "[INFO] Detected CPUs: $N_TOTAL"
echo "[INFO] Threads -> index: $INDEX_WORKERS; view: ${VIEW_WORKERS:-0}; split: $SPLIT_WORKERS"

# --- prep output dir ---
OUT_DIR="$OUT_DIR/bam_split"
mkdir -p $OUT_DIR

# --- index input BAM if missing ---
if [ ! -f "${FILE}.bai" ]; then
  echo "[INFO] Indexing BAM..."
  samtools index -@ "$INDEX_WORKERS" "$FILE"
else
  echo "[INFO] BAM already indexed; skipping."
fi

# --- split ---
# Output files will be written as <CB>.bam
if [ -n "$BARCODES" ]; then
  echo "[INFO] Splitting only barcodes from: $BARCODES"
  echo "[INFO] Bam file will be splitted into a total of $(cat "$BARCODES" | wc -l) files."
  samtools view -@"$VIEW_WORKERS" -b --tag-file CB:"$BARCODES" "$FILE" \
    | samtools split -@"$SPLIT_WORKERS" -d CB - -f "${OUT_DIR}/%!.bam" -u "${OUT_DIR}/unassigned.bam" --max-split -1
else
  echo "[INFO] Splitting ALL barcodes present in $FILE"
  samtools split -@"$SPLIT_WORKERS" -d CB --max-split -1 -f "${OUT_DIR}/%!.bam" \
    -u "${OUT_DIR}/unassigned.bam" "$FILE"
fi

# --- index outputs (parallel) ---
find "${OUT_DIR}/" -name "*.bam" ! -name "unassigned.bam" -print0 | xargs -0 -n16 -P "$N_TOTAL" -I{} samtools index -@ 1 "{}"

# --- compress ---
if [[ $COMPRESS == true ]]; then

  #now create a sha256 manifest with hashes for all files
  echo "[INFO] Creating sha256 hashes for each file..."
  find "${OUT_DIR}/" -type f -print0 | sort -z | xargs -0 -n16 -P "$N_TOTAL" sha256sum > "manifest.sha256"

  #compress the file
  file_name=$(basename "$OUT_DIR").tar.gzip
  echo "[INFO] compressing output directory..."
  cd ..
  tar -czf $file_name "$OUT_DIR" manifest.sha256
  rm -rf "$OUT_DIR" manifest.sha256
  echo "[DONE] Per-cell BAMs in: $PWD/$file_name"
else
  echo "[DONE] Per-cell BAMs in: $OUT_DIR"
fi  