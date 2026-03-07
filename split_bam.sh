#!/bin/bash
#SBATCH --ntasks=1 --cpus-per-task=64
#SBATCH --time=50:00:00
#SBATCH --job-name=split_bam
***REMOVED***
***REMOVED***
#SBATCH --mail-type=BEGIN,END,FAIL

# --- modules ---
module --force purge
module load calcua/2024a calcua/all
module load BEDTools/2.31.1-GCC-13.3.0 SAMtools/1.21-GCC-13.3.0

# --- inputs ---
FILE="/scratch/antwerpen/205/vsc20542/projects/scCNV_all_10X_data/cellranger_outputs/SuperMosaic/outs/possorted_bam.bam"
BARCODES="/scratch/antwerpen/205/vsc20542/projects/scCNV_all_10X_data/cellranger_outputs/SuperMosaic/outs/cellranger_barcodes.txt"  # comment this line to split ALL barcodes
OUT_DIR="/scratch/antwerpen/205/vsc20542/projects/scCNV_all_10X_data/cellranger_outputs/SuperMosaic/outs/bam_split"

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
mkdir -p $OUT_DIR
cd "$OUT_DIR" || exit 1

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
    | samtools split -@"$SPLIT_WORKERS" -d CB - -f '%!.bam' -u unassigned.bam --max-split -1
else
  echo "[INFO] Splitting ALL barcodes present in $FILE"
  samtools split -@"$SPLIT_WORKERS" -d CB --max-split -1 -f '%!.bam' \
    -u unassigned.bam "$FILE"
fi

# --- index outputs (parallel) ---
find . -name "*.bam" ! -name "unassigned.bam" -print0 | xargs -0 -n16 -P "$N_TOTAL" -I{} samtools index -@ 1 "{}"

#now create a sha256 manifest with hashes for all files
echo "[INFO] Creating sha256 hashes for each file..."
find . -type f -print0 | sort -z | xargs -0 -n16 -P "$N_TOTAL" sha256sum > ../manifest.sha256

#compress the file
file_name=$(basename "$OUT_DIR").tar.gzip
echo "[INFO] compressing output directory..."
cd ..
tar -czf $file_name "$OUT_DIR" manifest.sha256
rm -rf "$OUT_DIR"
echo "[DONE] Per-cell BAMs in: $PWD/$file_name"