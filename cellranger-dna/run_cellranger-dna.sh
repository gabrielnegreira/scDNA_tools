#!/bin/bash

#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=64
#SBATCH --time=48:00:00
#SBATCH --job-name=CELLRANGER_DNA
#SBATCH --mail-type=BEGIN,END,FAIL

#this is a simple bash script to run CELLRANGER_DNA in the Vlaams supercomputer.

#parameters

set -euo pipefail
IFS=$'\n\t'

#set usage helper
usage() {
  cat <<EOF
Usage: sbatch run_cellranger-dna.sh [options]

Options:
  -s <sample>             (REQUIRED) ID of the sample, i.e., a string that is shared in the names of the fastq files generated for that sample.
  -f <files path>         (REQUIRED) Path for the directory where the fastq files are located. It does not work recursively, so all files should be in the same directory.   
  -r <reference genome>   (REQUIRED) Path to reference genome directory. It should be a cellrangerDNA-formatted genome, made with `cellranger mkref`.
  -o <output path>        (REQUIRED) Output directory.
  -n <sample name>        (OPTIONAL) Sample name (optional; if omitted it is inferred from -s).
  -c <cell number>        (OPTIONAL) Force cellranger to use a specified number of cells. If not set, cellranger will determine it automatically.
  -h <help>               Show this help.
EOF
  exit 1
}

# parse options (note new -l)
while getopts ":s:f:r:o:n:c:h" opt; do
  case "$opt" in
    s) SAMPLE="$OPTARG" ;;
    f) FASTQ_FILES_DIR="$OPTARG" ;;
    r) REF_GENOME="$OPTARG" ;;
    o) OUTPUTS_DIR="$OPTARG" ;;
    n) SAMPLE_NAME="$OPTARG" ;;
    c) CELL_NUMBER="$OPTARG" ;;
    h) usage ;;
    \?) echo "Invalid option -$OPTARG" >&2; usage ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Basic validation
if [[ -z "${SAMPLE:-}" || -z "${FASTQ_FILES_DIR:-}" || -z "${REF_GENOME:-}" || -z "${OUTPUTS_DIR:-}" ]]; then
  echo "ERROR: -s <sample> -f <files path> -r <reference genome> and -o <output_dir> are required." >&2
  usage
fi

#setting sample name to sample ID if it was not provided
if [[ -z ${SAMPLE_NAME:-} ]]; then
    echo "setting sample name to sample id: ${SAMPLE}"
    SAMPLE_NAME=${SAMPLE}
fi

#setting cell number to `auto` if not set
if [[ -z ${CELL_NUMBER:-} || ! "$CELL_NUMBER" =~ ^[0-9]+$ ]]; then
    CELL_NUMBER=auto
fi

echo ">>> Parameters:"
echo "    Sample: $SAMPLE"
echo "    Fastq files path: $FASTQ_FILES_DIR"
echo "    Reference genome: $REF_GENOME"
echo "    Output directory: $OUTPUTS_DIR"
echo "    Sample name: $SAMPLE_NAME"
echo "    Cell number: $CELL_NUMBER"

#load modules
module --force purge
module load calcua/2024a calcua/all
module load SAMtools/1.21-GCC-13.3.0

#create a function to call cellranger from its install directory
cellranger-dna() {
    "/scratch/antwerpen/grp/aitg/tools/cellranger-dna-1.1.0/cellranger-dna" "$@"
}

#set output directory
mkdir -p "$OUTPUTS_DIR"
cd "$OUTPUTS_DIR"

#run cellranger
echo
echo These are the files stored in the reference genome folder:
tree "$REF_GENOME/"
echo 

#cellranger core command
CMD=(
    cellranger-dna cnv
    --id="$SAMPLE_NAME"
    --fastqs="$FASTQ_FILES_DIR"
    --reference="$REF_GENOME"
    --sample="$SAMPLE"
    --soft-min-avg-ploidy 1.5
    --soft-max-avg-ploidy 5
)

#append `force-cells` if `CELL_NUMBER` is set. 
if [[ "$CELL_NUMBER" != "auto" ]]; then
    CMD+=(--force-cells "$CELL_NUMBER")
fi

#run the command
"${CMD[@]}"

#create an extra file with the list of filtered barcodes only (useful to split the bam file later)
cd "$SAMPLE_NAME/outs"
awk -F',' 'NR>1 {print $1}' per_cell_summary_metrics.csv > cellranger_barcodes.txt