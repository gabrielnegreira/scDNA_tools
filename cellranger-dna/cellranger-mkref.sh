#!/bin/bash

#---------General parameters--------------------
#this script will take as input a reference genome (FASTA) and will format it to be compatible with cellrangerDNA.
#It checks for contigs < 500kb and fill them with Ns to reach the minimum required length of 500kb.
#It also outputs a contig definition (contig_defs.json) file to with dummy information to be used in cellrangerDNA. 

#parameters

set -euo pipefail
IFS=$'\n\t'

#set usage helper
usage() {
  cat <<EOF
Usage: sbatch cellranger-kmref.sh [options]

Options:
  -i <input>              (REQUIRED) Fasta file containing the reference genome sequence.
  -o <output path>        (OPTIONAL) Output directory. If not set, files will be saved in the currend working directory.
  -h <help>               Show this help.
EOF
  exit 1
}

# parse options (note new -l)
while getopts ":i:o:h" opt; do
  case "$opt" in
    i) REF_FASTA="$OPTARG" ;;
    o) OUTPUTS_DIR="$OPTARG" ;;
    h) usage ;;
    \?) echo "Invalid option -$OPTARG" >&2; usage ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Basic validation
if [[ ! ${REF_FASTA+x} ]]; then
  echo "ERROR: -i <input> is required." >&2
  usage
fi

if [[ ! ${OUTPUTS_DIR+x} ]]; then
  OUTPUTS_DIR="$(pwd)"
fi

#add a slash
OUTPUTS_DIR="${OUTPUTS_DIR}/"

echo ">>> Parameters:"
echo "    Reference Fasta file: $REF_FASTA"
echo "    Output directory: $OUTPUTS_DIR"

#load modules
module --force purge
module load calcua/2025a calcua/all
module load R/4.5.1-gfbf-2025a

#create a function to call the R script
#first define the path to the shell script
script_path=$(dirname $(realpath ${BASH_SOURCE[0]}))
#assuming both scripts are in the same directory, set the full path to the R script:
script_path="${script_path}/format_genome_for_cellranger-dna.R"

#now create the function `format`, which will simply call the R script.
format(){
   Rscript $script_path "$@"
}

#create a function to call mkref
mkref(){
     "/scratch/antwerpen/grp/aitg/tools/cellranger-dna-1.1.0/cellranger-dna" mkref "$@"
}

#first make sure R loads packages from my VSC_DATA directory
R_LIBS_DIR="${VSC_DATA}/Rlibs/${VSC_OS_LOCAL}/${VSC_ARCH_LOCAL}/R-${EBVERSIONR}"
mkdir -p "$R_LIBS_DIR"
echo "R_LIBS_USER=$R_LIBS_DIR" >> $VSC_HOME/.Renviron

#then use the R script to make sure the contigs are > 500kb and generate the contig_defs.json
format $REF_FASTA

#format it for cellrangerDNA
rm -rf refdata-genome_temp
echo 
echo formating the reference genome to work with CellRangerDNA...
mkref "genome_temp.fasta" "contig_defs_temp.json"

#move the generated files to the outputs folder
mkdir -p "$OUTPUTS_DIR"
mv refdata-genome_temp/* "$OUTPUTS_DIR"
rm -rf refdata-genome_temp
rm contig_defs_temp.json
rm genome_temp.fasta

echo "CellrangerDNA-compatible reference files stored in ${OUTPUTS_DIR}"