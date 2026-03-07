#create a function that will try to install libraries if they are not found
install_if_missing <- function(pkgs) {
  # Ensure a CRAN mirror is set (needed for non-interactive Rscript runs)
  if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }

  # Helper: attempt an install and return TRUE/FALSE without stopping the loop
  try_install <- function(expr) {
    ok <- TRUE
    tryCatch(
      { force(expr) },
      error = function(e) ok <<- FALSE,
      warning = function(w) { } # keep warnings, don't fail on them
    )
    ok
  }

  for (pkg in pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("[OK] %s is already installed", pkg))
      next
    }

    message(sprintf("[..] Installing %s from CRAN...", pkg))
    cran_ok <- try_install(install.packages(pkg, quiet = TRUE))

    if (!cran_ok || !requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("[..] CRAN install failed (or package not on CRAN). Trying Bioconductor for %s...", pkg))

      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        message("[..] Installing BiocManager from CRAN...")
        biocmgr_ok <- try_install(install.packages("BiocManager", quiet = TRUE))
        if (!biocmgr_ok) {
          warning(sprintf("[!!] Could not install BiocManager; cannot try Bioconductor for %s", pkg))
          next
        }
      }

      bioc_ok <- try_install(BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE))

      if (!bioc_ok || !requireNamespace(pkg, quietly = TRUE)) {
        warning(sprintf("[!!] Failed to install %s from both CRAN and Bioconductor", pkg))
      } else {
        message(sprintf("[OK] Installed %s from Bioconductor", pkg))
      }
    } else {
      message(sprintf("[OK] Installed %s from CRAN", pkg))
    }
  }
  invisible(TRUE)
}

#check if libraries are installed, if not, install them
libraries <- c("tools", "Biostrings", "jsonlite")
install_if_missing(libraries)

#load the libraries
for(lib in libraries){
  library(lib, character.only = TRUE)
}

#test arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 1){
  stop("Only one argument (the genome fasta file) should be provided")
}
file <- args[[1]]
file_dir <- normalizePath(dirname(file))
file_name <- basename(file)
name_prefix <- tools::file_path_sans_ext(file_name)

#get fasta file
contigs <- readDNAStringSet(file)

#clean contig names
names(contigs) <- gsub(" .*", "", names(contigs))

min_len <- 500000
lens <- width(contigs)

#ask if there are non_nuclear contigs
## Print numbered list
print(paste("Formating", file_name, "file to be compatible with CellrangerDNA."))
print("It will add Ns to contigs < 500kb so they reach the minimum requred size.")
print("It will also export the contig_defs.json configuration file with dummy information.")

contig_names <- names(contigs)
cat("\nThe provided fastq file contains the following contigs:\n")
for (i in seq_along(contig_names)) {
  cat(sprintf("  %d: %s\n", i, contig_names[i]))
}

# Capture user input
cat("\nEnter the numbers corresponding to non-nuclear contigs (comma-separated), or press Enter for none: ")
selection_input <- readLines(con = "stdin", n = 1)

# Parse selection
non_nuclear <- if (nzchar(selection_input)) {

  # Split and trim
  idx <- trimws(unlist(strsplit(selection_input, ",")))
  
  # Convert to numeric
  idx <- suppressWarnings(as.integer(idx))
  
  # Validate indices
  if (any(is.na(idx)) || any(idx < 1) || any(idx > length(contigs))) {
    stop("Invalid selection: please enter valid contig numbers.")
  }
  
  contigs[idx]

} else {
  character(0)
}

if(length(non_nuclear) > 0){
  print("These contigs will be marked as non-nuclear in the contig_defs.json file:")
  non_nuclear 
}else{
  print("All contigs will be considered nuclear contigs in the contig_defs.json file")
  non_nulear <- NULL
}

pad_lengths <- pmax(0, min_len - lens)

pads <- DNAStringSet(
  vapply(
    pad_lengths,
    function(n) if (n > 0) paste(rep("N", n), collapse="") else "",
    character(1)
  )
)

contigs_padded <- xscat(contigs, pads)
names(contigs_padded) <- names(contigs)

#prepared the contig_defs.json file
nuclear <- names(contigs_padded)
non_nuclear <- names(non_nuclear)
nuclear <- nuclear[which(!nuclear %in% non_nuclear)]
first_contig <- nuclear[1]

contig_defs <- list(
  primary_contigs = nuclear,
  non_nuclear_contigs = non_nuclear,
  sex_chromosomes = list(
    male   = setNames(list(2L), first_contig),
    female = setNames(list(2L), first_contig)
  )
)
#save the files
writeXStringSet(contigs_padded, "genome_temp.fasta")
write_json(
  contig_defs,
  path = "contig_defs_temp.json",
  pretty = TRUE,
  auto_unbox = TRUE
)

print(paste("files saved in", getwd()))