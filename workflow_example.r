#install packages####
install_packages <- FALSE #set to TRUE to install required packages
if(install_packages){
  # Install CRAN packages
  cran_packages <- c(
    "matrixStats",
    "dplyr",
    "tidyr",
    "tibble",
    "forcats",
    "ggplot2",
    "plotly",
    "mgcv",
    "ggrepel",
    "patchwork",
    "hues",
    "mixtools"
  )
  
  install.packages(cran_packages)
  
  # Install Bioconductor package manager if needed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Install Bioconductor packages
  bioc_packages <- c("rhdf5")
  BiocManager::install(bioc_packages)
}

#load libraries and functions####
#get the functions
source("scDNA_functions.R")

#get the data####
##get the count matrix####
count_matrix <- as.matrix(read.csv("count_matrix_example.csv", row.names = 1))
##get the bins metadata####
bins_meta <- read.csv("bins_meta_example.csv")

#build the scDNA obj####
obj <- build_scDNAobj(
  count_matrix = count_matrix, 
  bins_meta = bins_meta)

#perform the analysis####
##we start by distinguishing background droplets from true cells
obj <- tag_true_cells(obj, plot = TRUE)

##remove background dropplets from the object
obj <- subset_cells(obj, cell_or_background == "cell")

## now we do the somy analysis
obj <- tag_outlier_bins(obj) #will append to the bins_meta a flag stating if a bin is an outlier
obj <- correct_counts(obj, vars_to_correct = "gc_content") #will correct counts based on gc content (here shouldn't matter as gc_content is dummy values)
obj <- normalize_counts(obj)
obj <- calc_cells_ICF(obj)
obj <- tag_outlier_cells(obj, vars_to_check = c(ICCV = "upper", ICF_score = "upper"))
obj <- calc_somy(obj, int_method = "GMM") #for gaussian mixuture models set it to "GMM" (but it might get stuck in a loop sometimes, otherwise set it to "round")
obj <- summarise_karyotypes(obj)

#Visualize results####
##visualize somies####
## to visualize somies we can use the `plot_somies` function
plot_somies(obj, matrix_to_plot = "int_somy_matrix")

##visualize individual cells####
## we can use the `plot_cell` function to visualize a specific cell in an object
#let's visualize a cell with low and one with high ICF_score
low_ICF_cell <- obj$metadata$cells_meta %>%
  filter(cell_or_background == "cell" & ICF_score < quantile(ICF_score, 0.2)) %>%
  rownames() %>%
  sample(size = 1)

high_ICF_cell <- obj$metadata$cells_meta %>%
  filter(cell_or_background == "cell" & ICF_score > quantile(ICF_score, 0.95)) %>%
  rownames() %>%
  sample(size = 1)

p1 <- plot_cell(obj, low_ICF_cell)
p2 <- plot_cell(obj, high_ICF_cell)

cowplot::plot_grid(p1, p2, nrow = 1)

##visualize ICF models####
#we can also visualize the model for the ICF_score
plot_ICF_model(obj, cell = high_ICF_cell)
