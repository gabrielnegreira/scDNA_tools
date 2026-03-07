#scDNA pipeline
#author: Gabriel Heringer Negreira
#Description: This is a group of functions created to analyze the outputs from cellranger-dna of 10X Genomics(https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/what-is-cell-ranger-dna)
# The main input of this script is the cnv_data.h5 file from the cellranger-dna pipeline, which is a HDF5 file containing most of the outputs of the pipeline.
# This HDF5 file is used to extract the count matrix, the cell metadata and the bins metadata. These are combined in a list object (scDNAobj) using the function build_scDNAobj.
# Alternatively, a scDNAobj can be built with a previously built count matrix, cells metadata data frame, and bins metadata data frame. 

#TO DO:
#check dependencies.
#check packages required in each function.
#indicate in the 'plot_karyotypes()' function, which karyotypes are only found in noisy cells.
#indicate in the 'plot_cell()' function if the cells being plotted is noisy or not.

#dependencies
#rhdf5: this is a BioConductor package used to handle HDF5 files. Source: https://bioconductor.org/packages/release/bioc/html/rhdf5.html
#ggplot2
#dplyr
#matrixStats

#libraries
library(matrixStats)
#start by importing the required packages
library(rhdf5)
library(dplyr)
library(tidyr)
library(tibble)
library(forcats)
library(ggplot2)
library(plotly)
library(mgcv)
library(ggrepel)
library(patchwork)
library(hues)
#internal functions####
#check_object
#this function will check if the scDNAobj has the expected structure and the expected columns in the metadata dataframes.
#still not implemented
#this function checks if required packages are installed and return an error if not. Used inside functions.
check_required_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but is not installed. Please install it."))
    }
  }
}

#unlist_2 is an alternative to the base 'unlist' function which allows to chose a character for separating values names. It is used to
#name the bins appropriately in the function build_scDNAobj
unlist_2 <- function(x, sep = "_"){
  names <- names(x)
  names <- rep(names, times = sapply(x, length))
  x <- unlist(x)
  names(x) <- mapply(gsub, names, paste0(names, sep), names(x))
  return(x)
}

#normalize####
#normalize divides a vector by its mean, median or density peak. there is an option to log2 transform it.
normalize <- function(x, method = c("mean", "median", "density_peak"), log2 = FALSE, set_zeros_to = 0){
  method = match.arg(method)
  
  x[x==0] <- set_zeros_to
  
  if(method == "mean"){
    x <- x/mean(x)
  }
  if(method == "median"){
    x <- x/median(x)
  }
  if(method == "density_peak"){
    peak <- density(x)
    peak <- peak$x[which(peak$y == max(peak$y))]
    x <- x/peak
  }
  
  if(log2){
    x <- log2(x)
  }
  return(x)
}


#Mode####
#Mode is a simple function that returns the statistical mode of a vector.
Mode <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux <- ux[tab == max(tab)]
  ux <- ux[1]
  return(ux)
}

#Ratio####
#This is a simple function that will return the ratio between two numbers. 
ratio <- function(x){
  return(x[1]/x[2])
}

#find_knee####
#this function finds the knee point of a ordered graph
find_knee <- function(x, y, val_to_return = c("x", "y", "index")){
  
  val_to_return <- match.arg(val_to_return)
  # Sort the values in decreasing order
  if(missing(y)){
    x <- sort(x, decreasing = TRUE)
    y <- seq_along(x)  
  }else{
    y <- y[order(x, decreasing = TRUE)]
    x <- x[order(x, decreasing = TRUE)]
  }
  
  # Define start and end points of the line
  line_start_x <- x[1]
  line_start_y <- y[1]
  line_end_x <- x[length(x)]
  line_end_y <- y[length(y)]
  
  #Find A, B, and C for the line equation A*x + B*y + C = 0
  A <- line_end_y - line_start_y
  B <- -(line_end_x - line_start_x)
  C <- line_end_x*line_start_y - line_end_y*line_start_x
  
  #basically Ax+By+C describes the *line* between the start and end. Any point (x,y) of that line will return a value of 0. 
  #so if we plug the points of the *curve*, we get actually the distance between that point and the line.
  distances <- vector(length = length(y))
  for(i in seq_along(y)){
    point_x <- x[i]
    point_y <- y[i]
    distance <- (point_x*A)+(point_y*B)+C
    distance <- abs(distance)
    distance <- distance/sqrt(A^2 + B^2)
    distances[i] <- distance
  }
  knee <- which.max(distances)[1] #in case of multiple values, return the first
  knee <- list(x = x[knee], y = y[knee], index = knee)
  return(knee[[val_to_return]])
}

#find_scale_factor####
#this function will calculate a number by which when a vector is multiplied it reaches values closest possible to integers.
find_scale_factor <- function(x, range = c(1.8,5), resolution = 500){
  require(dplyr)
  vec_length <- (range[2]-range[1]) * resolution
  #x <- unlist(x) 
  values_to_try <- seq(range[1], range[2], length = vec_length) #creates a vector containing 'length' values between 'range'
  

  #create a subfunction to calculate distances to integers
  integer_dist <- function(x){
    x <- mean(abs(x-round(x)))
    return(x)
  }
  
  #create the vectors where information will be stored
  tried_factors <- vector(length = vec_length)
  mean_distance_to_integers <- tried_factors
  #try each value and calculate the mean distance to integers. 
  for(i in c(1:length(values_to_try))){ 
    factor_to_try <- values_to_try[i] 
    mean_distance_to_integers[i] <- integer_dist(x*factor_to_try)
    tried_factors[i] <- factor_to_try
  }
  mean_distance_to_integers <- round(mean_distance_to_integers, digits = 3)
  results <- data.frame(tried_factors, mean_distance_to_integers)
  
  #select the scale factor that results in the minimum mean distance to integers.
  #if there are ties, select the lowest scale factor.
  scale_factor <- results %>%
    filter(mean_distance_to_integers == min(mean_distance_to_integers, na.rm = TRUE)) %>%
    filter(tried_factors == min(tried_factors, na.rm = TRUE)) %>%
    getElement("tried_factors")
  
  scale_factor <- list(attempts = results, factor = scale_factor)
  return(scale_factor)
}

#create colors
#this function takes a vector as input and returns a named vector with colors as values and name as elements from the input vector
create_colors <-function(x, 
                         palette = c("default", "all colors", "colorblind friendly", "fancy light", 
                                     "fancy dark", "shades", "tarnish", "pastel", "pimp", 
                                     "intense", "fluo", "red roses", "ochre sand", "yellow lime",
                                     "green mint", "ice cube", "blue ocean", "indigo night", "purple wine"), 
                         seed.use = 123){
  #get required libraries  
  check_required_packages("hues")
  #test inputs
  palette <- match.arg(palette)
  #create pallete parameters
  palette_param <- list(default = c(hmin = 0, hmax = 360, cmin = 30, cmax = 80, lmin = 35, lmax = 80),
                        `all colors` = c(hmin = 0, hmax = 360, cmin = 0, cmax = 100, lmin = 0, lmax = 100),
                        `colorblind friendly` = c(hmin = 0, hmax = 360, cmin = 40, cmax = 70, lmin = 15, lmax = 85),
                        `fancy light` = c(hmin = 0, hmax = 360, cmin = 15, cmax = 40, lmin = 70, lmax = 100),
                        `fancy dark` = c(hmin = 0, hmax = 360, cmin = 8, cmax = 40, lmin = 7, lmax = 40),
                        shades = c(hmin = 0, hmax = 240, cmin = 0, cmax = 15, lmin = 0, lmax = 100),
                        tarnish = c(hmin = 0, hmax = 360, cmin = 0, cmax = 15, lmin = 30, lmax = 70),
                        pastel = c(hmin = 0, hmax = 360, cmin = 0, cmax = 30, lmin = 70, lmax = 100),
                        pimp = c(hmin = 0, hmax = 360, cmin = 30, cmax = 100, lmin = 25, lmax = 70),
                        intense = c(hmin = 0, hmax = 360, cmin = 20, cmax = 100, lmin = 15, lmax = 80),
                        fluo = c(hmin = 0, hmax = 300, cmin = 35, cmax = 100, lmin = 75, lmax = 100),
                        `red roses` = c(hmin = 330, hmax = 20, cmin = 10, cmax = 100, lmin = 35, lmax = 100),
                        `ochre sand` = c(hmin = 20, hmax = 60, cmin = 20, cmax = 50, lmin = 35, lmax = 100),
                        `yellow lime` = c(hmin = 60, hmax = 90, cmin = 10, cmax = 100, lmin = 35, lmax = 100),
                        `green mint` = c(hmin = 90, hmax = 150, cmin = 10, cmax = 100, lmin = 35, lmax = 100),
                        `ice cube` = c(hmin = 150, hmax = 200, cmin = 0, cmax = 100, lmin = 35, lmax = 100),
                        `blue ocean` = c(hmin = 220, hmax = 260, cmin = 8, cmax = 80, lmin = 0, lmax = 50),
                        `indigo night` = c(hmin = 260, hmax = 290, cmin = 40, cmax = 100, lmin = 35, lmax = 100),
                        `purple wine` = c(hmin = 290, hmax = 330, cmin = 0, cmax = 100, lmin = 0, lmax = 40))
  
  x <- sort(unique(x))
  set.seed(seed.use)
  colors <- hues::iwanthue(n = length(x), 
                           hmin = palette_param[[palette]]["hmin"], 
                           hmax = palette_param[[palette]]["hmax"], 
                           lmin = palette_param[[palette]]["lmin"], 
                           lmax = palette_param[[palette]]["lmax"], 
                           cmin = palette_param[[palette]]["cmin"], 
                           cmax = palette_param[[palette]]["cmax"])
  x <- setNames(colors, x)
  return(x)
}

#tag_true_cells####
#this function discriminates true cells from background based on total number of reads associated with each barcode
tag_true_cells <- function(scDNAobj, plot = TRUE, interactive_plot = TRUE, return_obj = c("scDNAobj","plot")){
  return_obj <- match.arg(return_obj)
  #get cells meta 
  cells_meta <- scDNAobj$metadata$cells_meta %>%
    arrange(desc(n_reads)) %>%
    mutate(rank = row_number())
  #find elbow point of rank plot
  n_read_threshold <- find_knee(x = log10(cells_meta$rank), y = log10(cells_meta$n_reads), val_to_return = "y")
  n_read_threshold <- 10^n_read_threshold
  
  cells_meta <- cells_meta %>%
    mutate(cell_or_background = ifelse(n_reads > n_read_threshold, "cell", "background"))
  
  if(plot|return_obj == "plot"){
    plot <- cells_meta %>%
    mutate(cell_or_background = fct_rev(cell_or_background)) %>%
    ggplot(aes(x = rank, y = n_reads, color = cell_or_background))+
      geom_line(linewidth = 1)+
      scale_x_continuous(transform = "log10")+
      scale_y_continuous(transform = "log10")+
      scale_color_manual(values = c(cell = "darkblue", background = "grey"))+
      guides(color = guide_legend(title = NULL))+
      labs(x = "Barcode Rank", y = "Read Counts")+
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 12))
    if(interactive_plot){
      plot <- ggplotly(plot)
    }else{
      last_cell <- cells_meta %>%
        filter(cell_or_background == "cell") %>%
        filter(rank == max(rank)) %>%
        mutate(label = paste0("cell: ", barcode, "\nrank: ", rank, "\nread count: ", n_reads))
      plot <- plot+geom_label_repel(data = last_cell, aes(label = label), show.legend = FALSE, nudge_x = -0.5, nudge_y = -0.5, size = 2)
    }
    show(plot)
  }
  scDNAobj$metadata$cells_meta <- cells_meta %>%
    select(-rank)
  if(return_obj == "plot"){
    return(plot)
  }else{
    return(invisible(scDNAobj))  
  }
}

#subset_cells####
#this function allows to remove cells based on their metadata 
subset_cells <- function(scDNAobj, ...){
  #get cells metadata
  cell_meta <- scDNAobj$metadata$cells_meta
  #filter cells
  scDNAobj$metadata$cells_meta <- scDNAobj$metadata$cells_meta %>%
    filter(...) 
  cells_to_keep <- rownames(scDNAobj$metadata$cells_meta)
  
  #remove cells from count matrices
  scDNAobj$counts <- lapply(scDNAobj$counts, function(x)x <- x[,cells_to_keep])
  
  #remove cells_from models
  if("ICF_models" %in% names(scDNAobj$models)){
    scDNAobj$models$ICF_models <- scDNAobj$models$ICF_models[cells_to_keep]
  }
  return(invisible(scDNAobj))
}

#main functions####
#build_scDNAobj####
#this function will extract the count_matrix, the cells metadata and the bins metadata from the cnv_data.h5 file from cellranger-dna. 
#Alternativelly, the count matrix, cells metadata dataframe, and the bins metadata dataframe can be provided directly. In this case they will be combined in a list object with the structure of a scDNAobj. 
build_scDNAobj <- function(h5, count_matrix = NULL, cells_meta = NULL, bins_meta = NULL, transpose_count_matrix = FALSE, mean_read_length, genome_size){

  #if a h5 file is provided, all the other arguments are ignored. 
  if(!missing(h5)){
    h5file <- H5Fopen(h5)
    #extract the cells metadata
    cells_meta <- h5file$per_cell_summary_metrics %>% 
      as.data.frame()
    rownames(cells_meta) <- cells_meta$barcode
    cells_meta$file_name <- basename(h5)
    #build the bins metadata
    bins_meta <- data.frame(bin = names(unlist_2(h5file$genome_tracks$n_fraction)),
                            gc_content = unlist(h5file$genome_tracks$gc_fraction),
                            mappability = unlist(h5file$genome_tracks$mappability),
                            n_fraction = unlist(h5file$genome_tracks$n_fraction),
                            is_mappable = as.logical(unlist(h5file$genome_tracks$is_mappable))) %>%
                 mutate(chromosome = sub("_.*", "", bin)) %>%
                 mutate(has_Ns = n_fraction > 0) %>%
                 select(bin, chromosome, gc_content, mappability, is_mappable, n_fraction, has_Ns)
    rownames(bins_meta) <- bins_meta$bin
    #extract the count matrix
    count_matrix <- as.matrix(do.call(rbind, h5file$raw_counts))
    count_matrix <- count_matrix[,c(1:nrow(cells_meta))] #this is done to remove nodes and keep only single-cells
    colnames(count_matrix) <- rownames(cells_meta)
    rownames(count_matrix) <- rownames(bins_meta)
  }else{
    if(is.null(cells_meta)){
      cells_meta <- data.frame(barcode = colnames(count_matrix)) %>%
        mutate(cell_id = row_number())
    }
    #test if the provided objects have the expected information
    expected <- c("barcode", "cell_id")
    if(sum(is.na(match(expected, colnames(cells_meta)))) > 0){
      missing <- expected[which(!expected %in% colnames(cells_meta))]
      stop(paste("cells_meta dataframe is missing the following columns:", paste(missing, collapse = ", ")))
    }
    expected <- c("bin", "chromosome", "start", "end", "gc_content", "mappability", "is_mappable")
    if(sum(is.na(match(expected, colnames(bins_meta)))) > 0){
      missing <- expected[which(!expected %in% colnames(bins_meta))]
      stop(paste("bins_meta dataframe is missing the following columns:", paste(expected, collapse = ", ")))
    }
  }
  
  
  #adding extra information to the cells_meta
  cells_meta$n_reads <- colSums(count_matrix)[rownames(cells_meta)]
  #calculate effective_depth_of_coverage if not present in cells metadata
  if(!"effective_depth_of_coverage" %in% colnames(cells_meta)){
    warning("`effective_depth_of_coverage` not present in the cells metadata. Calculating it now...")
    cells_meta$effective_depth_of_coverage <- (cells_meta$n_reads*mean_read_length)/genome_size
  }
  
  #adding extra information to the bins_meta
  bins_meta$mean_raw_counts <- rowMeans(count_matrix)[rownames(bins_meta)]
  
  scDNAobj <- list(counts = c(), metadata = c())
  scDNAobj$counts$raw_counts <- count_matrix
  scDNAobj$metadata$cells_meta <- cells_meta
  scDNAobj$metadata$bins_meta <- bins_meta
  fun_name <- as.character(match.call())[[1]]
  scDNAobj$misc$workflow <- fun_name
  return(invisible(scDNAobj))
}

#correct_counts####
#this function will model the relationship between the specified explanatory variables and the normalized depth of bins and
#will use this model to calculate a correction factor for each bin
correct_counts <- function(scDNAobj, vars_to_correct = c("gc_content", "mappability")){
 
  #Warn if a corrected_counts matrix is already present
  if("corrected_counts" %in% names(scDNAobj$counts)){
    warning("The scDNAobj already contain a corrected_counts matrix, which will be overwritten...")
  }
   
  #get the count matrix
  count_matrix <- scDNAobj$counts$raw_counts
  
  #get the bins metadata
  bins_meta <- scDNAobj$metadata$bins_meta[rownames(count_matrix),]
  
  #create a separate dataframe for the modelling
  to_model <- bins_meta[,vars_to_correct, drop = FALSE]
  
  #remove outliers for each explanatory variable (this is independent from `tag_outlier_bins()`)
  for(var in vars_to_correct){
      outliers <- filter(to_model, !!sym(var) %in% boxplot.stats(!!sym(var))$out) %>%
        rownames()
      to_model <- to_model[which(!rownames(to_model) %in% outliers), , drop = FALSE]
  }
  
  #append the mean count to the model data frame
  to_model$normalized_mean_count <- rowMeans(apply(count_matrix, 2, normalize, method = "mean", log2 = FALSE), na.rm = TRUE)[rownames(to_model)]
  
  #mount the formula
  formula_to_use <- paste0("s(", vars_to_correct, ")")
  formula_to_use <- paste(formula_to_use, collapse = " + ")
  formula_to_use <- as.formula(paste0("normalized_mean_count ~ ", formula_to_use))
  
  #build the model
  message(paste("using GAM modeling to compensate the effect of", paste(vars_to_correct, collapse = " and ")), " on raw counts.")
  model <- mgcv::gam(formula = formula_to_use, data = to_model)
  
  #calculate correction factor for all bins
  bins_meta <- bins_meta %>% 
    mutate(correction_factor = predict(model, bins_meta[,vars_to_correct, drop = FALSE])) %>%
    mutate(correction_factor = median(correction_factor, na.rm = TRUE)/correction_factor)

  #apply correction factors to the count matrix
  count_matrix <- sweep(count_matrix, MARGIN = 1, STATS = bins_meta[rownames(count_matrix), "correction_factor"], FUN = `*`)
  count_matrix[count_matrix < 0] <- 0 #set counts corrected to negative values back to 0
  bins_meta$mean_corrected_counts <- rowMeans(count_matrix)[rownames(bins_meta)]
  
  #return object
  scDNAobj$counts$corrected_counts <- count_matrix
  scDNAobj$metadata$bins_meta <- bins_meta
  scDNAobj$models$bins_correction_factor <- model
  fun_name <- as.character(match.call())[[1]]
  scDNAobj$misc$workflow <- c(scDNAobj$misc$workflow[which(scDNAobj$misc$workflow != fun_name)], fun_name)
  return(invisible(scDNAobj))
}

#normalize_counts####
normalize_counts <- function(scDNAobj){
  #get the count matrix (the matrix to use will depend on if counts were corrected or not)
  matrix_to_use <- c("corrected_counts", "raw_counts") #order determines the priority
  matrix_to_use <- matrix_to_use[which(matrix_to_use %in% names(scDNAobj$counts))][1]
  message(paste0("Using the `", matrix_to_use, "` as count matrix for per-cell normalization..."))
  count_matrix <- scDNAobj$counts[[matrix_to_use]]
  
  #normalize by column 
  count_matrix <- apply(count_matrix, 2, normalize, method = "mean", log2 = FALSE)
  scDNAobj$counts$normalized_counts <- count_matrix
  
  fun_name <- as.character(match.call())[[1]]
  scDNAobj$misc$workflow <- c(scDNAobj$misc$workflow[which(scDNAobj$misc$workflow != fun_name)], fun_name)
  
  return(invisible(scDNAobj))
}


#tag_outlier_bins####
#this function will tag bins with outlier count values. Empty bins will also be flagged as outliers. Additional metadata can be included.
tag_outlier_bins <- function(scDNAobj, additional_vars_to_check = NULL){
  
  #check if `additional_vars_to_check` specifify if outilers detection is for `bottom`, `top`, or `both`
  wrong <- additional_vars_to_check[which(!additional_vars_to_check %in% c("upper", "lower", "both"))]
  if(length(wrong) != 0 | length(names(additional_vars_to_check)) != length(additional_vars_to_check)){
    stop("`additional_vars_to_check` must be a named vector with names corresponding to bin metadata and values being either 'upper', 'lower', or 'both'.")
  }
  
  #get bins_meta
  bins_meta <- scDNAobj$metadata$bins_meta
  
  #check if the vars_to_check are present in the bins metadata
  if(!is.null(additional_vars_to_check)){
    #invert names in addition_vars_to_check to make things easier downstream
    additional_vars_to_check <- setNames(names(additional_vars_to_check), additional_vars_to_check)
    missing_data <- additional_vars_to_check[which(!additional_vars_to_check %in% colnames(bins_meta))]
    if(length(missing_data != 0)){
      stop(paste0("Could not find the following variables in the bins metadata: ",paste0(missing_data, collapse = ", ")))
    }  
  }
  
  #get the count matrix (the matrix to use will depend on if counts were corrected or not)
  matrix_to_use <- c("normalized_counts", "corrected_counts", "raw_counts") #order determines the priority
  matrix_to_use <- matrix_to_use[which(matrix_to_use %in% names(scDNAobj$counts))][1]
  message(paste0("Using the `", matrix_to_use, "` to identify outlier bins..."))
  count_matrix <- scDNAobj$counts[[matrix_to_use]]
  
  #get the bins metadata
  bins_meta <- scDNAobj$metadata$bins_meta  
    
  #tag empty bins based on counts (bins with total counts lower than 1% the number of cells are tagged))
  empty <- rowMeans(count_matrix)
  empty <- empty[empty <= 0.1]  
  col_name <- c(lower = paste0("mean_", matrix_to_use))
  bins_meta$is_empty <- ifelse(rownames(bins_meta) %in% rownames(count_matrix), rownames(bins_meta) %in% names(empty), NA) 
  bins_meta[[col_name]] <- rowMeans(count_matrix)[rownames(bins_meta)]
  
  #now tag outlier bins for each variable (including counts)
  #since tibbles do not store rownames, first we store them
  bins_meta <- rownames_to_column(bins_meta, "rowname")
  vars_to_check <- c(col_name, additional_vars_to_check)
  for(i in seq_along(vars_to_check)){
    var <- vars_to_check[i]
    out_type <- names(vars_to_check)[i]
    out_col_name <- paste0("is_outlier_", var)
    bins_meta <- bins_meta %>%
      group_by(chromosome) %>%
      mutate(!!out_col_name := case_when(out_type == "both" ~ sqrt(!!sym(var)) < boxplot.stats(sqrt(!!sym(var)))$stats[1] | !!sym(var) > boxplot.stats(!!sym(var))$stats[5], #this way NA values will return NA as instead of FALSE
                                         out_type == "lower" ~ sqrt(!!sym(var)) < boxplot.stats(sqrt(!!sym(var)))$stats[1],
                                         out_type == "upper" ~ !!sym(var) > boxplot.stats(!!sym(var))$stats[5])) %>%
      ungroup()
  }  
  bins_meta <- column_to_rownames(bins_meta, "rowname")
  #now we tag any bin which is outlier for at least one variable as outlier
  bins_meta$is_outlier <- rowSums(bins_meta[,grep("is_outlier_|^is_empty$", colnames(bins_meta)), drop = FALSE], na.rm = TRUE) > 0
  scDNAobj$metadata$bins_meta <- bins_meta
  #apend the run to the workflow
  fun_name <- as.character(match.call())[[1]]
  scDNAobj$misc$workflow <- c(scDNAobj$misc$workflow[which(scDNAobj$misc$workflow != fun_name)], fun_name)
  
  message(paste("found", sum(bins_meta$is_outlier), "outlier bins, of which", sum(bins_meta$is_empty), "were empty."))
  return(invisible(scDNAobj))
}

#calc_cell_ICF####
#this function calculates the intra-cromosomal-fluctuation (ICF) of a cell
calc_cells_ICF <- function(scDNAobj, ICF_loess_span = 0.75, max_chromo = 5){
  
  #get the count matrix
  matrix_to_use <- c("normalized_counts", "corrected_counts", "raw_counts") #order determines the priority
  matrix_to_use <- matrix_to_use[which(matrix_to_use %in% names(scDNAobj$counts))][1]
  message(paste0("Using the `", matrix_to_use, "` to calculate cells' ICF..."))
  count_matrix <- scDNAobj$counts[[matrix_to_use]]  
  
  if(matrix_to_use != "normalized_counts"){
   count_matrix <- apply(count_matrix, 2, normalize, method = "mean", log2 = FALSE)
  }
  
  #get the needed metadata
  cells_meta <- scDNAobj$metadata$cells_meta
  bins_meta <- scDNAobj$metadata$bins_meta
  
  #check if tag_outlier_bins was run and if not, run it
  if(!"is_outlier" %in% colnames(bins_meta)){
    warning("outlier bins are not tagged. Running `tag_outlier_bins()` before estimating ICFs...")
    scDNAobj <- tag_outlier_bins(scDNAobj)
    bins_meta <- scDNAobj$metadata$bins_meta
  }
  #remove outlier bins from the count matrix
  bins <- rownames(filter(bins_meta, !is_outlier & is_mappable))
  bins <- bins[bins %in% rownames(count_matrix)]
  count_matrix <- count_matrix[bins,]
  #make sure bins_meta only contains bins found in the count_matrix
  bins_meta <- bins_meta[rownames(count_matrix),]
  
  #ICF is calculated per chromosome. Thus we need to append the chromosome to the count matrix
  count_matrix <- cbind(data.frame(chromosome = bins_meta$chromosome), count_matrix) %>%
    group_by(chromosome) %>%
    mutate(bin_position = row_number()) %>%
    relocate(bin_position, .after = chromosome)
  
  #save the ICF models in the scDNA object
  ICF_models <- list()
  cells_meta$ICF_score <- NA
  cells_meta$ICCV <- NA
  cat("Making ICF models...\n")
  pb <- txtProgressBar(min = 3, max = ncol(count_matrix))
  for(i in c(3:ncol(count_matrix))){
    setTxtProgressBar(pb, i)
    cell <- colnames(count_matrix[,i ,drop = FALSE])
    
    model_data <- count_matrix[,c("chromosome", "bin_position", cell)] %>%
      rename(count := !!sym(cell)) %>%
      group_by(chromosome) %>%
      mutate(ICCV = sd(count)/mean(count))
    
    model_data <- model_data %>%
      group_by(chromosome) %>%
      mutate(loess_fitted = loess(count ~ bin_position, span = ICF_loess_span)$fitted) %>%
      mutate(lm_fitted = mean(count)) %>%
      ungroup()
    ICF_score <- model_data %>%
      group_by(chromosome) %>%
      summarise(ICF = mean(abs(loess_fitted - lm_fitted))) %>%
      slice_max(order_by = ICF, n = max_chromo) %>%
      getElement("ICF") %>%
      median()
      
    #append data to object.
    ICF_models[[cell]]$matrix_used <- matrix_to_use
    ICF_models[[cell]]$model_data <- model_data
    ICF_models[[cell]]$loess_span <- ICF_loess_span
    cells_meta[cell,]$ICF_score <- ICF_score
    cells_meta[cell,]$ICCV <- mean(unique(model_data$ICCV), na.rm = TRUE)
  }
  cat("\n")
  
  #now compensate the ICF by the ICCV
  ##build a model with the non-outlier cells
  to_model <- cells_meta %>%
    filter(!ICF_score %in% boxplot.stats(ICF_score)$out) %>%
    select(ICF_score, ICCV)
  gam_model <- mgcv::gam(formula = ICF_score ~ ICCV, data = to_model)
  
  #now extrapolate the model to all cells
  cells_meta$gam_fit <- predict(gam_model, cells_meta[,"ICCV", drop = FALSE])
  
  #calculate a correction factor and use it to fix the ICF_score
  cells_meta <- cells_meta %>%
    mutate(raw_ICF_score = ICF_score) %>%
    mutate(ICF_score = raw_ICF_score * median(gam_fit)/gam_fit)
  
  cells_meta$gam_fit <- NULL
  #store the outputs in the scDNAobject
  scDNAobj$models$ICF_models <- ICF_models
  scDNAobj$metadata$cells_meta <- cells_meta
  fun_name <- as.character(match.call())[[1]]
  scDNAobj$misc$workflow <- c(scDNAobj$misc$workflow[which(scDNAobj$misc$workflow != fun_name)], fun_name)
  return(invisible(scDNAobj))
}

#tag outlier cells####
#this function will tag cells which have high intra-chromosomal variation (ICF) as 'noisy'. 
tag_outlier_cells <- function(scDNAobj, vars_to_check){
  
  #check if `additional_vars_to_check` specifify if outilers detection is for `bottom`, `top`, or `both`
  wrong <- vars_to_check[which(!vars_to_check %in% c("upper", "lower", "both"))]
  if(length(wrong) != 0 | length(names(vars_to_check)) != length(vars_to_check)){
    stop("`vars_to_check` must be a named vector with names corresponding to bin metadata and values being either 'upper', 'lower', or 'both'.")
  }
  
  #get cells metadata
  cells_meta <- scDNAobj$metadata$cells_meta
  
  #invert names in addition_vars_to_check to make things easier downstream
  vars_to_check <- setNames(names(vars_to_check), vars_to_check)
  
  #check if the requested variables are found in the cells metadata
  missing_data <- vars_to_check[which(!vars_to_check %in% colnames(cells_meta))]
  if(length(missing_data != 0)){
    stop(paste0("Could not find the following variables in the cells metadata: ",paste0(missing_data, collapse = ", ")))
  }  
  
  #now calculate outliers for each variable
  for(i in seq_along(vars_to_check)){
    var <- vars_to_check[i]
    out_type <- names(vars_to_check)[i]
    out_col_name <- paste0("is_outlier_", var)
    cells_meta <- cells_meta %>%
      mutate(!!out_col_name := case_when(out_type == "both" ~ sqrt(!!sym(var)) < boxplot.stats(sqrt(!!sym(var)))$stats[1] | !!sym(var) > boxplot.stats(!!sym(var))$stats[5], #this way NA values will return NA as instead of FALSE
                                         out_type == "lower" ~ sqrt(!!sym(var)) < boxplot.stats(sqrt(!!sym(var)))$stats[1],
                                         out_type == "upper" ~ !!sym(var) > boxplot.stats(!!sym(var))$stats[5]))
  }  
  
  #now we tag any bin which is outlier for at least one variable as outlier
  cells_meta$is_outlier <- rowSums(cells_meta[,grep("is_outlier_", colnames(cells_meta)), drop = FALSE], na.rm = TRUE) > 0
  
  #append the updated cells metadata
  scDNAobj$metadata$cells_meta <- cells_meta
  #append the run to workflow
  fun_name <- as.character(match.call())[[1]]
  scDNAobj$misc$workflow <- c(scDNAobj$misc$workflow[which(scDNAobj$misc$workflow != fun_name)], fun_name)
  return(invisible(scDNAobj))
}

#adjust_somy_distributions####
#this is an internal function which will move the distribution of read counts of each chromosome in such a way that the peak of the distribution is at an integer value.
#this function will calculate a density distribution of somy values for a given chromosome, and adjust that distribution so peaks suround the closest integer
#This is an additional step used to remove chromosome-specific biases.
adjust_somy_distributions <- function(somy_matrix){
  for(i in c(1:nrow(somy_matrix))){
    x <- somy_matrix[i,]
    density <- density(x)
    y_max <- which(density$y == max(density$y))
    peak1 <- density$x[y_max]
    
    adjust_factor <- round(peak1)/peak1 #will calculate how far the peak is from its closest integer
    x <- x*adjust_factor #will multiply x to make the peak of the distribution to became an integer
    somy_matrix[i,] <- x
  }
  return(somy_matrix)
}

#calc_somy####
#this function will calculate the copy number of each chromosome in each cell. 
calc_somy <- function(scDNAobj, matrix_to_use = "corrected_counts", ploidy_limits = c(1.8, 5), int_method = c("GMM", "round"), GMM_type = c("per_chromo", "per_cell"), adjust_distributions = TRUE){
  
  #libraries
  require(dplyr)
  require(matrixStats)
  
  smooth_counts <- FALSE #temporarily smooth_counts is not implemented so set it to false
  
  #check inputs
  int_method <- match.arg(int_method)
  GMM_type <- match.arg(GMM_type)
  
  message(paste("Calculating somies for", nrow(scDNAobj$metadata$cells_meta), "cells..."))
  
  #get the count matrix
  matrix_to_use <- c("normalized_counts", "corrected_counts", "raw_counts") #order determines the priority
  matrix_to_use <- matrix_to_use[which(matrix_to_use %in% names(scDNAobj$counts))][1]
  message(paste0("Using the `", matrix_to_use, "` to calculate cells' ICF..."))
  
  if(matrix_to_use != "normalized_counts"){
    warning("Count matrix was not normalized, normalizing it now...")
    scDNAobj <- normalize_counts(scDNAobj)
  }
  
  count_matrix <- scDNAobj$counts[[matrix_to_use]]
  
  #check if outlier_bins were tagged
  if(!"is_outlier" %in% colnames(scDNAobj$metadata$bins_meta)){
    message("outlier bins were not determined. Running `tag_outlier_bins()` now.")
    scDNAobj <- tag_outlier_bins(scDNAobj)
  }
  
  #get bins and cells metadata
  bins_meta <- scDNAobj$metadata$bins_meta %>%
    filter(!is_outlier, is_mappable)
  cells_meta <- scDNAobj$metadata$cells_meta
  
  #get the count matrix, remove outlier bins and normalize it
  count_matrix <- count_matrix[rownames(bins_meta),]
  
  #smooth the count matrix. Currently not implemented
  if(smooth_counts){
    message("smoothing count matrix with loess...")
    for(chromo in unique(bins_meta$chromosome)){
      bins <- rownames(filter(bins_meta, chromosome == chromo))
      count_matrix[bins,] <- apply(count_matrix[bins,, drop = FALSE], 2, function(x){
        x <- loess(x ~ c(1:length(x)))$fitted
        return(x)
      })
    }
    scDNAobj$counts$smoothed_counts <- count_matrix
  }
  
  #calculate raw somies.
  #raw somies are simply the normalized mean count of each chromosome.
  message("calculating chromosomes' mean...")
  #get the chromosomes
  chromosomes <- unique(bins_meta$chromosome)
  
  #create an empty matrix where values will be stored 
  raw_somy_matrix <- matrix(nrow = length(chromosomes), ncol = ncol(count_matrix))
  rownames(raw_somy_matrix) <- chromosomes
  colnames(raw_somy_matrix) <- colnames(count_matrix)
  
  #this is a second matrix where the standard deviation of each chromosome in each cell is stored
  #not sure if this is needed. Might remove it aftewards
  chromo_sd_matrix <- raw_somy_matrix
  
  #for each chromosome, calculate their mean
  for(i in c(1:length(chromosomes))){
    chromo <- chromosomes[i]
    bins <- filter(bins_meta, chromosome == chromo)$bin
    bins <- bins[bins %in% rownames(count_matrix)]
    raw_somy_matrix[i,] <- colMeans(count_matrix[bins,], na.rm = TRUE)
    chromo_sd_matrix[i,] <- colSds(count_matrix[bins,], na.rm = TRUE)
  }
  
  #add an average sd for each cell in the cells_meta
  cells_meta$mean_chromo_sd <- colMeans(chromo_sd_matrix, na.rm = TRUE)

  #create a column in the cells metadata to store the cells' scale factor
  cells_meta$scale_factor <- NA
  
  #for each cell, normalize the values and multiply by the scale factor
  message("finding cells scale factor...")
  for(cell in colnames(raw_somy_matrix)){
    scale_factor <- find_scale_factor(raw_somy_matrix[,cell], range = ploidy_limits)
    #store the calculation in the scDNA obj
    scDNAobj$scale_factors[[cell]] <- scale_factor
    cells_meta[cell,]$scale_factor <- scale_factor$factor
  }
  #multiply the matrix by the scale factors.
  raw_somy_matrix <- t(t(raw_somy_matrix)*cells_meta$scale_factor)
  
  #adjust somy distributions if needed 
  if(adjust_distributions){
   raw_somy_matrix <- adjust_somy_distributions(raw_somy_matrix)
  }
  
  #now calculate the integer somies
  if(int_method == "GMM"){
    message("applying GMMs to solve integer somies...")
    
    #create the matrix where the integer somies will be stored
    int_somy_matrix <- matrix(ncol = ncol(raw_somy_matrix), nrow = nrow(raw_somy_matrix))
    colnames(int_somy_matrix) <- colnames(raw_somy_matrix)
    rownames(int_somy_matrix) <- rownames(raw_somy_matrix)
    
    if(GMM_type == "per_chromo"){
      #create the list where the GMMs will be stored
      somy_GMMs <- vector(mode = "list", length = nrow(int_somy_matrix))
      names(somy_GMMs) <- rownames(int_somy_matrix)
      #build GMMs for each chromosome
      for(chromo in rownames(int_somy_matrix)){
        GMM <- GMM_to_integers(raw_somy_matrix[chromo,], add_jitter = TRUE)
        int_somy_matrix[chromo,] <- GMM$integers
        somy_GMMs[[chromo]] <- GMM$GMM_model
      }
    }
    
    if(GMM_type == "per_cell"){
      #create the list where the GMMs will be stored
      somy_GMMs <- vector(mode = "list", length = ncol(int_somy_matrix))
      names(somy_GMMs) <- colnames(int_somy_matrix)
      #build GMMs for each chromosome
      for(cell in colnames(int_somy_matrix)){
        
        #get the needed info for that cell
        scale_factor <- cells_meta[cell,]$scale_factor
        possible_somies <- unique(round(raw_somy_matrix[,cell]))
        values <- count_matrix[, cell] * scale_factor
        
        #build the GMM for that cell
        GMM <- GMM_to_integers(values, centers = possible_somies, add_jitter = FALSE)
        
        #calculate the somy of each chromosome
        integers <- data.frame(integer = GMM$integers, chromosome = bins_meta[names(values),]$chromosome) %>%
          group_by(chromosome) %>%
          summarise(integer = round(median(integer))) %>%
          column_to_rownames("chromosome")
        
        int_somy_matrix[,cell] <- integers[rownames(int_somy_matrix),]
        somy_GMMs[[cell]] <- GMM$GMM_model 
      }
    }
    scDNAobj$models$somy_GMMs <- somy_GMMs
  }
  
  if(int_method == "round"){
    int_somy_matrix <- round(raw_somy_matrix)
  }
  
  #add the mean distance to integers to the cells
  cells_meta$mean_int_dist <- colMeans(abs(int_somy_matrix - raw_somy_matrix))[rownames(cells_meta)]
  cells_meta$max_int_dist <- colMaxs(abs(int_somy_matrix - raw_somy_matrix))[rownames(cells_meta)]
  
  #append the outputs to the sDNA object
  scDNAobj$somies$raw_somy_matrix <- raw_somy_matrix
  scDNAobj$somies$chromo_sd_matrix <- chromo_sd_matrix
  scDNAobj$somies$int_somy_matrix <- int_somy_matrix
  scDNAobj$metadata$cells_meta <- cells_meta
  fun_name <- as.character(match.call())[[1]]
  scDNAobj$misc$workflow <- c(scDNAobj$misc$workflow[which(scDNAobj$misc$workflow != fun_name)], fun_name)
  return(invisible(scDNAobj))
}

#summarise_karyotypes####
#this function will summarise the karyotypes found in the object. Only works after a successful run of the calc_somy() function
summarise_karyotypes <- function(scDNAobj, ignore_outlier_cells = TRUE){
  
  #check if int_somy_matrix is available
  if(!"int_somy_matrix" %in% names(scDNAobj$somies)){
    stop("Can't find 'int_somy_matrix' in the scDNA object. Did you run `calc_somy()`?")
  }
  
  #get the needed info
  cells_meta <- scDNAobj$metadata$cells_meta
  int_somy_matrix <- scDNAobj$somies$int_somy_matrix
  
  #remove outlier cells
  if(ignore_outlier_cells){
    if(!"is_outlier" %in% colnames(cells_meta)){
      warning("Outlier cells were not tagged and thus, will not be removed from the karyotype summary. Consider running `tag_outlier_cells` summarising karyotypes.")
    }
    cells <- filter(cells_meta, !is_outlier)
  }else{
    cells <- cells_meta
  }
  
  #add the cells karyotypes 
  cells$karyotype <- apply(int_somy_matrix, 2, paste, collapse = "_")[rownames(cells)]
  
  #create the list of karyotypes
  karyo_list <- cells %>%
    group_by(karyotype) %>%
    summarise(ncells = n(),
              noutliers = sum(is_outlier)) %>%
    ungroup() %>%
    mutate(proportion = ncells/sum(ncells), 
           #sample_name = as.factor(sample_name),
           outlier_only = ncells - noutliers == 0) %>%
    arrange(outlier_only, desc(ncells)) %>%
    mutate(karyo_id = paste0("kar" ,row_number())) %>%
    mutate(karyo_id = factor(karyo_id, levels = karyo_id)) %>%
    ungroup() %>%
    as.data.frame()
    
  rownames(karyo_list) <- karyo_list$karyo_id
  
  #now create a matrix with the somies of each karyotype
  karyo_somies <- matrix(nrow = nrow(karyo_list), ncol = nrow(int_somy_matrix))
  colnames(karyo_somies) <- rownames(int_somy_matrix)
  rownames(karyo_somies) <- karyo_list$karyo_id
  for(i in c(1:nrow(karyo_list))){
    karyo_somies[i,] <- as.integer(unlist(strsplit(karyo_list$karyotype[i], split = "_")))
  }
  #bind_the karyo_somies to the karyo_list
  karyo_list <- cbind(karyo_list, karyo_somies)
  #add the baseline ploidy of the karyotype
  karyo_list$baseline_ploidy <- apply(karyo_somies, 1, Mode)
  #set the cells karyotypes to the karyo_id instead
  cells_meta$karyotype <- NA
  cells_meta[rownames(cells),]$karyotype <- cells$karyotype
  cells_meta$karyotype <- karyo_list$karyo_id[match(cells_meta$karyotype, karyo_list$karyotype)]
  
  #return the objects
  scDNAobj$karyotypes$karyo_list <- karyo_list
  scDNAobj$metadata$cells_meta <- cells_meta
  print(paste0("found ", nrow(karyo_list), " unique karyotypes, of which ", sum(karyo_list$outlier_only), " were found only in outlier cells."))
  fun_name <- as.character(match.call())[[1]]
  scDNAobj$misc$workflow <- c(scDNAobj$misc$workflow[which(scDNAobj$misc$workflow != fun_name)], fun_name)
  return(invisible(scDNAobj))
}

#GMM_to_integers####
#internal function
#this function will apply gaussian mixture models to convert floating points to integers.
#it expects a numeric vector and will use the rounded values to determine k and mu. 
GMM_to_integers <- function(x, centers = "auto", add_jitter = TRUE, show_plot = FALSE, limit_sd = TRUE, use_heuristic = TRUE, heuristic_threshold = 0.3){
  require(mixtools)
  #infer GMM parameters
  counts <- sort(table(round(x)), decreasing = TRUE)
  if(is.numeric(centers)){
    mu <- as.integer(centers)
  }else{
    if(tolower(centers) == "auto"){
      mu <- as.integer(names(counts))
    }else{
      stop("centers must be either a vector of integers or a string 'auto'")
    }
  }
  
  k <- length(mu)
  if(k == 1){
    object <- list(GMM_model = NULL, integers = rep(mu, times = length(x)))
    return(object)
  }
  #to avoid issues with components with a single value, we add some noise to the data:
  to_model <- x
  if(add_jitter){
    to_model <- c(to_model, jitter(to_model))
  }
  if(limit_sd){
    sd.constr <- rep("a", times = k)
  }else{
    sd.constr <- NULL
  }
  #now we create the model
  model <- list(all.loglik = 1, norm_sigma = 100000000)
  #sometimes the EM algorithm does only 1 iteration and this creates a lot of errors. Here, if few iterations were performed, it will repeat the algorithm
  #it will also repeat if one of the gaussians has a standard deviation too high.
  while(length(model$all.loglik) < 5){ 
    warning("Number of iterations too low. Repeating mixture model algortithm.")
    model <- try(normalmixEM(to_model, k = k, mean.constr = mu, sd.constr = sd.constr, epsilon = 1e-08, maxit = 10000, maxrestarts = 1000, ECM = TRUE))
    #if the model returns an error message, re-run it
    if(class(model) == "try-error"){
      model$all.loglik <- 1
    }
  }
  if(show_plot){
    plot_GMM(model)
  }
  #getting the rounded values 
  post <- model$posterior[c(1:length(x)),]
  colnames(post) <- round(model$mu)
  
  #creating 0 columns for integers non-mu integers.
  missing_int <- c(min(mu):max(mu))
  missing_int <- missing_int[!missing_int %in% mu]
  non_mu <- matrix(ncol = length(missing_int), nrow = nrow(post))
  colnames(non_mu) <- missing_int
  non_mu[is.na(non_mu)] <- 0
  post <- cbind(post, non_mu)
  
  integers <- vector(length = length(x))
  for(i in c(1:length(x))){
    floor <- floor(x[i])
    ceiling <- ceiling(x[i])
    #if floor or ceiling are outside the range of possible MUs, adjust
    if(floor < min(round(model$mu))){
      floor <- sort(unique(round(model$mu)))[1]
      ceiling <- sort(unique(round(model$mu)))[2]
    }
    if(ceiling > max(round(model$mu))){
      floor <- sort(unique(round(model$mu)), decreasing = TRUE)[2]
      ceiling <- sort(unique(round(model$mu)), decreasing = TRUE)[1]
    }
    if(floor == ceiling){
      integers[i] <- value
      next
    }else{
      floor <- as.character(floor)
      ceiling <- as.character(ceiling)
      value <- as.integer(names(which.max(post[i,c(floor, ceiling)])))
    }
    integers[i] <- value
  }
  
  #heurstically determine that numbers with distances greater than a threshold to their rounded value should be rounded instead of using the GMMs
  if(use_heuristic){
    #calculate distance to integers
    int_dist <- abs(x - round(x))
    indexes <- which(int_dist < heuristic_threshold)
    integers[indexes] <- round(x[indexes])
  }
  object <- list(GMM_model = model, integers = as.integer(integers))
  return(object)
}


#plot_GMM####
#this is a function to plot a GMM model created with the function normalmixEM from the mixtools package.
plot_GMM <- function(model, histo_bins =100){
  #custom function
  fun_prop <-function(x, mean, sd, proportion){
    proportion * dnorm(x = x, mean = mean, sd = sd)
  }
  plot <- data.frame(x = model$x) %>%
    ggplot(aes(x = x))+
    #geom_density(aes(y = ..density..), bins = histo_bins, fill = "orange", color = "orange", alpha = 0.5)
    geom_histogram(aes(y = after_stat(density)), bins = histo_bins)
    #add the gaussian for each cluster
  for(i in c(1:length(model$mu))){
    plot <- plot + stat_function(geom = "line", fun = fun_prop, args = list(mean = model$mu[i], sd = model$sigma[i], proportion = model$lambda[i]), color = i)
  }
  plot(plot)
}

#plot_loess####
#This function plots a loess model object. 
plot_loess <-function(model){
  require(ggplot2)
  to_plot <- data.frame(x = as.numeric(model$x), y = model$y, fitted = model$fitted)
  plot <- ggplot(to_plot, aes(x = x, y = y))+
    geom_point()+
    geom_line(aes(y = fitted), color = "red", linewidth = 2)+
    labs(x = model$terms[[3]], y = model$terms[[2]])
  return(plot)
}

#plot_cell####
#This function will make a plot summarizing main results for a cell.
plot_cell <- function(scDNAobj, cell, karyotype, metadata_to_show = NULL){
  require(dplyr)
  require(ggplot2)
  require(patchwork)
  
  #check which analysis were done in the object
  outlier_bins_were_tagged = "is_outlier" %in% colnames(scDNAobj$metadata$bins_meta)
  outlier_cells_were_tagged = "is_outlier" %in% colnames(scDNAobj$metadata$cells_meta)
  reads_were_corrected = "corrected_counts" %in% names(scDNAobj$counts)
  reads_were_smoothed = "smoothed_counts" %in% names(scDNAobj$counts)
  cells_were_scaled = "scale_factor" %in% colnames(scDNAobj$metadata$cells_meta)
  raw_somy_was_calculated = "raw_somy_matrix" %in% names(scDNAobj$somies)
  integer_somy_was_calculated = "int_somy_matrix" %in% names(scDNAobj$somies)

  #get the cells metadata
  cells_meta <- scDNAobj$metadata$cells_meta
  #if a karyotype is specified instead of a cell, randomly choose a cell with that karyotype
  if(missing(cell)){
    if(missing(karyotype)){
      message("randomly choosing a cell from the provided scDNA object")
      cell <- sample(cells_meta$barcode, 1)
    }else{
      cell <- cells_meta$barcode[which(cells_meta$karyotype == karyotype)]
      message(paste("randomly choosing a cell from the provided scDNA object with karyotype equal to", karyotype))
      cell <- sample(cell, 1)
    }
  }
  
  #if cell is numeric, convert it to the barcode of the cell in that index
  if(is.numeric(cell)){
    cell <- cells_meta[cell,]$barcode
  }
  
  #get the metadata for the specified cell
  cell_meta <- cells_meta[cell, ]
  #get the bins metadata
  bins_meta <- scDNAobj$metadata$bins_meta 
  
  #Create the plot title if metadata_to_show is in the metadata
  plot_title <- paste0("cell: ", cell_meta$barcode)
  #append additional metadata to title
  if(!is.null(metadata_to_show)){
    if(!metadata_to_show %in% colnames(cells_meta)){
      stop(paste0("Can't find `", paste0(metadata_to_show, collapse = ", "), "` in the cells metadata"))
    }else{
      plot_title <- paste0(plot_title, "; ", paste0(paste(metadata_to_show, cell_meta[metadata_to_show], sep = ": "), collapse = "; "))
    }  
  }
  
  #get the count matrix (the matrix to use will depend on if counts were corrected or not)
  matrix_to_use <- c("normalized_counts", "corrected_counts", "raw_counts") #order determines the priority
  matrix_to_use <- matrix_to_use[which(matrix_to_use %in% names(scDNAobj$counts))][1]
  count_matrix <- scDNAobj$counts[[matrix_to_use]]
  
  cell_reads <- data.frame(bin = rownames(count_matrix), value = count_matrix[,cell])
  colnames(cell_reads)[2] <- matrix_to_use
  
  #bind bin metadata
  to_plot <- cbind(cell_reads, bins_meta[rownames(cell_reads), -which(colnames(bins_meta) %in% colnames(cell_reads))])
  to_plot$bin_pos <- c(1:nrow(to_plot))
  
  #remove non mappable bins
  to_plot <- to_plot[which(to_plot$is_mappable),]
  #remove outlier bins
  if(outlier_bins_were_tagged){
    to_plot <- to_plot[which(to_plot$is_outlier == FALSE),]
  }
  
  #create the base plot
  x_axis <- to_plot %>%
    group_by(chromosome) %>%
    summarise(start = min(bin_pos), end = max(bin_pos), middle = median(bin_pos))
  p1 <- to_plot  %>%
    ggplot(., aes(x = bin_pos, y = .data[[colnames(.)[2]]], color = chromosome))+
    geom_point(size = 0.75)+
    scale_color_manual(values = rep(c("orange", "black"), times = 1000))+
    scale_x_continuous(breaks = x_axis$middle, labels = x_axis$chromosome)+
    scale_y_continuous(breaks = c(1:1000))+
    guides(color = "none")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor.x = element_blank())
  
  
  #now add raw_somies or integer_somies if they are available
  if(integer_somy_was_calculated | raw_somy_was_calculated){
    p1$data[[matrix_to_use]] <- p1$data[[matrix_to_use]]*cell_meta$scale_factor 
    somies <- c("int_somy_matrix", "raw_somy_matrix")
    somies <- somies[which(somies %in% names(scDNAobj$somies))][1]
    somies <- scDNAobj$somies[[somies]][,cell, drop = TRUE]
    somies <- p1$data %>%
      mutate(somy = somies[chromosome])
    p1 <- p1+scale_color_manual(values = rep(c("#F3D18F", "#292929"), times = 1000))
    p1 <- p1+geom_line(data = somies, aes(y = somy, group = chromosome), color = "red")
  }
  
  #now if the cell scale factor was determined, plot it
  if(cells_were_scaled){
    p2 <- scDNAobj$scale_factors[[cell]]$attempts %>%
      ggplot(aes(x = tried_factors, y = mean_distance_to_integers))+
      geom_path()+
      geom_vline(xintercept = scDNAobj$scale_factors[[cell]]$factor, color = "red", linetype = "dashed")+
      labs(x = "scale_factor")+
      ggtitle(paste0("scale factor: ", round(scDNAobj$scale_factors[[cell]]$factor, 2)))+
      theme(plot.title = element_text(size = 11))
  }else{
    p2 <- ggplot()+
      theme_void()+
      ggtitle("No scale factor determined")
  }
  
  
  p3 <- cells_meta %>%
    mutate(color = ifelse(barcode == cell, "cell", 
                          ifelse(is_outlier,  "outlier",  "other"))) %>%
    mutate(color = factor(color, levels = c("cell", "other", "outlier"))) %>%
    arrange(desc(color)) %>%
    ggplot(aes(x = ICCV, y = ICF_score, color = color))+
    geom_point()+
    guides(color = "none", size = "none")+
    ggtitle(paste0("Coverage: ", round(cell_meta$effective_depth_of_coverage,2), "x; ICF_score: ", round(cell_meta$ICF_score, 2)))+
    scale_color_manual(values = c(cell = "red", outlier = "lightgrey", other = "#c9a9a9"))+
    theme(plot.title = element_text(size = 11))
  
 
  plot <- (p3+p2)/p1
  plot <- plot +plot_annotation(title = plot_title)
  
  return(plot)
}

#plot_ICF_model
plot_ICF_model <- function(scDNAobj, cell){
  
  #get the cells metadata
  cells_meta <- scDNAobj$metadata$cells_meta
  #if a karyotype is specified instead of a cell, randomly choose a cell with that karyotype
  if(missing(cell)){
      message(paste("randomly choosing a cell from the provided scDNA object"))
      cell <- sample(rownames(cells_meta), size = 1)
  }
  
  #if cell is numeric, convert it to the barcode of the cell in that index
  if(is.numeric(cell)){
    cell <- cells_meta[cell,]$barcode
  }
  
  cell_meta <- cells_meta[cell,]
  
  plot <- scDNAobj$models$ICF_models[[cell]]$model_data %>%
    mutate(bin = row_number())%>%
    ggplot(aes(x = bin, y = count, group = chromosome))+
    geom_point()+
    geom_line(aes(y = loess_fitted), color = "blue")+
    geom_line(aes(y = lm_fitted), color = "red")+
    ggtitle(label = paste0("Cell: ", cell_meta$barcode, "; ICF_score: ", round(cell_meta$ICF_score, digits = 3)))+
    guides(col = "none")+
    #scale_color_manual(values = rep(c("black", "orange"), times = 100))+
    facet_wrap(vars(chromosome), scales = "free")
  
  return(plot)
}


#plot_somies####
#this function will plot the matrix of raw or intenger somies
plot_somies <- function(scDNAobj, matrix_to_plot = "int_somy_matrix", annotations = NULL, annotation_colors = NULL, remove_outliers = FALSE){
  require(ComplexHeatmap)
  require(circlize)
  require(dplyr)
  
  data_to_plot <- scDNAobj$somies[[matrix_to_plot]]
  cells_meta <- scDNAobj$metadata$cells_meta[colnames(data_to_plot),]
  
  #getting the unique integers
  breaks <- sort(unique(as.integer(data_to_plot)))
  color_palette <- heat_col(breaks)
  color_palette <- colorRamp2(breaks, color_palette)
  
  if(remove_outliers){
    data_to_plot <- data_to_plot[,scDNAobj$metadata$cells_meta$barcode[which(!scDNAobj$metadata$cells_meta$is_outlier)]]
    total_cells <- nrow(scDNAobj$metadata$cells_meta)
    outlier_cells <- nrow(filter(scDNAobj$metadata$cells_meta, is_outlier))
    title <- paste0(ncol(data_to_plot), " cells. ", outlier_cells, " outlier cells not shown.")
    
  }else{
    total_cells <- nrow(scDNAobj$metadata$cells_meta)
    outlier_cells <- nrow(filter(scDNAobj$metadata$cells_meta, is_outlier))
    title <- paste0(total_cells, " cells, of which ", outlier_cells, " are tagged as outliers")
  }
  
  if(!is.null(annotations)){
    if(!remove_outliers){
      annotations <- c(annotations[which(annotations != "is_outlier")], "is_outlier")
    }
    
    annotations <- scDNAobj$metadata$cells_meta[colnames(data_to_plot),annotations, drop = FALSE]
    
    if(is.null(annotation_colors)){
      #select only non-numeric columns
      cols <- sapply(annotations, is.numeric)
      cols <- names(cols[which(cols == FALSE)])
      annotations_unique <- lapply(annotations[,cols, drop = FALSE], unique)
      annotation_colors <- create_colors(unlist(annotations_unique))
      annotation_colors <- lapply(annotations_unique, function(x)annotation_colors[as.character(x)])
    }
    annotations <- HeatmapAnnotation(df = annotations, col = annotation_colors)
  }
  
  plot <- Heatmap(data_to_plot,
          show_column_names = FALSE,
          cluster_rows = FALSE,
          name = "somy",
          #cluster_column_slices = FALSE, 
          top_annotation = annotations,
          col = color_palette,
          use_raster = FALSE,
          column_title = title,
          column_title_side = "bottom",
          na_col = "grey")
  
  return(plot)
}

#plot karyotypes####
#this function will plot the karyotypes found in the scDNAobj.
plot_karyotypes <- function(scDNAobj){
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(patchwork)
  
  #get the needed objects
  karyo_list <- scDNAobj$karyotypes$karyo_list
  int_somy_matrix <- scDNAobj$somies$int_somy_matrix
  sample_name <- unique(scDNAobj$metadata$cells_meta$sample_name)
  
  #convert the matrix to tibble
  to_plot <- karyo_list %>%
    pivot_longer(cols = c(rownames(int_somy_matrix)), names_to = "chromosome", values_to = "somy") %>%
    mutate(somy = as.integer(somy)) %>%
    mutate(chromosome = factor(chromosome, levels = rev(levels(factor(chromosome)))))
  
  #set the colors
  colors <- heat_col(as.integer(unique(to_plot$somy)))
  
  #create a title for the plot with additional information
  title <- paste0("Sample: ", sample_name, "; ", nrow(karyo_list), " karyotypes found of which ", sum(karyo_list$outlier_only), " are only found in outlier cells.")
  #now first build the bar plot
  max_y <- max(karyo_list$proportion) * 1.23
  p1 <- ggplot(karyo_list, aes(x = karyo_id, y = proportion, label = ncells))+
    geom_col(fill = "black")+
    geom_text(vjust = -0.3, size = 3.7, color = "white")+
    geom_text(vjust = -0.3, size = 3.5, color = "black")+
    theme_minimal()+
    scale_x_discrete(expand = c(0.025,0))+
    scale_y_continuous(limits = c(0, max_y), breaks = seq(0.1, 1, by = 0.1))+
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = margin(0, 0, 0, 0))
  
  #now make the heatmap
  p2 <- ggplot(to_plot, aes(x = karyo_id, y = chromosome, fill = as.factor(somy)))+
    geom_tile(width = 0.98, col = "#2a5686")+
    scale_x_discrete(expand = c(0.025,0))+
    scale_y_discrete()+
    labs(fill = "Somy")+
    scale_fill_manual(values = colors)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.margin = margin(0,0,0,0),
          panel.background = element_blank())
  
  #aling both plots
  plot <- p1/p2+plot_layout(ncol = 1, heights = c(0.2,0.8))+plot_annotation(title = title)
  
  return(plot)

}

#heat_col####
#this is a function used to generate the colors for the heatmaps
heat_col <- function(breaks){
  require(colorspace)
  #set the color palette
  colors <- c("#001221", "#002342", "#002342", 
              "#014175", "#035ba3", "#00c3ff", 
              "#00ffee", "#33ff00", "#ccff00", 
              "#fffa00","#ffa600", "#D73027", 
              "#A50026", "#541b1b", "#4d0600")
  
  #get the maximum value
  breaks <- as.integer(breaks)
  n <- max(breaks)
  #by generating 9 colors with this palete, we make sure that values between 0 and 8 are mapped to discrete colors.
  colors <- colorRampPalette(colors)(9)
  #if it asks for more than 9 colors, it will generate the additional colors by darkening the last color
  if(n >= 9){
    for(i in c(10:(n+1))){
      colors[i] <- darken(colors[i-1], amount = 0.4)
    }
  } 
  #now covert values to indices. value 0 should be mapped to the first color, and so on.
  breaks <- breaks+1
  #map the colors to the values
  colors <- colors[breaks]
  names(colors) <- breaks-1
  return(colors)
}


#karyo_dist####
#this is an internal function which will calculate the copy number distance between karyotypes from a karyo_list data frame.
karyo_dist <- function(karyo_somies){
  
  karyo_somies <- as.matrix(karyo_somies)
  
  #get the karytotypes and their ids
  ids <- rownames(karyo_somies)
  
  #creating the matrix where the information will be stored
  diff_matrix <- matrix(nrow = length(ids), ncol = length(ids))
  colnames(diff_matrix) <- ids
  rownames(diff_matrix) <- ids
  
  #comparing each karyotype vs each other
  for(i in c(1:nrow(karyo_somies))){
    diff_matrix[i,] <- rowSums(t(abs(t(karyo_somies) - as.integer(karyo_somies[i,]))))
  }
  
 #return the matrix
  return(diff_matrix)
}


#karyo_network####
#this function will build a minimum spanning tree to visualize the similarity relationship between karyotypes in a scDNA object.
#obs: root not implemented yet
karyo_network <- function(scDNAobj, root = NULL, method = "mst"){
  #libraries
  require(igraph)
  require(visNetwork)
  require(pegas)
  require(hues)
  
  #get the needed data
  karyo_list <- scDNAobj$karyotypes$karyo_list
  int_somy_matrix <- scDNAobj$somies$int_somy_matrix
  karyo_somies <- karyo_list[,rownames(int_somy_matrix)]
  
  print(paste0("Building network for ", nrow(karyo_list), " karyotypes using method: ", method, "..."))
  #create distance matrix between karyotypes
  karyo_dist_matrix <- karyo_dist(karyo_somies)
  
  #now build a minimum spanning tree with this distance matrix
  if(tolower(method) == "msn"){
    network <- pegas::msn(karyo_dist_matrix)
  }else{
    if(tolower(method) == "rmst"){
      network <- pegas::rmst(karyo_dist_matrix)
    }else{
      network <- pegas::mst(karyo_dist_matrix)
    }
  }
  
  #extract the nodes and edge lists
  nodes <- cbind(data.frame(name = karyo_list$karyo_id), karyo_list)
  edges <- data.frame(from = network[,1], to = network[,2], somy_changes = network[,3], alt_link = FALSE)
  
  #getting alternative links (only valid for randomized minimum spanning trees)
  if(!is.null(attr(network, "alter.links"))){
    edges <- rbind(edges, data.frame(from = attr(network, "alter.links")[,1], to = attr(network, "alter.links")[,2], somy_changes = attr(network, "alter.links")[,3], alt_link = TRUE))
  }
  #the mst, rmst, and msn functions create an edge list with the indexes of each node. Here it will replace them by their labels.
  edges$from <- attr(network, "labels")[edges$from]
  edges$to <- attr(network, "labels")[edges$to]
  
  #calculate somy change per chromosome and append it to the edge list
  edges <- cbind(edges, abs(karyo_somies[edges$to,] - karyo_somies[edges$from,]))
  
  #return scDNA object
  scDNAobj$karyotypes$karyonet_nodes <- nodes
  scDNAobj$karyotypes$karyonet_edges <- edges
  scDNAobj$karyotypes$karyonet_dist_matrix <- karyo_dist_matrix
  
  return(invisible(scDNAobj))
}

#plot_karyonet####
#this function will plot the karyotype network created by the karyo_network() function.

plot_karyonet <- function(scDNAobj, color_by = NULL){
  #get required packages
  require(visNetwork)
  require(dplyr)
  require(colorspace)
  require(hues) 
  
  #get needed information 
  nodes <- scDNAobj$karyotypes$karyo_list
  edges <- scDNAobj$karyotypes$karyonet_edges
  chromo_cols <- rownames(scDNAobj$somies$int_somy_matrix)
  sample_name <- unique(scDNAobj$metadata$cells_meta$sample_name)
  
  #formating the nodes list
  nodes <- nodes %>%
    mutate(id = karyo_id,
           value = ifelse(ncells >1, round(proportion*10000), 1),
           value = sqrt(value),
           size = value,
           fixed = size == max(size),
           mass = (sqrt(proportion)*10),
           label = ifelse(ncells == 1, NA, as.character(karyo_id)),
           #label = paste0(gsub("kar", " ", label), " "),
           shape = ifelse(ncells == 1, "dot", "circle"))
  nodes$title <- create_nodes_labels(nodes, chromo_cols)
 
  
  #color the groups
  if(missing(color_by) | is.null(color_by)){
    nodes$group <- "1"
    show_legend <- FALSE
  }else{
    show_legend <- TRUE
    nodes$group <- as.factor(nodes[,color_by])
  }
  #get all groups
  groups <- sort(unique(nodes$group))
  #check which is the group with the largest size
  main_group <- nodes %>%
    group_by(group) %>%
    summarise(size = sum(size)) %>%
    slice_max(size) %>%
    getElement("group")
  #set an individual color for each group
  colors <- iwanthue(length(groups), lmin = 70, lmax = 100)
  names(colors) <- groups
  #make sure that the largest group get the dark grey color
  colors[main_group] <- "#454545"
  #add the color to the nodes and adjust label color as well
  nodes <- nodes %>%
    mutate(color = colors[group],
           font.color = ifelse(color == "#454545", "white", "black"))
  
  #create legend for the nodes
  nodes_legend <- data.frame(group = groups, label = as.character(groups), color = colors[groups]) %>%
    mutate(font.color = ifelse(color == "#454545", "white", "black")) %>%
    arrange(group)
  
  #formating edge list 
  edges <- edges %>%
    mutate(color = ifelse(somy_changes <= 1, "black", "orange"),
           width = ifelse(somy_changes <= 1, 10, 1),
           label = ifelse(somy_changes <= 1, NA, as.character(somy_changes)),
           length = ifelse(somy_changes <= 1, somy_changes, somy_changes*20),
           physics = !alt_link,
           dashes = alt_link,
           color = ifelse(alt_link, "grey", color),
           label = ifelse(alt_link, NA, label),
           width = ifelse(alt_link, 1, width),
           smooth = physics)
           
  
  #adding html popup to edges
  edges$title <- create_edge_labels(edges, nodes, chromo_cols)
  edges <- mutate(edges, title = ifelse(alt_link, NA, title))
  
  #making the title of the plot
  title <- paste0(sample_name, ": ", nrow(karyo_list), " karyotypes linked by a total of ", sum(edges$somy_changes[which(!edges$alt_link)]), " minimum somy changes.")
  
  #making the plot
  network <- visNetwork(nodes, edges, main = title, height = "95vh", width = "100%") %>%
    visNodes(scaling = list(label = list(min = 1, max = 100, drawThreshold = 0))) %>%
    visEdges(font = list(size = 30), scaling = list(label = list(drawThreshold = 0))) %>%
    visIgraphLayout(physics = TRUE, layout = "layout_with_fr") %>%
    visPhysics(stabilization = FALSE, solver = "barnesHut", barnesHut = list(avoidOverlap = 0.5, centralGravity = 0.1, springConstant = 0.1, damping = 0.5), timestep = 0.3) %>%
    visLegend(enabled = show_legend, main = color_by, useGroups = FALSE, addNodes = nodes_legend, zoom = FALSE, stepY = 50) %>%
    visInteraction(tooltipStyle = 'position: fixed;visibility:hidden;padding: 5px;line-height:1; white-space: nowrap;
    font-size:10px;background-color: #f2edbd;') %>%
    visOptions(selectedBy = color_by, highlightNearest = list(enabled = FALSE), nodesIdSelection = TRUE , collapse = FALSE) 
    
  return(network)
  
}

#internal functions for the network plot####
#create_nodes_labels####
#this function will create the HTML formated text to be displayed as the tooltip of the nodes. 
create_nodes_labels <- function(nodes, chromo_cols){
  
  labels <- as.character(nodes$karyo_id)
  karyotypes <- nodes[,chromo_cols]
  
  colors <- heat_col(unique(as.numeric(unlist(karyotypes))))
  
  for(i in c(1:nrow(nodes))){
    karyotype <- as.numeric(karyotypes[i,])
    karyo_id <- paste("<b>", labels[i], "</b><br>", sep = "")
    number_of_cells <- paste(nodes$ncells[i],
                             ifelse(nodes$ncells[i] > 2, "cells", "cell"),
                             paste0("<br>", round(nodes$proportion[i]*100, 2), "% of population"),
                             "<br><br><b>Somies:</b><br>")
    karyo_colors <- colors[as.character(karyotype)]
    chromos <- c(1:length(karyotype))
    chromos[which(as.numeric(chromos) < 10)] <- paste("0", chromos[which(as.numeric(chromos) < 10)], sep = "")
    chromos <-  paste("chr", chromos, ":", sep = "")
    karyotype <- paste("<b><span style=\"color:", karyo_colors, "\">", karyotype, "</span></b>", sep = "") #adding color to each somy value
    karyotype <- paste(chromos, " <b>",karyotype,"</b>", sep = "")
    karyotype <- paste(karyotype, "<br>", sep = "")#adding a line break at the end of each somy value
    karyotype <- paste(karyotype, collapse = "")
    label <- paste(karyo_id, number_of_cells, karyotype,  collapse = "")
    labels[i] <- label
  }
  return(labels)
}


#create_edge_labels####
#this function will create the HTML formated text to be displayed as the tooltip of the edges.
create_edge_labels <- function(edges, karyo_list, chromo_cols){
  from <- edges$from
  to <- edges$to
  
  #comparing each pair ok karyotypes
  labels <- character(nrow(edges))
  
  for(i in c(1:nrow(edges))){
    #converting to vectors
    kar1_id <- edges$from[i]
    kar2_id <- edges$to[i]
    kar1 <- karyo_list[kar1_id, chromo_cols]
    kar2 <- karyo_list[kar2_id, chromo_cols]
    
    #getting which are the different chromosomes
    diff_chromo <- which(kar1 != kar2)
    somy_change <- as.numeric(kar1[diff_chromo]) - as.numeric(kar2[diff_chromo])
    #marking in red the chromosome that changed between both karyotypes (using html code)
    kar1[diff_chromo] <- paste("<b><span style=\"color:red\">", kar1[diff_chromo], "</span></b>", sep = "")
    kar2[diff_chromo] <- paste("<b><span style=\"color:red\">", kar2[diff_chromo], "</span></b>", sep = "")
    #adding a space between somies of kar 1 and kar 2
    label <- paste(kar1, kar2, sep = " <-> ")
    #adding number the total number of somy changes
    #label <- paste()
    #making a vector with chromosomes names
    chromos <- colnames(karyo_list[,chromo_cols])
    #adding chromosomes names to the label
    label <- paste(chromos, label)
    #adding the total number of somy changes between the two karyotypes
    somy_change <- paste(sum(abs(somy_change)), "somy changes<br>")
    label <- c(somy_change, label)
    #adding the name of the karyotypes that the edge links
    karyos <- paste0("<b>", kar1_id, " <-> " , kar2_id, "</b>")
    label <- c(karyos, label)
    #adding a html line brek between each chromosome
    label <- paste0(label, "<br>")
    label <- paste(label, collapse = "")
    #making karyotypes names
    
    labels[i] <- label
  }
  return(labels)
} 


#aggregate_objects####
#this is a function that will combine multiple scDNAobjects into a single one
#it is incomplete for now but should work with objects which were not processed yet (no calc_somy, etc)
#objects which are lost during aggregation: models, 
aggregate_objects <- function(scDNAobj_list){
  
  #libraries
  require(dplyr)
  
  #bins_meta of aggregated object will be the bins_meta of the first object
  bins_meta <- scDNAobj_list[[1]]$metadata$bins_meta
  
  #for the other object, bind them from the list
  agg_count_matrix <- c()
  agg_cells_meta <- c()
  for(i in c(1:length(scDNAobj_list))){
    cells <- gsub("-[0-9]", paste0("-", i), rownames(scDNAobj_list[[i]]$metadata$cells_meta))
    count_matrix <- scDNAobj_list[[i]]$counts$raw_counts
    cells_meta <- scDNAobj_list[[i]]$metadata$cells_meta
    colnames(count_matrix) <- cells
    rownames(cells_meta) <- cells
    cells_meta$original_barcode <- cells_meta$barcode
    cells_meta$barcode <- cells
    agg_count_matrix <- cbind(agg_count_matrix, count_matrix)
    agg_cells_meta <- bind_rows(agg_cells_meta, cells_meta)
  }
  
  #now we update the bins_meta with the data of aggregated count_matrix 
  bins_meta$mean_count <- rowMeans(agg_count_matrix)
  bins_meta$normalized_mean_count <- rowMeans(apply(agg_count_matrix[rownames(bins_meta),], 2, normalize, method = "mean"))
  bins_meta$correction_factor <- NULL
  bins_meta$is_empty <- NULL
  bins_meta$is_outlier <- NULL
  
  scDNAobj <- list(counts = c(), metadata = c())
  scDNAobj$counts$raw_counts <- agg_count_matrix
  scDNAobj$metadata$cells_meta <- agg_cells_meta
  scDNAobj$metadata$bins_meta <- bins_meta
  return(invisible(scDNAobj))
}

#bind_dataframes####
#this is a custom function used to retrieve a specific dataframe from multiple scDNA objects stored in a list
pull_and_bind <- function(list, to_pull){
  require(purrr)
  retrieved_data <- c()
  for(i in c(1:length(list))){
    df <-  flatten(list[[i]])[[to_pull]]
    retrieved_data <- rbind(retrieved_data, df)
  }
   return(retrieved_data)
}

