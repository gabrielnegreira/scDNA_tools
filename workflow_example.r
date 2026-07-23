#load libraries and functions####
#get the functions
source("scDNA_functions.R")

#get the data####
##get the count matrix####
count_matrix <- as.matrix(read.csv("count_matrix_example.csv", row.names = 1))
#remove kDNA from count matrix 

##make a dummy cells_metadata####
cells_meta <- data.frame(barcode = colnames(count_matrix))

#for now cells_meta is required to have a `cell_id` column with a numeric identifier of that cells. 
# This was based on the cellranger-dna metadata.
#I think this is actually not necessary so I might remove this requirement in the future. 
#But for now, we create one:
cells_meta$cell_id <- 1:nrow(cells_meta)

#the cells metadata must have rownames matching the colnames of the count matrix!
#this is important because plotting functions use the row names of the cells meta to match metadata with count values in the count matrix!
rownames(cells_meta) <- cells_meta$barcode

#confirm they match
all(rownames(cells_meta) == colnames(count_matrix))

#let's inspect cells metadata:
head(cells_meta)

##nmake a dummy bins metadata####
##If we don't have information about the bins gc content, mappability, etc, we can create dummy values
bins_meta <- data.frame(
  start = NA,
  end = NA,
  bin = rownames(count_matrix),
  chromosome = gsub("_.*", "", rownames(count_matrix)), #this gets the chromosome from the bin name
  gc_content = rnorm(n = nrow(count_matrix), mean = 0.5, sd = 0.2),
  mappability = 1,
  is_mappable = TRUE
)

#we also need to set bins_meta rownames to the bins names in the count matrix
##OBS: currently there will be no warning if this is not met. So I should fix it later.
rownames(bins_meta) <- bins_meta$bin
#confirm they match
all(rownames(cells_meta) == colnames(count_matrix))

#let's inspect bins_meta
head(bins_meta)

#build the scDNA obj####
obj <- build_scDNAobj(
  count_matrix = count_matrix, 
  cells_meta = cells_meta,
  bins_meta = bins_meta,
  mean_read_length = 400, #the average read length of the library. We can get that from the TapeStation profile. Only used to calculate the effective depth of coverage.
  genome_size = 33000000 #also only used for calculating effective depth of coverage. Might remove these two in the future.
)

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
obj <- calc_somy(obj, int_method = "round") #for gaussian mixuture models set it to "GMM" (but it might get stuck in a loop sometimes. Need to fix it)
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
