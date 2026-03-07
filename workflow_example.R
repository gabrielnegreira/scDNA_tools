#workflow example
#this is a simple script with an example of the workflow of scDNA

#libraries####
library(tidyverse)

#get the files which we will analyze.
files <- list.files(path = "inputs", pattern = "B12|G3|A9|D6|E4|G6|G8|H12", full.names = TRUE)
#files <- files[c(2:6)]
#files <- files[7]
#files <- "inputs/BPK282_cnv_data.h5"

#files <- "inputs/BPK081_cnv_data.h5"
#files <- "inputs/BPK081_cnv_data.h5"

#make a list of scDNA objects 
scDNAobj_list <- lapply(files, build_scDNAobj)
names(scDNAobj_list) <- basename(files)

scDNAobj_list <- scDNAobj_list %>%
  lapply(calc_somy) %>%
  lapply(tag_outlier_cells, ICV_loess_span = 0.75) %>%
  lapply(summarize_karyotypes) %>%
  lapply(karyo_network, show_plot = FALSE)


write_rds(scDNAobj_list, file = "scDNAobj_list.RDS")

plot_karyotypes(scDNAobj_list[[1]])