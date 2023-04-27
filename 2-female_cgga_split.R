# Female CGGA Level_Split Comparison

library(tidyverse)
library(edgeR)
source("process_dge_fun.R")

counts_file <- "data/read_counts.csv"
metadata_file <- "data/metadata.csv"

metadata <- read_csv(metadata_file) |>
  select(CGGA_ID, Level_Split, Gender) |>
  mutate(diag = Level_Split) |>
  filter(!is.na(diag), Gender == "Female")

counts <- read_csv(counts_file) |>
  select(HGNC_Symbol, any_of(metadata$CGGA_ID)) |>
  column_to_rownames("HGNC_Symbol")

groups <- metadata$diag
design_matrix <- model.matrix( ~ 0 + groups)

dge <- DGEList(counts = counts, group = groups, genes = rownames(counts))

keep <- filterByExpr(dge)
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

out_filtered <-
  process_dge(dge_filtered, design_matrix, num_genes = 1000)

out_filtered$complete_table %>%
  rownames_to_column("HGNC_Symbol") %>%
  write_csv("results/female_cgga_split_dge.csv")

saveRDS(out_filtered, file = "data_store/female_cgga_split_dge.rds")
