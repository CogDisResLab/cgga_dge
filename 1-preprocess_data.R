# Process the raw data into usable form.

library(tidyverse)

waystation <- tempdir()

raw_files <- list.files("raw")

unzipped <- raw_files |>
  map( ~ unzip(file.path("raw", .x), exdir = waystation)) |>
  map_chr( ~ keep(.x,  ~ str_detect(.x, "MACOSX", negate = TRUE))) |>
  set_names(basename)

metadata_325 <- read_tsv(unzipped[1]) |>
  select(CGGA_ID, Age, Gender, PRS_type, Histology, Grade) |>
  filter(!is.na(Grade), !is.na(Gender)) |>
  mutate(
    Level_Split = case_when(Grade == "WHO II" ~ "Low",
                            .default = "High"),
    Level_Compare = case_when(Grade == "WHO II" ~ "Low",
                              Grade == "WHO IV" ~ "High",
                              .default = NA)
  ) |>
  write_csv(file.path("data", "metadata_325.csv"))

read_counts_325 <- read_tsv(unzipped[2]) |>
  select(HGNC_Symbol = gene_name, any_of(metadata_325$CGGA_ID)) |>
  write_csv(file.path("data", "read_counts_325.csv"))

metadata_693 <- read_tsv(unzipped[3]) |>
  select(CGGA_ID, Age, Gender, PRS_type, Histology, Grade) |>
  filter(!is.na(Grade), !is.na(Gender)) |>
  mutate(
    Level_Split = case_when(Grade == "WHO II" ~ "Low",
                            .default = "High"),
    Level_Compare = case_when(Grade == "WHO II" ~ "Low",
                              Grade == "WHO IV" ~ "High",
                              .default = NA)
  ) |>
  write_csv(file.path("data", "metadata_693.csv"))

read_counts_693 <- read_tsv(unzipped[4]) |>
  select(HGNC_Symbol = gene_name, any_of(metadata_693$CGGA_ID)) |>
  write_csv(file.path("data", "read_counts_693.csv"))

metadata <- bind_rows(metadata_325, metadata_693) |>
  write_csv(file.path("data", "metadata.csv"))

read_counts <- read_counts_325 |>
  inner_join(read_counts_693) |>
  write_csv(file.path("data", "read_counts.csv"))

unlink(waystation, recursive = TRUE)
