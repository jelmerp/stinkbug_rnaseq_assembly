# SET-UP -----------------------------------------------------------------------
# Load packages
library(tidyverse)
library(here)

# Define input files
annot_file <- "/fs/ess/PAS0471/frederico/results/EnTap/taurulus/final_results/final_annotations_no_contam_lvl0.tsv"

# Define output files
outdir <- here("results/KEGG")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
kegg_map_file <- here(outdir, "KEGG_map.tsv")


# MAIN -------------------------------------------------------------------------
# Get the KEGG pathway descriptions
kegg_url <- "https://rest.kegg.jp/list/pathway"
kegg_info <- read.delim(kegg_url, header = FALSE, col.names = c("term", "description"))

# Prep the lookup table
kegg_map <- read.delim(annot_file, sep = "\t", header = TRUE) |> 
  select(trans_id = Query.Sequence, term = EggNOG.KEGG.Terms) |>
  as_tibble() |> 
  filter(!is.na(term), term != "") |> 
  separate_longer_delim(term, delim = ",") |>
  filter(!is.na(term)) |>
  mutate(
    term = paste0("map", sprintf("%05d", as.integer(term))),
    gene = sub("t\\d+", "", trans_id)
  ) |>
  # Add pathway descriptions
  left_join(kegg_info, by = "term") |>
  # Remove human disease pathways (map05) and drug development pathways (map07)
  filter(!grepl("map05|map07", term)) |>
  arrange(term) |>
  select(term, gene, description)

# Write the output file
write_tsv(kegg_map, kegg_map_file)
