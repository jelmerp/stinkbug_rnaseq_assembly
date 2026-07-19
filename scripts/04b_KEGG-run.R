# SET-UP -----------------------------------------------------------------------
# Load packages and function-files
library(tidyverse)
library(clusterProfiler)
library(here)
source(here("mcic-scripts/rnaseq/rfuns/enrich_funs.R"))

# Define input files
DE_file <- here("results/DE/DE_taurulus_adults.tsv")
kegg_map_file <- here("results/KEGG/KEGG_map.tsv")

# Define output files
kegg_outfile <- here("results/KEGG/KEGG_res.tsv")

# Read input files
DE <- read_tsv(DE_file, show_col_types = FALSE)
kegg_map <- read_tsv(kegg_map_file, show_col_types = FALSE)


# KEGG ENRICHMENT ANALYSIS -----------------------------------------------------
# For upregulated DEGs
up <- run_ora(
  df = DE,
  term_map = kegg_map, 
  contrast = "long_short",
  DE_direction = "up",
  simplify_terms = TRUE,
  return_df = TRUE
  )

# For downregulated DEGs
down <- run_ora(
  df = DE,
  term_map = kegg_map, 
  contrast = "long_short",
  DE_direction = "down",
  simplify_terms = TRUE,
  return_df = TRUE
)

# Combine the two tables, keep only significant ones
kegg <- bind_rows(up, down) |> filter(sig == TRUE)

# Write output to file
write_tsv(kegg, kegg_outfile)
