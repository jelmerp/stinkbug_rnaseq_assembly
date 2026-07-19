# SET-UP -----------------------------------------------------------------------
# Load packages and function-files
library(tidyverse)
library(clusterProfiler)
library(here)
source(here("mcic-scripts/rnaseq/rfuns/enrich_funs.R"))

# Define input files
#DE_day_file <- "/fs/ess/PAS0471/frederico/results/DE/taurulus_dl.txt" # File with DE results for both ages, no longer appropriate
DE_file <- here("results/DE/DE_taurulus_adults.tsv")
GO_map_file <- here("results/GO/GO_map.tsv")

# Define output files
outdir <- here("results/GO")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
GO_outfile <- here(outdir, "GO_dl.tsv")

# Read input files
DE <- read_tsv(DE_file, show_col_types = FALSE)
GO_map <- read_tsv(GO_map_file, show_col_types = FALSE)

# Check number of up- and downregulated genes
DE |> filter(padj < 0.05) |> count(lfc > 0)  # 2535, 2685 (Both ages (earlier results): 1460, 1220)


# GO ENRICHMENT ANALYSIS -------------------------------------------------------
# For upregulated DEGs
up <- run_ora(
  df = DE,
  term_map = GO_map, 
  contrast = "long_short",
  DE_direction = "up",
  simplify_terms = TRUE,
  return_df = TRUE
  )

# For downregulated DEGs
down <- run_ora(
  df = DE,
  term_map = GO_map, 
  contrast = "long_short",
  DE_direction = "down",
  simplify_terms = TRUE,
  return_df = TRUE
)

# Combine the two tables, keep only significant ones
GO <- bind_rows(up, down) |> filter(sig == TRUE)

# Write output to file
write_tsv(GO, GO_outfile)
