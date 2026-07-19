# SET-UP -----------------------------------------------------------------------
# Load packages
library(tidyverse)
library(clusterProfiler)

# Source script with functions
source("mcic-scripts/rnaseq/rfuns/enrich_funs.R")

# Define input files
GO_map_file <- "/fs/ess/PAS0471/frederico/results/EnTap/taurulus/final_results/final_annotations_no_contam_lvl1_enrich_geneid_go.tsv"

# Define output files
outdir <- "results/GO"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
GO_mapfile <- file.path(outdir, "GO_map.tsv")

# Read the GO file
GO_init <- read_tsv(GO_map_file, col_types = "cc") |>
  dplyr::rename(gene = gene_id, term = go_term)


# MAIN -------------------------------------------------------------------------
# Get GO info: ontology and description
GO_info <- get_GO_info()

# Get GO info: hierarchy levels
GO_levels <- get_GO_levels()
lvl2_terms <- GO_levels |> filter(GO_lvl < 3) |> pull(GO_ID)

# Extract all valid GO terms mapped to the fly genome
library(org.Dm.eg.db)
insect_whitelist <- AnnotationDbi::keys(org.Dm.eg.db, keytype = "GO")

# Create the GO map
GO_map <- GO_init |>
  mutate(gene = sub("(.*)(t\\d+)$", "\\1", gene)) |>
  distinct() |>
  left_join(GO_info, by = "term") |>
  dplyr::select(term, gene, description, ontology) |>
  # Remove terms with no description
  filter(!is.na(description)) |> 
  # Remove terms not annotated to fruitfly
  filter(term %in% insect_whitelist) |>
  # Remove levels 1 and 2:
  filter(!term %in% lvl2_terms) |>
  arrange(term)

# Write the output file
write_tsv(GO_map, GO_mapfile)
