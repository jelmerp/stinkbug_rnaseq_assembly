# SET-UP -----------------------------------------------------------------------
## Install/load packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse",       # Misc. data manipulation and plotting
              "clusterProfiler") # Enrichment testing
pacman::p_load(char = packages)

## Source script with functions
source("scripts/kegg_fun.R")

## Input files
DE_age_file <- "results/DE/taurulus_age.txt"
DE_day_file <- "results/DE/taurulus_daylength.txt"
kegg_dir <- "results/kegg"
kegg_map_file <- file.path(kegg_dir, "kegg_map.txt")
## TODO - GET PATHWAY DESCRIPTIONS

## Output file
kegg_res_file <- file.path(kegg_dir, "kegg_all.txt")

## Read input files
DE_age <- read_tsv(DE_age_file, show_col_types = FALSE)
DE_day <- read_tsv(DE_day_file, show_col_types = FALSE)
kegg_map <- read_tsv(kegg_map_file, show_col_types = FALSE)
#kegg_descr <- read_tsv(kegg_descr_file, show_col_types = FALSE)

## Combine DE results
DE <- bind_rows(DE_age, DE_day)


# KEGG ENRICHMENT TEST ---------------------------------------------------------
## Get all pairwise contrasts
contrasts <- unique(DE$contrast)

## Run enrichment test for all pairwise contrasts
kegg_both <- map_dfr(.x = contrasts, .f = run_kegg,
                     direction = "both", DE_res = DE, kegg_map = kegg_map)

kegg_up <- map_dfr(.x = contrasts, .f = run_kegg, "up", DE, kegg_map)
kegg_down <- map_dfr(.x = contrasts, .f = run_kegg, "down", DE, kegg_map)

kegg_res <- bind_rows(kegg_both, kegg_up, kegg_down) %>%
  select(contrast,
         direction,
         padj = p.adjust,
         count = Count,
         category = ID,
         description = Description,
         gene_ids = geneID)

## Write output file
write_tsv(kegg_res, outfile_kegg_res)
