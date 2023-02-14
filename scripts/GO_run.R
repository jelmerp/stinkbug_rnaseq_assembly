# SET-UP -----------------------------------------------------------------------
## Load packages
library(tidyverse)
library(goseq)

## Source script with functions
source("scripts/GO_fun.R")

## Define input files
DE_age_file <- "results/DE/taurulus_age.txt"
DE_day_file <- "results/DE/taurulus_daylength.txt"
GO_map_file <- "/fs/project/PAS0471/frederico/results/EnTap/taurulus/final_results/final_annotations_no_contam_lvl1_enrich_geneid_go.tsv"
gene_length_file <- "fred_taurulus_genelen.tsv"

## Define output files
outdir <- "results/GO"
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
GO_res_file <- file.path(outdir, "GO_all.txt")


# PREP INPUT -------------------------------------------------------------------
## Read input files
DE_age <- read_tsv(DE_age_file, show_col_types = FALSE)
DE_day <- read_tsv(DE_day_file, show_col_types = FALSE)
GO_map <- read_tsv(GO_map_file, show_col_types = FALSE)
gene_lengths <- read_tsv(gene_length_file,
                         show_col_types = FALSE,
                         col_names = c("trans_id", "length"))

## Combine DE results
DE <- bind_rows(DE_age, DE_day)

## Get mean transcript length for each gene
gene_lengths <- gene_lengths %>%
  mutate(gene_id = sub("(.*)(t\\d+)$", "\\1", trans_id)) %>%
  select(-trans_id) %>% 
  group_by(gene_id) %>%
  summarize(length = mean(length))

## Process the GO map
GO_map <- GO_map %>%
  mutate(gene_id = sub("(.*)(t\\d+)$", "\\1", gene_id)) %>%
  arrange(gene_id) %>%
  distinct() %>%
  as.data.frame()       # goseq fails with a tibble


# RUN GO ANALYSIS --------------------------------------------------------------
GO_res <- map_dfr(.x = unique(DE$contrast),
                  .f = GO_wrap,
                  DE_res = DE, GO_map = GO_map, gene_lens = gene_lengths,
                  p_DE = 0.01, lfc_DE = 2,
                  DE_direction = "both") # One of: 'either', 'up', 'down', 'both'
GO_sig <- GO_res %>% filter(sig == 1) 
write_tsv(GO_sig, GO_res_file)

## Significant categories:
GO_sig %>% pull(description)

## Check virus-related categories:
GO_res[grep("virus", GO_res$description, ignore.case = TRUE), ]

## Check categories that are no longer significant after multiple testing correction:
GO_res %>% filter(p < 0.05, padj > 0.05, ontology == "BP")


# CREATE PLOTS -----------------------------------------------------------------
## Plot for either DE direction:
GO_plot(GO_res = GO_res,
        contrasts = unique(DE$contrast),
        xlabs = c("age", "daylength"),
        title = "taurulus")

## Plot with up- and downregulated genes separately:
GO_plot(GO_res = GO_res,
        facet_var = "contrast", x_var = "DE_direction",
        contrasts = unique(DE$contrast),
        xlabs = c("age", "daylength"),
        title = "taurulus")
