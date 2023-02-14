# SET-UP -----------------------------------------------------------------------
## Load packages
if(!"pacman" %in% installed.packages()) install.packages("pacman")
library(pacman)
packages <- c("tidyverse")
p_load(char = packages, install = TRUE)

## Define input files
annot_file <- "results/entap/final_ed/final_annotations_no_contam_lvl0.tsv"

## Define output files
outdir <- "results/kegg"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
kegg_map_file <- file.path(outdir, "kegg_map.txt")


# Create KEGG map --------------------------------------------------------------
kegg_map <- read.delim(annot_file, sep = "\t", header = TRUE) %>% 
  select(trans_id = Query.Sequence,
         kegg_pathway = EggNOG.KEGG.Terms) %>%
  as_tibble() %>% 
  filter(!is.na(kegg_pathway), kegg_pathway != "") %>% 
  separate(kegg_pathway, sep = ",", fill = "right", 
           into = paste0("k", 1:max(str_count(.$kegg_pathway, ",") + 1))) %>%
  pivot_longer(cols = -trans_id, names_to = NULL, values_to = "kegg_pathway") %>%
  filter(!is.na(kegg_pathway)) %>%
  mutate(kegg_pathway = paste0("map", sprintf("%05d", as.integer(kegg_pathway))),
         gene_id = sub("t\\d+", "", trans_id)) %>% 
  arrange(gene_id) %>%
  select(kegg_pathway, gene_id)

write_tsv(kegg_map, kegg_map_file)
