# SET-UP -----------------------------------------------------------------------
## Load packages
library(tidyverse)
library(here)

## Load script with functions
source("scripts/DE_fun.R")

## Define input files
DE_age_file <- "results/DE/taurulus_age.txt"
kallisto_dir <- "/fs/project/PAS0471/frederico/results/kallisto/taurulus"
dds_in <- here(kallisto_dir, "deseq_object.rds")
annot_file <- "/fs/project/PAS0471/frederico/results/EnTap/taurulus/final_ed/final_annotations_no_contam_lvl0.tsv"

## Read DESeq object
dds <- readRDS(dds_in)

## Load annotation
annot <- read_tsv(annot_file, show_col_types = FALSE) %>%
  mutate(geneID = sub("t1$", "", `Query Sequence`)) %>%
  select(geneID, Description, Species, Frame, subj = `Subject Sequence`)

## To plot 1 gene:
plot_counts(dds = dds[, dds$age == "adult"],
            geneID = "NonamEVm004236",
            annot = annot)
