## Load packages 
if (!"pacman" %in% installed.packages()) install.packages("pacman")
if (!"BiocManager" %in% installed.packages()) install.packages("BiocManager")
if (!"tximport" %in% installed.packages()) BiocManager::install("tximport")
if (!"rhdf5" %in% installed.packages()) BiocManager::install("rhdf5")
if (!"DESeq2" %in% installed.packages()) BiocManager::install("DESeq2")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
packages <- c("tidyverse", "here", "PCAtools", "tximport", "rhdf5", "DESeq2",
              "apeglm", "BiocParallel", "ashr","vsn", "pheatmap")
pacman::p_load(char = packages, install = TRUE)

## Define input files
annot_file <- "/fs/project/PAS0471/frederico/results/EnTap/taurulus/final_ed/final_annotations_no_contam_lvl0.tsv"
kallisto_dir <- "/fs/project/PAS0471/frederico/results/kallisto/taurulus"
dds_in <- here(kallisto_dir, "deseq_object.rds")

## Define output files
res_age_file <- "results/DE/taurulus_age.txt"
res_day_file <- "results/DE/taurulus_daylength.txt"

## Load annotation
annot <- read_tsv(annot_file, show_col_types = FALSE) %>%
  mutate(geneID = sub("t1$", "", `Query Sequence`)) %>%
  select(geneID, Description, Species, Frame, subj = `Subject Sequence`)

## Read DESeq object
dds <- readRDS(dds_in)


# ANALYSIS BY DAYLENGTH --------------------------------------------------------
dds_day <- dds
design(dds_day) <- formula(~age + daylength)
dds_day <- DESeq(dds_day)

res_day <- results(dds_day, name = "daylength_short_vs_long") %>%
  as.data.frame() %>%
  arrange(padj) %>%
  rownames_to_column("gene_id") %>%
  mutate(contrast = "daylength")

write_tsv(res_day, res_day_file)


# ANALYSIS BY AGE --------------------------------------------------------------
dds_age <- dds
design(dds_age) <- formula(~daylength + age)
dds_age <- DESeq(dds_age)

res_age <- results(dds_age, name = "age_nymph_vs_adult") %>%
  as.data.frame() %>%
  arrange(padj) %>%
  rownames_to_column("gene_id") %>%
  mutate(contrast = "age")

write_tsv(res_age, res_age_file)
