# 2026-07-18, Jelmer Poelstra, R 4.5.2 at OSC-cardinal using Positron
# (Made using /fs/ess/PAS0471/frederico/scripts/DESeq_taurulus_dl.R as template)


# SETUP ------------------------------------------------------------------------
# Load packages and function-files
library(tidyverse)
library(DESeq2)
library(here)
source(here("mcic-scripts/rnaseq/rfuns/DE_funs.R"))

# Define input files
dds_file <- "/fs/ess/PAS0471/frederico/results/kallisto/taurulus/deseq_object.rds"

# Define output files
outdir <- here("results/DE")
DE_file <- here(outdir, "DE_taurulus_adults.tsv")


# MAIN -------------------------------------------------------------------------
# Read input files
dds <- readRDS(dds_file)

# DE analysis for adults only, by daylength
dds <- dds[, dds$age == "adult"]
design(dds) <- formula(~ daylength)
dds <- DESeq(dds)

# Extract the results
DE <- extract_DE(dds, fct = "daylength", comp = c("long", "short"))

#DE <- results(dds) %>%
#  as.data.frame() %>%
#  arrange(padj) %>%
#  rownames_to_column("gene_id")


# WRITE OUTPUTS ----------------------------------------------------------------
write_tsv(DE, DE_file)
