#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --output=slurm-tximport-%j.out


# SET-UP -----------------------------------------------------------------------
message("\n## Starting script tximport.R")
Sys.time()
message()

## Load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
if (!"BiocManager" %in% installed.packages()) install.packages("BiocManager")
if (!"tximport" %in% installed.packages()) BiocManager::install("tximport")
if (!"rhdf5" %in% installed.packages()) BiocManager::install("rhdf5")
if (!"DESeq2" %in% installed.packages()) BiocManager::install("DESeq2")
packages <- c("tidyverse", "here", "tximport", "rhdf5", "DESeq2")
pacman::p_load(char = packages, install = TRUE)

## Input dirs and files
kallisto_dir <- here("results/kallisto/beetle")
meta_file <- here("metadata/samples.txt")
transmap_file <- here("results/trinotate/GENE_TRANS_MAP")

## Output file
dds_out <- here(kallisto_dir, "deseq_object.rds")

## Settings
dirname_column <- "id_long"    # Column in metadata with Kallisto dirnames
sampleid_column <- "id_short"  # Column in metadata with desired sample names


# PROCESS METADATA AND KALLISTO COUNTS -----------------------------------------
## Read metadata
meta <- read.delim(meta_file) %>% arrange(.data[[sampleid_column]])
rownames(meta) <- meta[[sampleid_column]]

## Read transcript-to-gene map
transmap <- read_tsv(transmap_file,
                     col_names = c("GENEID", "TXNAME"),
                     col_types = "cc") %>%
  select(TXNAME, GENEID)

## Import Kallisto transcript counts --
## create gene-level count estimates normalized by library size and transcript length
## See https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
kallisto_files <- here(kallisto_dir, meta[[dirname_column]], "abundance.h5")
names(kallisto_files) <- meta[[sampleid_column]]

txi <- tximport(kallisto_files,
                type = "kallisto",
                tx2gene = transmap,
                countsFromAbundance = "lengthScaledTPM")


# CREATE DESEQ OBJECT ----------------------------------------------------------
## Check that sample names are the same and samples are in same order
stopifnot(all(rownames(meta) == colnames(txi$counts)))
message("\n## Sample names:")
print(rownames(meta))
message("\n## Dimensions of count matrix:")
dim(txi$counts)

## Create DESeq object
dds <- DESeqDataSetFromTximport(txi, meta, ~1)

## Save DESeq object
saveRDS(dds, dds_out)


# WRAP UP ----------------------------------------------------------------------
## List output
message("\n## Listing output file:")
system(paste("ls -lh", dds_out))

message("\n## Done with script tximport.R")
Sys.time()
message()
