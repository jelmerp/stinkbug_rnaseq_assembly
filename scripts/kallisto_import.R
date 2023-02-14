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
kallisto_dir <- here("results/kallisto/")
meta_file <- here("data/meta/samples.txt")
transcript_file <- here("results/evigene/okayset/transcriptIDs.txt")

## Output file
dds_out <- here(kallisto_dir, "deseq_object.rds")

## Settings
sampleid_column <- "sampleID"  # Column in metadata with desired sample names


# PROCESS METADATA AND KALLISTO COUNTS -----------------------------------------
## Read metadata
meta <- read.delim(meta_file) %>%
  arrange(.data[[sampleid_column]]) %>%
  mutate(daylength = factor(daylength, levels = c("long", "short")),
         age = factor(age))
rownames(meta) <- meta[[sampleid_column]]

## Read transcript-to-gene map
transcript <- readLines(transcript_file)
geneID <- sub("t\\d+$", "", transcript)
#transcriptIDs <- sub(".*(t\\d+$)", "\\1", transcript)
transmap <- data.frame(TXNAME = transcript, GENEID = geneID) %>%
  arrange(GENEID, TXNAME)

## Import Kallisto transcript counts --
## create gene-level count estimates normalized by library size and transcript length
## See https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
kallisto_dirs <- list.files(kallisto_dir, pattern = "fastq.gz")
kallisto_files <- here(kallisto_dir, kallisto_dirs, "abundance.h5")
names(kallisto_files) <- sub("_R1.fastq.gz", "", kallisto_dirs)

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
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = meta,
                                design = ~1)

## Save DESeq object
saveRDS(dds, dds_out)

# FOR SEPARATE DE SCRIPT -----------------------
## To load:
# dds <- readRDS(dds_out)
## To add a design:
# dds_bothfactors <- dds
# design(dds_bothfactors) <- formula(~ age + daylength) # Without interaction
# design(dds_bothfactors) <- formula(~ age * daylength) # With interaction
# ## To run the DE analysis
# dds_bothfactors <- DESeq(dds_bothfactors)
# ## To extract results
# results(dds_bothfactors)
# -----------------------------------------------

# WRAP UP ----------------------------------------------------------------------
## List output
message("\n## Listing output file:")
system(paste("ls -lh", dds_out))

message("\n## Done with script tximport.R")
Sys.time()
message()
