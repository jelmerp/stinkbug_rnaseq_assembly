## Packages
library(tidyverse)

## Input files
fnames <- list.files("results/kallisto", pattern = "fastq.gz")

## Output files
outdir <- "data/meta"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## Infer sample names and factors
sampleID <- sub("_R1.fastq.gz", "", fnames)
daylength <- ifelse(grepl("1-", sampleID), "long", "short")
age <- ifelse(grepl("NST", sampleID), "nymph", "adult")

## Create dataframe and write to file
meta <- data.frame(sampleID, daylength, age)
write_tsv(meta, file.path(outdir, "samples.txt"))

# "NST" = nymph, "ST" = adult
# "1-" = long day, "2-" = short day
