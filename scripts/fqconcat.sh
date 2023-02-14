#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=5
#SBATCH --out=slurm-fqconcat-%j.out


# SETUP ------------------------------------------------------------------------
## Bash strict settings
set -euo pipefail

## Command-line args
indir1=""
indir2=""
outdir=""

indir1=$1
indir2=$2
outdir=$3

## Hardcoded variables
SBATCH_ACCOUNT=PAS0471

## Report
echo "## Starting script fqconcat.sh"
date
echo
echo "## Input dir 1:     $indir1"
echo "## Input dir 2:     $indir2"
echo "## Output dir:      $outdir"
echo -e "-------------------\n"

## Input checks
[[ "$indir1" = "$outdir" ]] && echo "## ERROR: Input dir 1 should not be the same as the output dir" && exit 1
[[ "$indir2" = "$outdir" ]] && echo "## ERROR: Input dir 2 should not be the same as the output dir" && exit 1
[[ ! -d "$indir1" ]] && echo "## ERROR: Input dir $indir1 does not exist" && exit 1
[[ ! -d "$indir2" ]] && echo "## ERROR: Input dir $indir2 does not exist" && exit 1

## Create output dir, if necessary
mkdir -p "$outdir"

## Make input dir read-only
echo "## Making input files read-only..."
chmod a-w "$indir1"/*fastq.gz "$indir2"/*fastq.gz


# CONCAT FASTQ FILES -----------------------------------------------------------
echo -e "\n## Looping through files in $indir1..."
for fq1 in "$indir1"/*fastq.gz; do
    ## Determine sample ID, read direction etc
    fq1_base=$(basename "$fq1")
    sample_id=$(echo "$fq1_base" | sed -E 's/(.*)NGS312_.*/\1/')
    sample_id2=$(echo "$sample_id" | sed -E 's/-/_/')
    read_dir=$(echo "$fq1_base" | sed -E 's/.*_L001_(R[12])_001.fastq.gz/\1/')
    
    ## Determine name of secondary file for same sample and read direction
    fq2=$(find "$indir2"/ -name "${sample_id2}*${read_dir}*")
    
    ## Define name for concatenated output file
    fq_concat=$outdir/"$sample_id"_"$read_dir"_001.fastq.gz
    
    ## Report
    echo
    echo "## Sample ID:                $sample_id"
    echo "## FASTQ file 1:             $fq1"
    echo "## FASTQ file 2:             $fq2"
    echo "## Concatenated FASTQ file:  $fq_concat"

    ## Submit a job to concatenate the two files
    cmd="zcat $fq1 $fq2 | gzip > $fq_concat; du -h $fq1 $fq2 $fq_concat; chmod a-w $fq_concat"
    sbatch --account="$SBATCH_ACCOUNT" --wrap="$cmd"
    echo "------------------"
done

## Report
echo -e "\n## Done with script fqconcat.sh"
date