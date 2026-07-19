#!/bin/bash

# SOFTWARE ---------------------------------------------------------------------
mamba create -y -p /fs/project/PAS0471/jelmer/conda/soapdenovo-trans-1.0.4 -c bioconda soapdenovo-trans=1.04
mamba create -y -p /fs/project/PAS0471/jelmer/conda/transabyss-2.0.1 -c bioconda transabyss=2.0.1

## Evigene
wget https://sourceforge.net/projects/evidentialgene/files/evigene22may07.tar
mamba create -y -p /fs/project/PAS0471/jelmer/conda/evigene -c bioconda exonerate=2.4.0 cd-hit=4.8.1 blast=2.12.0
cp evigene/scripts/* /fs/project/PAS0471/jelmer/conda/evigene/bin

# ENTAP ------------------------------------------------------------------------
## Install EnTAP and dependencies
mamba create -p /fs/project/PAS0471/jelmer/conda/entap-0.10.8 -c conda-forge matplotlib=3.5.2 
mamba install -c bioconda interproscan diamond=0.9.10 transdecoder=5.3.0 rsem=1.3.3

mkdir -p software/entap && cd software/entap  
wget https://github.com/harta55/EnTAP/archive/refs/tags/v0.10.8-beta.tar.gz # https://github.com/harta55/EnTAP/tags
tar -xzvf v0.10.8-beta.tar.gz
cd EnTAP-0.10.8-beta
module load cmake/3.17.2
cmake CMakeLists.txt -DCMAKE_INSTALL_PREFIX=/fs/project/PAS0471/jelmer/conda/entap-0.10.8
make install
cp src/entap_graphing.py /fs/project/PAS0471/jelmer/conda/entap-0.10.8/bin

mkdir -p results/entap/
cp software/entap/EnTAP-0.10.8-beta/entap_config.ini results/entap/

mkdir -p software/genemark-st && cd software/genemark-st
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_UDnOI/gmst_linux_64.tar.gz
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_UDnOI/gm_key_64.gz
cp * /fs/project/PAS0471/jelmer/conda/entap-0.10.8/bin/

## Get NCBI NR database
cd /fs/project/PAS0471/jelmer/refdata/nr
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz
tar -xvzf nr.*
cat nr.* > nr_database.fasta

## Get NCBI RefSeq protein db
cd /fs/project/PAS0471/jelmer/refdata/refseq
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*protein*faa.gz
cat complete* | gunzip -c > refseq_complete.faa

## Get Uniprot/Swissprot db
cd /fs/project/PAS0471/jelmer/refdata/uniprot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# ORGANIZE DATA FILES ----------------------------------------------------------
## Dirs
dir_fq_org1=data/original/NGS312_AlbertoFredHickmann_RNAseq-314740426
dir_fq_org2=data/original/COMP_NGS312_AlbertoFredHickmann_RNAseq-314740426
dir_fq1=data/fastq/run1 && mkdir -p "$dir_fq1"
dir_fq2=data/fastq/run2 && mkdir -p "$dir_fq2"
dir_fq_concat=data/fastq/concat

## Unzip files in one of the original FASTQ dirs
for zip in "$dir_fq_org2"/*zip; do
    unzip "$zip" -d "$dir_fq_org2"
done

## Copy FASTQ files into a single dir
find "$dir_fq_org1" -name "*fastq.gz" -exec cp -v {} "$dir_fq1" \;
find "$dir_fq_org2" -name "*fastq.gz" -exec cp -v {} "$dir_fq2" \;

## Concatenate files from same sample
sbatch scripts/fqconcat.sh "$dir_fq1" "$dir_fq2" "$dir_fq_concat"

## Download virus genome
# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1410467
conda activate /users/PAS0471/jelmer/miniconda3/envs/blast-env
ref_virus=data/ref/Hhvirus/Hhvirus_NC_022611.1.fa
esearch -db nuccore -query "NC_022611.1" | efetch -format fasta > "$ref_virus"
#grep -v ">" data/ref/Hhvirus_NC_022611.1.fa | tr -d "\n" | wc -c # 9271 bp

# Download heros genome
mkdir -p data/ref/heros
wget -P data/ref/heros https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/667/255/GCA_003667255.1_E_heros_v1.0_sep_2018/GCA_003667255.1_E_heros_v1.0_sep_2018_genomic.fna.gz

## Build Kraken2 db with fungi
kraken_db=/fs/project/PAS0471/jelmer/refdata/kraken/std-plus-fungi
sbatch mcic-scripts/metagenomics/kraken-build-custom-db.sh -d "$kraken_db"

## Subset FASTQ files for testruns
mcic-scripts/misc/fqsub_dir.sh -i "$dir_fq_concat" -o data/fastq/concat/subset
