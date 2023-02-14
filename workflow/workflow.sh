# SETUP ------------------------------------------------------------------------
# Settings
minq=5               # Min. qual for TrimGalore trimming - follows https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
minlen=36            # Min. read length for TrimGalore trimming - follows https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html

# Files and dirs
kraken_db=/fs/project/PAS0471/jelmer/refdata/kraken/std-plus-fungi


# QC AND ASSEMBLY PRE-PROCESSING -----------------------------------------------
# Run FastQC
for fq in data/fastq/concat/*fastq.gz; do
    sbatch mcic-scripts/qc/fastqc.sh -i "$fq" -o results/fastqc
done
sbatch mcic-scripts/qc/multiqc.sh -i results/fastqc -o results/multiqc

# Run rcorrector
sbatch mcic-scripts/assembly_trans/rcorrector.sh -i data/fastq/concat -o results/rcorrector

# Remove reads with unfixable errors, according to rcorrector
for R1 in results/rcorrector/*R1*.cor.fq.gz; do
    sbatch mcic-scripts/assembly_trans/rcorrfilter.sh -i "$R1" -o results/rcorrfilter
done

# Run TrimGalore
for R1 in results/rcorrfilter/*R1*fastq.gz; do
    sbatch mcic-scripts/trim/trimgalore.sh -i "$R1" -o results/trimgalore/rcorr \
        -O results/fastqc/trimmed_rcorr -q "$minq" -l "$minlen"
done

# Run sortmerna
repo_dir=results/sortmerna/sortmerna_repo && mkdir -p "$repo_dir"
git clone https://github.com/biocore/sortmerna "$repo_dir"
for R1 in results/trimgalore/rcorr/*_R1_*fastq.gz; do
    sbatch mcic-scripts/rnaseq/sortmerna.sh  -i "$R1" -o results/sortmerna -r "$repo_dir"
done


# CONTAMINANT DETECTION & VIRUS ASSEMBLY --------------------------------------- 
# Run Kraken2 on db with fungi
for R1 in results/sortmerna/unmapped/*_R1_*fastq.gz; do 
    sbatch mcic-scripts/metagenomics/kraken-run.sh -i "$R1" -o results/kraken -d "$kraken_db" -W -w -c 0.5
done
#grep "Halyomorpha halys virus" "$dir_krak"/*report.txt

# Make Krona plots
for kraken_out in results/kraken/*_main.txt; do
    sampleID=$(basename "$kraken_out" _main.txt)
    sbatch mcic-scripts/metagenomics/krona.sh "$kraken_out" results/krona/"$sampleID".html
done

# Extract viral-assigned Kraken2 reads
taxids="699189 1111709 1410467" # Iflaviridae, unclassified Iflaviridae, Halyomorpha halys virus
krakreads_dir=results/kraken/conf0/"${taxids// /_}"_reads
for kraken_main in results/kraken/conf0/*_main.txt; do 
    R1_name=$(basename "$kraken_main" _main.txt | sed 's/_001/_R1_001/')
    R1_in=results/sortmerna/unmapped/"$R1_name".fastq.gz
    ls "$kraken_main" "$R1_in"
    sbatch mcic-scripts/metagenomics/extract_kraken_reads.sh -i "$R1_in" -I "$kraken_main" -t "$taxids" -o "$krakreads_dir"
done

# Assemble virus genome from Kraken reads - heros
find "$krakreads_dir" -name "*ST[12]-*_R1*fastq.gz" -size +1k > config/spades/heros_krakenreads_conf0_R1s.txt
bash mcic-scripts/assembly/spades_R1files2yaml.sh config/spades/heros_krakenreads_conf0_R1s.txt config/spades/heros_krakenreads_conf0.yaml
sbatch -t60 --mem=20G -c4 mcic-scripts/assembly/spades.sh \
    -y config/spades/heros_krakenreads_conf0.yaml -o results/vir_assm/kraken0_heros -m rnaviral

# Map to viral genome
viral_asm=results/vir_assm/kraken0_heros/contigs.fasta
sbatch mcic-scripts/rnaseq/star_index.sh -i "$viral_asm" -o results/viral_asm/heros/star_index

for R1 in results/sortmerna/unmapped/*ST*_R1_*.gz; do
    sbatch mcic-scripts/rnaseq/star_align.sh -i "$R1" \
        -o results/map2viral/heros/ -r results/viral_asm/heros/star_index
done

# Run kraken again on unmapped reads
kraken_db=/fs/project/PAS0471/jelmer/refdata/kraken/std-plus-fungi
for R1 in results/map2virus/heros/unmapped/*_R1.fastq.gz; do
    sbatch mcic-scripts/metagenomics/kraken-run.sh -i "$R1" -o results/kraken/heros_postviral/ -d "$kraken_db" -W -w -c 0.5 -q 0
done

# Read normalization for SOAPdenovo and SPAdes
for R1 in "$fqdir_nonorm"/*R1*fastq.gz; do
    sbatch mcic-scripts/assembly_trans/orna.sh -i "$R1" -o "$fqdir_norm"
done

# ASSEMBLE STINKBUG TRANSCRIPTOME ----------------------------------------------
fqdir_nonorm=results/kraken/heros_postviral/unclassified
fqdir_norm=results/orna/heros

# Run Trinity - de novo
sbatch mcic-scripts/assembly_trans/trinity.sh -i "$fqdir_nonorm" -o results/trinity/heros/trinity_denovo

# Run Trinity - reference-guided for heros
ref=data/ref/heros/GCA_003667255.1_E_heros_v1.0_sep_2018_genomic.fna.gz
sbatch mcic-scripts/rnaseq/star_index.sh -i "$ref" -o data/ref/heros/star_index
for R1 in results/kraken/heros_postviral/unclassified/*_R1*gz; do
    sbatch mcic-scripts/rnaseq/star_align.sh -i "$R1" -o results/map2heros -r data/ref/heros/star_index
done
sbatch mcic-scripts/convert/merge-bam.sh results/map2heros results/map2heros/merged/merged.bam
sbatch mcic-scripts/assembly_trans/trinity_guided.sh -i results/map2heros/merged/merged.bam -o results/trinity/heros/trinity_guided

# Run Trans-ABySS -- separately for multiple kmer size
for k in 21 23 25 27 29 31 33 35; do
    sbatch mcic-scripts/assembly_trans/transabyss.sh -i "$fqdir_norm" -I heros_transabyss_k"$k" -o results/transabyss/heros/k"$k" -k "$k"
done

# Run SPAdes -- will automatically run for multiple kmer sizes
mkdir -p results/spades/heros
find results/orna/heros -name "*R1*fastq.gz" -size +1k > results/spades/heros/R1.txt
bash mcic-scripts/assembly/spades_R1files2yaml.sh results/spades/heros/R1.txt results/spades/heros/R1.yaml
sbatch mcic-scripts/assembly/spades.sh -y results/spades/heros/R1.yaml -o results/spades/heros/ -m rna -a "--ss rf"


# ASSEMBLY QC -------------------------------------------------------------------
assembly=results/trinity/taurulus/Trinity.fasta

# Run Trinity-stats -- for 1 assembly at a time
source activate /users/PAS0471/jelmer/miniconda3/envs/trinity-env
TrinityStats.pl "$assembly" > results/trinity/taurulus/Trinity_assembly.metrics

# Run BUSCO -- for 1 assembly at a time
sbatch mcic-scripts/assembly_trans/busco.sh -i "$assembly" -o results/busco/taurulus -d insecta_odb10

# rnaQUAST -- for all assemblies at once
assemblies=$(find results/trinity -name "Trinity*fasta" | head | tr '\n' ' ')   # Get a space-separated list of assembly FASTA files
sbatch mcic-scripts/assembly_trans/rnaquast.sh -i "$assemblies" -o results/rnaquast

# Run Detonate -- for 1 assembly at a time
# (Include the same FASTQ files that were also used for assembly with the -I option)
sbatch mcic-scripts/assembly_trans/detonate.sh -i "$assembly" -o results/detonate/taurulus -I "$fqdir_norm"


# ASSEMBLY MERGING AND FILTERING -----------------------------------------------
# Add assembly-IDs to contigs and concatenate assemblies
mkdir -p results/assemblies/merged/

asm_trin_gg=results/trinity/heros/trinity_guided/Trinity-GG.fasta
asm_trin=results/trinity/heros/trinity_denovo/Trinity.fasta
asm_spades=results/spades/heros/transcripts.fasta
asm_ab_21=results/transabyss/heros_transabyss_k21-final.fa
asm_ab_23=results/transabyss/heros_transabyss_k23-final.fa

sed 's/>/>tringg_/' "$asm_trin_gg" > results/assemblies/trinity_gg.fasta
sed 's/>/>trin_/' "$asm_trin" > results/assemblies/trinity.fasta
sed 's/>/>spades_/' "$asm_spades" > results/assemblies/spades.fasta
sed 's/>/>trab21_/' "$asm_ab_21" > results/assemblies/transabyss_k21.fasta
sed 's/>/>trab23_/' "$asm_ab_23" > results/assemblies/transabyss_k23.fasta
cat results/assemblies/*fasta > results/assemblies/merged/merged.fasta

# Merge and filter transcripts with EviGene
sbatch mcic-scripts/assembly_trans/evigene.sh -i results/assemblies/merged/merged.fasta -o results/evigene


# ASSEMBLY ANNOTATION ----------------------------------------------------------
# Only take _primary_ transcript for each gene
assembly_all=results/evigene/okayset/merged.okay.mrna
assembly_primary=results/evigene/okayset/merged.okay.primarytrans.fasta
awk -v RS='>' '/t1 type/ { print ">" $0 }' "$assembly_all" > "$assembly_primary"

# Configure EnTap
refseq_db=/fs/project/PAS0471/jelmer/refdata/refseq/complete/refseq_complete.faa
uniprot_db=/fs/project/PAS0471/jelmer/refdata/uniprot/uniprot_sprot.fasta
dbs="$refseq_db $uniprot_db"
config=results/entap/entap_config.ini
sbatch mcic-scripts/assembly_trans/entap_config.sh -d "$dbs" -c "$config" -o results/entap

# Run Entap
#> Note: When rerunning Entap in a new dir without rerunning the config script, copy: (1) the config.ini file, (2) the `bin` dir, (3) the `databases` dir
#> Note: file paths need to be absolute, hence addition of $PWD
assembly="$PWD"/results/evigene/okayset/merged.okay.primarytrans.fasta
refseq_db="$PWD"/results/entap/bin/refseq_complete.dmnd
uniprot_db="$PWD"/results/entap/bin/uniprot_sprot.dmnd
dbs="$refseq_db $uniprot_db"
config="$PWD"/results/entap/entap_config.ini
sbatch mcic-scripts/assembly_trans/entap.sh -i "$assembly" -d "$dbs" -c "$config" -o "$PWD"/results/entap

# Process and check EnTap results (script will also produce a gene-to-transcript map)
entap_dir=results/entap
in_1trans=results/evigene/okayset/merged.okay.primarytrans.fasta
in_alltrans=results/evigene/okayset/merged.okay.mrna
assembly_out=results/entap/final_ed/assembly.fa
sbatch mcic-scripts/assembly_trans/entap_process.sh -i "$entap_dir" -o "$assembly_out" -a "$in_1trans" -A "$in_alltrans"

# 2022-08-22: Get gene lengths for GO
module load miniconda3; source activate /fs/project/PAS0471/jelmer/conda/bioawk
bioawk -c fastx '{ print $name, length($seq) }' results/entap/final_ed/assembly.fa > results/entap/final_ed/genelen.tsv

# 2022-09-10: Check KEGG annotations
ann=results/entap/final_ed/final_annotations_no_contam_lvl0.tsv
head -n1 $ann | tr "\t" "\n" | cat -n
awk -F "\t" '$14 ~ "a" {print $1, $14, $34}' $ann
awk -F "\t" -v OFS="\t" '$34 != "" {print $1, $14, $34}' $ann | head

    
# GENE AND TRANSCRIPT QUANTIFICATION -------------------------------------------
assembly=results/entap/final_ed/assembly.fa

# Transcript list for making gene-to-transcript map for Kallisto import
grep ">" "$assembly" | sed -E 's/>(.*) type=.*/\1/' > results/evigene/okayset/transcriptIDs.txt

# Run Kallisto
sbatch mcic-scripts/rnaseq/kallisto_index.sh -i "$assembly" -o results/kallisto/assembly.idx
for R1 in results/kraken/heros_postviral/unclassified/*R1*fastq.gz; do
    sbatch mcic-scripts/rnaseq/kallisto_quant.sh -i "$R1" -o results/kallisto -r results/kallisto/assembly.idx
done
