
#######################################################################
############     Xenome	(First try)					############# 
#######################################################################

# ---- Sequence -----
#cp /share/apps/STAR/indexes/Mus_musculus/Mus_musculus.GRCm38.89.gtf  Mus_musculus.GRCm38.89.gtf
#cp /share/apps/STAR/indexes/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly.fa .

# Change header of the file in place to avoid bug met before with faidx " Ignoring duplicate sequence "
#gawk -i inplace  '/^>(.*)/{print $1"_M"; next}{print }' Mus_musculus.GRCm38.dna.primary_assembly.fa

#cat /share/apps/STAR/indexes/Homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa Mus_musculus.GRCm38.dna.primary_assembly.fa > H.GRCh38.83_M.GRCm38.89.primary_assembly.fa

#rm Mus_musculus.GRCm38.dna.primary_assembly.fa

#samtools faidx H.GRCh38.83_M.GRCm38.89.primary_assembly.fa
#faidx H.GRCh38.83_M.GRCm38.89.primary_assembly.fa -i chromsizes > yourgenome.genome.sizes
#faToTwoBit H.GRCh38.83_M.GRCm38.89.primary_assembly.fa H.GRCh38.83_M.GRCm38.89.primary_assembly.fa.2bit


# ---- Annotation -----
#cp /share/apps/STAR/indexes/Mus_musculus/Mus_musculus.GRCm38.89.gtf  Mus_musculus.GRCm38.89.gtf

# Change also chromosome id for gtf 
#gawk -i inplace  'OFS="\t" {if (NR > 5) $1=$1"_M"; print}'  Mus_musculus.GRCm38.89.gtf

#cat /share/apps/STAR/indexes/Homo_sapiens/Homo_sapiens.GRCh38.83.gtf Mus_musculus.GRCm38.89.gtf > H.GRCh38.83_M.GRCm38.89.gtf

#rm Mus_musculus.GRCm38.89.gtf

# ---- STAR INDEX -----

#/share/apps/STAR/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir /data/villemin/genome/xenome/ --genomeFastaFiles /data/villemin/genome/xenome/H.GRCh38.83_M.GRCm38.89.primary_assembly.fa  --runThreadN 30 --sjdbGTFfile /data/villemin/genome/xenome/H.GRCh38.83_M.GRCm38.89.gtf

# ------ Super fast parallel fastqc ------  #
 
#ls -1 /data/villemin/data/toulouse/scRNAseqPDX/*.fastq.gz | grep "R1\|R2" | parallel  fastqc  {}
#mv /data/villemin/data/toulouse/scRNAseqPDX/*.html ./report
#mv /data/villemin/data/toulouse/scRNAseqPDX/*.zip ./report


#####################################################################
############     Xenome	(10X Genomicss)					############# 
##################################################################### 

####  Cell Ranger #####
#conda install -c conda-forge rust
#conda install -c conda-forge go
#conda install -c anaconda bzip2
#conda install -c anaconda make
#conda install -c anaconda clang

#/data/villemin/soft/cellranger-6.0.0/bin/cellranger

####eCellRanger Count (quantitates a single run)
# cellranger count --id=COURSE \
#--transcriptome=/bi/apps/cellranger/references/GRCh38/ \
#--fastqs=/bi/home/andrewss/10X/ \
#--localcores=8 \
#--localmem=32

#### CellRanger aggr (merges multiple runs)
# cellranger aggr --id=MERGED \
#csv=merge_me.csv \
#normalize=mapped
#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation
#https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#GRCh38mm10_2020A

#https://github.com/alexdobin/STAR/issues/346
export PATH=/data/villemin/soft/cellranger-6.0.0:$PATH

# Genome metadata
genome="GRCh38"
version="2020-A"


# Set up source and build directories
build="GRCh38-2020-A_build"
mkdir -p "$build"


# Download source files if they do not exist in reference_sources/ folder
source="reference_sources"
mkdir -p "$source"


fasta_url="http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fasta_in="${source}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"
gtf_in="${source}/gencode.v32.primary_assembly.annotation.gtf"


if [ ! -f "$fasta_in" ]; then
    curl -sS "$fasta_url" | zcat > "$fasta_in"
fi
if [ ! -f "$gtf_in" ]; then
    curl -sS "$gtf_url" | zcat > "$gtf_in"
fi


# Modify sequence headers in the Ensembl FASTA to match the file
# "GRCh38.primary_assembly.genome.fa" from GENCODE. Unplaced and unlocalized
# sequences such as "KI270728.1" have the same names in both versions.
#
# Input FASTA:
#   >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
#
# Output FASTA:
#   >chr1 1
fasta_modified="$build/$(basename "$fasta_in").modified"
# sed commands:
# 1. Replace metadata after space with original contig name, as in GENCODE
# 2. Add "chr" to names of autosomes and sex chromosomes
# 3. Handle the mitochrondrial chromosome
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"


# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages
#
# Input GTF:
#     ... gene_id "ENSG00000223972.5"; ...
# Output GTF:
#     ... gene_id "ENSG00000223972"; gene_version "5"; ...
gtf_modified="$build/$(basename "$gtf_in").modified"
# Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"


# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
# - These filters are relevant only to GTF files from GENCODE. The GTFs from
#   Ensembl release 98 have the following differences:
#   - The names "gene_biotype" and "transcript_biotype" are used instead of
#     "gene_type" and "transcript_type".
#   - Readthrough transcripts are present but are not marked with the
#     "readthrough_transcript" tag.
#   - Only the X chromosome versions of genes in the pseudoautosomal regions
#     are present, so there is no "PAR" tag.
BIOTYPE_PATTERN=\
"(protein_coding|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""


# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "PAR" tag (only present for Y chromosome PAR)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"


# Filter the GTF file based on the gene allowlist
gtf_filtered="${build}/$(basename "$gtf_in").filtered"
# Copy header lines beginning with "#"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"
# Filter to the gene allowlist
grep -Ff "${build}/gene_allowlist" "$gtf_modified" \
    >> "$gtf_filtered"


# Create reference package
cellranger mkref --ref-version="$version" \
    --genome="$genome" --fasta="$fasta_modified" --genes="$gtf_filtered"

Mouse reference, mm10
# Genome metadata
genome="mm10"
version="2020-A"


# Set up source and build directories
build="mm10-2020-A_build"
mkdir -p "$build"


# Download source files if they do not exist in reference_sources/ folder
source="reference_sources"
mkdir -p "$source"


fasta_url="http://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
fasta_in="${source}/Mus_musculus.GRCm38.dna.primary_assembly.fa"
gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz"
gtf_in="${source}/gencode.vM23.primary_assembly.annotation.gtf"


if [ ! -f "$fasta_in" ]; then
    curl -sS "$fasta_url" | zcat > "$fasta_in"
fi
if [ ! -f "$gtf_in" ]; then
    curl -sS "$gtf_url" | zcat > "$gtf_in"
fi


# Modify sequence headers in the Ensembl FASTA to match the file
# "GRCm38.primary_assembly.genome.fa" from GENCODE. Unplaced and unlocalized
# sequences such as "GL456210.1" have the same names in both versions.
#
# Input FASTA:
#   >1 dna:chromosome chromosome:GRCm38:1:1:195471971:1 REF
#
# Output FASTA:
#   >chr1 1
fasta_modified="$build/$(basename "$fasta_in").modified"
# sed commands:
# 1. Replace metadata after space with original contig name, as in GENCODE
# 2. Add "chr" to names of autosomes and sex chromosomes
# 3. Handle the mitochrondrial chromosome
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"


# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages
#
# Input GTF:
#     ... gene_id "ENSMUSG00000102693.1"; ...
# Output GTF:
#     ... gene_id "ENSMUSG00000102693"; gene_version "1"; ...
gtf_modified="$build/$(basename "$gtf_in").modified"
# Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"


# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
# - These filters are relevant only to GTF files from GENCODE. The GTFs from
#   Ensembl release 98 have the following differences:
#   - The names "gene_biotype" and "transcript_biotype" are used instead of
#     "gene_type" and "transcript_type".
#   - Readthrough transcripts are present but are not marked with the
#     "readthrough_transcript" tag.
BIOTYPE_PATTERN=\
"(protein_coding|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""


# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"


# Filter the GTF file based on the gene allowlist
gtf_filtered="${build}/$(basename "$gtf_in").filtered"
# Copy header lines beginning with "#"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"
# Filter to the gene allowlist
grep -Ff "${build}/gene_allowlist" "$gtf_modified" \
    >> "$gtf_filtered"


# Create reference package
cellranger mkref --ref-version="$version" \
    --genome="$genome" --fasta="$fasta_modified" --genes="$gtf_filtered"

Human and mouse reference dataset, GRCh38 and mm10
#################### SETUP ####################


human_genome="GRCh38"
mouse_genome="mm10"
version="2020-A"


build="GRCh38_and_mm10-2020-A_build"
mkdir -p "$build"


# Download source files if they do not exist in reference_sources/ folder
source="reference_sources"
mkdir -p "$source"


human_fasta_url="http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
human_fasta_in="${source}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
human_gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"
human_gtf_in="${source}/gencode.v32.primary_assembly.annotation.gtf"
mouse_fasta_url="http://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
mouse_fasta_in="${source}/Mus_musculus.GRCm38.dna.primary_assembly.fa"
mouse_gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz"
mouse_gtf_in="${source}/gencode.vM23.primary_assembly.annotation.gtf"


if [ ! -f "$human_fasta_in" ]; then
    curl -sS "$human_fasta_url" | zcat > "$human_fasta_in"
fi
if [ ! -f "$human_gtf_in" ]; then
    curl -sS "$human_gtf_url" | zcat > "$human_gtf_in"
fi
if [ ! -f "$mouse_fasta_in" ]; then
    curl -sS "$mouse_fasta_url" | zcat > "$mouse_fasta_in"
fi
if [ ! -f "$mouse_gtf_in" ]; then
    curl -sS "$mouse_gtf_url" | zcat > "$mouse_gtf_in"
fi


# String patterns used for both genomes
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"


BIOTYPE_PATTERN=\
"(protein_coding|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""


#################### HUMAN ####################
# Please see the GRCh38-2020-A build documentation for details on these steps.


# Process FASTA -- translate chromosome names
human_fasta_modified="$build/$(basename "$human_fasta_in").modified"
cat "$human_fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$human_fasta_modified"


# Process GTF -- split Ensembl IDs from version suffixes
human_gtf_modified="$build/$(basename "$human_gtf_in").modified"
cat "$human_gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$human_gtf_modified"


# Process GTF -- filter based on gene/transcript tags
cat "$human_gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"


human_gtf_filtered="${build}/$(basename "$human_gtf_in").filtered"
grep -E "^#" "$human_gtf_modified" > "$human_gtf_filtered"
grep -Ff "${build}/gene_allowlist" "$human_gtf_modified" \
    >> "$human_gtf_filtered"


#################### MOUSE ####################
# Please see the mm10-2020-A build documentation for details on these steps.


# Process FASTA -- translate chromosome names
mouse_fasta_modified="$build/$(basename "$mouse_fasta_in").modified"
cat "$mouse_fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$mouse_fasta_modified"


# Process GTF -- split Ensembl IDs from version suffixes
mouse_gtf_modified="$build/$(basename "$mouse_gtf_in").modified"
cat "$mouse_gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$mouse_gtf_modified"


# Process GTF -- filter based on gene/transcript tags
cat "$mouse_gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"


mouse_gtf_filtered="${build}/$(basename "$mouse_gtf_in").filtered"
grep -E "^#" "$mouse_gtf_modified" > "$mouse_gtf_filtered"
grep -Ff "${build}/gene_allowlist" "$mouse_gtf_modified" \
    >> "$mouse_gtf_filtered"


#################### MKREF ####################


cellranger mkref --ref-version="$version" \
    --genome="$human_genome" --fasta="$human_fasta_modified" --genes="$human_gtf_filtered" \
    --genome="$mouse_genome" --fasta="$mouse_fasta_modified" --genes="$mouse_gtf_filtered"