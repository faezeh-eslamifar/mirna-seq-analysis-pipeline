#!/bin/bash

# Root directories
ROOT_DIR="$(pwd)"
INPUT_DIR="${ROOT_DIR}/input"
OUTPUT_DIR="${ROOT_DIR}/output"
GENOME_DIR="${ROOT_DIR}/genomes"
BOWTIE_INDEX="${GENOME_DIR}/INDEX/hsa_genome"
HAIRPIN_FILE="${GENOME_DIR}/hairpin/hairpin.fa"
MIRDEEP2_DIR="${MIRDEEP2_DIR:?Please set MIRDEEP2_DIR to your miRDeep2 installation directory}"

# Create necessary directories
mkdir -p "${OUTPUT_DIR}/fastqc_results/raw"
mkdir -p "${OUTPUT_DIR}/fastqc_results/trimmed"
mkdir -p "${OUTPUT_DIR}/trimmed"
mkdir -p "${OUTPUT_DIR}/bowtie"
mkdir -p "${OUTPUT_DIR}/mirdeep2_input"
mkdir -p "${OUTPUT_DIR}/mirdeep2_output"
mkdir -p "${GENOME_DIR}/INDEX"
mkdir -p "${GENOME_DIR}/mirna"
mkdir -p "${GENOME_DIR}/hairpin"


# Library-specific adapter sequences: replace these placeholders before running.
ADAPTER_5P="REPLACE_WITH_YOUR_5_PRIME_ADAPTER_SEQUENCE"
ADAPTER_3P="REPLACE_WITH_YOUR_3_PRIME_ADAPTER_SEQUENCE"

# Step 1: Download miRNA reference sequences and annotations
echo "Checking for existing miRNA reference and annotations..."

# Download and decompress the human GRCh38 reference genome
if [ ! -f "${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
    echo "Downloading Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz..."

    wget -O "${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
        "https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
        || { echo "Failed to download the reference genome"; exit 1; }

    gunzip "${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
        || { echo "Failed to decompress the reference genome"; exit 1; }

    echo "Reference genome downloaded and decompressed successfully."
else
    echo "Homo_sapiens.GRCh38.dna.primary_assembly.fa already exists."
fi

if [ ! -f "${GENOME_DIR}/mirna/mature_miRNA.fa" ]; then
    echo "Downloading mature_miRNA.fa..."
    wget -O "${GENOME_DIR}/mirna/mature_miRNA.fa" "https://www.mirbase.org/download/mature.fa" || { echo "Failed to download mature_miRNA.fa"; exit 1; }
else
    echo "mature_miRNA.fa already exists."
fi

if [ ! -f "${HAIRPIN_FILE}" ]; then
    echo "Downloading hairpin.fa..."
    wget -O "${HAIRPIN_FILE}" "https://www.mirbase.org/download/hairpin.fa" || { echo "Failed to download hairpin.fa"; exit 1; }
else
    echo "hairpin.fa already exists."
fi

# Step 2: Filter mature_miRNA.fa for human sequences
echo "Filtering mature_miRNA.fa for human sequences..."
if [ ! -f "${GENOME_DIR}/mirna/mature_hsa.fa" ]; then
    perl "${MIRDEEP2_DIR}/bin/extract_miRNAs.pl" \
        "${GENOME_DIR}/mirna/mature_miRNA.fa" hsa > "${GENOME_DIR}/mirna/mature_hsa.fa"
else
    echo "mature_hsa.fa already exists."
fi

# Step 3: Filter hairpin.fa for human sequences
echo "Filtering hairpin.fa for human sequences..."
if [ ! -f "${GENOME_DIR}/hairpin/hairpin_hsa.fa" ]; then
    perl "${MIRDEEP2_DIR}/bin/extract_miRNAs.pl" \
        "${HAIRPIN_FILE}" hsa > "${GENOME_DIR}/hairpin/hairpin_hsa.fa" || { echo "Failed to filter human sequences from hairpin.fa"; exit 1; }
else
    echo "hairpin_hsa.fa already exists."
fi

# Step 4: Preprocess genome file
echo "Preprocessing genome file for miRDeep2 compatibility..."
if [ ! -f "${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly_cleaned.fa" ]; then
    perl "${MIRDEEP2_DIR}/bin/remove_white_space_in_id.pl" \
        "${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa" > "${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly_cleaned.fa" || { echo "Failed to clean whitespaces in genome file"; exit 1; }
else
    echo "Filtered genome file already exists."
fi

# Step 5: Build Bowtie index
echo "Checking for existing Bowtie index..."
if [ ! -f "${BOWTIE_INDEX}.1.ebwt" ]; then
    echo "Building Bowtie index for miRNA reference..."
    bowtie-build "${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly_cleaned.fa" "${BOWTIE_INDEX}" || { echo "Failed to build Bowtie index"; exit 1; }
else
    echo "Bowtie index already exists."
fi

# Step 6: Process FASTQ files
for FASTQ_FILE in "${INPUT_DIR}"/*.fastq.gz; do
    SAMPLE_NAME=$(basename "$FASTQ_FILE" .fastq.gz)
    echo "Processing sample: $SAMPLE_NAME"

    # Step 6.1: Run FastQC on raw data
    fastqc "$FASTQ_FILE" -o "${OUTPUT_DIR}/fastqc_results/raw" || { echo "FastQC failed for raw data"; exit 1; }

    # Step 6.2: Remove poly-A tails
    cutadapt -a "A{10}" -a "G{10}" -o "${OUTPUT_DIR}/trimmed/${SAMPLE_NAME}_trimmed1.fastq.gz" "$FASTQ_FILE" || { echo "Cutadapt failed"; exit 1; }

    # Step 6.3: Trim reads with adapters
    cutadapt -a "$ADAPTER_3P" -g "$ADAPTER_5P" \
        -o "${OUTPUT_DIR}/trimmed/${SAMPLE_NAME}_trimmed.fastq.gz" \
        --minimum-length 17 --trim-n --quality-cutoff 30,30 \
        "${OUTPUT_DIR}/trimmed/${SAMPLE_NAME}_trimmed1.fastq.gz" || { echo "Cutadapt failed"; exit 1; }
    rm "${OUTPUT_DIR}/trimmed/${SAMPLE_NAME}_trimmed1.fastq.gz"

    # Step 6.4: Run FastQC on trimmed data
    fastqc "${OUTPUT_DIR}/trimmed/${SAMPLE_NAME}_trimmed.fastq.gz" -o "${OUTPUT_DIR}/fastqc_results/trimmed" || { echo "FastQC failed for trimmed data"; exit 1; }

    # Step 7: Align reads with Bowtie
    gunzip -c "${OUTPUT_DIR}/trimmed/${SAMPLE_NAME}_trimmed.fastq.gz" > "${OUTPUT_DIR}/trimmed/${SAMPLE_NAME}_trimmed.fastq"
    bowtie -q -p 4 -l 17 -n 1 -x "${BOWTIE_INDEX}" "${OUTPUT_DIR}/trimmed/${SAMPLE_NAME}_trimmed.fastq" "${OUTPUT_DIR}/bowtie/${SAMPLE_NAME}.sam" || { echo "Bowtie alignment failed"; exit 1; }
    

    # Step 8: Convert SAM to ARF
    grep -v "^@" "${OUTPUT_DIR}/bowtie/${SAMPLE_NAME}.sam" > "${OUTPUT_DIR}/bowtie/${SAMPLE_NAME}_cleaned.sam"
    perl "${MIRDEEP2_DIR}/bin/convert_bowtie_output.pl" \
        "${OUTPUT_DIR}/bowtie/${SAMPLE_NAME}_cleaned.sam" > "${OUTPUT_DIR}/mirdeep_input/${SAMPLE_NAME}.arf"
    rm "${OUTPUT_DIR}/bowtie/${SAMPLE_NAME}.sam"

    # Step 9: Rewrite ARF file read IDs
    awk '{ unique_id = NR; count = "x1"; $1 = "seq_" unique_id "_" count; print }' OFS='\t' \
        "${OUTPUT_DIR}/mirdeep_input/${SAMPLE_NAME}.arf" > "${OUTPUT_DIR}/mirdeep_input/${SAMPLE_NAME}_fixed.arf"
    rm "${OUTPUT_DIR}/mirdeep_input/${SAMPLE_NAME}.arf"

    # Step 10: Convert trimmed FASTQ to FASTA
    perl "${MIRDEEP2_DIR}/bin/fastq2fasta.pl" \
        "${OUTPUT_DIR}/trimmed/${SAMPLE_NAME}_trimmed.fastq" > "${OUTPUT_DIR}/mirdeep_input/${SAMPLE_NAME}.fa"
        
    rm "${OUTPUT_DIR}/trimmed/${SAMPLE_NAME}_trimmed.fastq"

    # Rewrite FASTA file read IDs with counts
    awk '
    BEGIN { delete sequence_counts; }
    {
        if ($0 ~ /^>/) { header = $0; }
        else {
            sequence = $0;
            sequence_counts[sequence]++;
        }
    }
    END {
        for (seq in sequence_counts) {
            count = sequence_counts[seq];
            print ">seq_" NR "_x" count;
            print seq;
        }
    }' "${OUTPUT_DIR}/mirdeep_input/${SAMPLE_NAME}.fa" > "${OUTPUT_DIR}/mirdeep_input/${SAMPLE_NAME}_fixed.fa"
    rm "${OUTPUT_DIR}/mirdeep_input/${SAMPLE_NAME}.fa"

    # Step 11: Run miRDeep2 analysis
    miRDeep2.pl "${OUTPUT_DIR}/mirdeep_input/${SAMPLE_NAME}_fixed.fa" \
        "${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly_cleaned.fa" \
        "${OUTPUT_DIR}/mirdeep_input/${SAMPLE_NAME}_fixed.arf" \
        "${GENOME_DIR}/mirna/mature_hsa.fa" none "${GENOME_DIR}/hairpin/hairpin_hsa.fa" \
        -t Human -p > "${OUTPUT_DIR}/mirdeep2/${SAMPLE_NAME}_mirdeep2_output" || { echo "miRDeep2 analysis failed"; }
done

echo "miRNA analysis completed."



