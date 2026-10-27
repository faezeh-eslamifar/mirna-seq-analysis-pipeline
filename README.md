# mirna-seq-analysis-pipeline
A Bash-based pipeline for preprocessing and analyzing human miRNA-seq data using Cutadapt, Bowtie, and miRDeep2.

"Overview"

This pipeline automates the preprocessing and analysis of human miRNA sequencing data, starting from raw FASTQ files and generating input files and analysis outputs for miRDeep2.

The workflow includes quality control, adapter trimming, reference genome preparation, genome alignment, and miRNA analysis. Downstream differential expression analysis is performed separately.

"Workflow"

Reference genome preparation: Downloads the human GRCh38 primary assembly genome from Ensembl if it is not already available, then preprocesses it for miRDeep2 compatibility.

miRNA reference preparation: Downloads mature miRNA and hairpin sequences from miRBase if they are not already available, then extracts human-specific sequences.

Genome indexing: Builds a Bowtie index for the human reference genome if it does not already exist.

Raw data quality control: Runs FastQC on raw FASTQ files.

Read preprocessing: Removes specified poly-A/poly-G sequences, trims library-specific adapters, applies quality trimming, and discards reads shorter than 17 nucleotides using Cutadapt.

Genome alignment: Aligns the processed reads to the human reference genome using Bowtie.

Input preparation: Converts and reformats alignment and read files into ARF and FASTA formats for miRDeep2.

miRNA analysis: Runs miRDeep2 separately for each sample.

"Requirements"

The pipeline requires the following software and tools:

Bash - FastQC - Cutadapt - Bowtie - miRDeep2

Perl and the required Perl modules - Wget

AWK and standard Unix command-line utilities

A Linux or compatible Unix-like environment is required.

"Input Data"

Place the compressed FASTQ files (*.fastq.gz) in the input/ directory.

Each FASTQ file is processed individually. The sample name is derived from the filename by removing the .fastq.gz extension.

"Reference Resources"

"Human Reference Genome"

The pipeline downloads the human GRCh38 primary assembly genome from Ensembl (release 116) if the uncompressed reference file is not already present.

The expected file path is:

genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa

The downloaded genome is decompressed automatically. The pipeline then preprocesses the genome file and builds the Bowtie index if the corresponding output files are not already available.

"miRNA Reference Sequences"

The pipeline downloads the mature miRNA and hairpin reference files from miRBase if they are not already present. Human-specific sequences are extracted using the extract_miRNAs.pl script included with miRDeep2.

The resulting reference files are stored in:

genomes/mirna/mature_hsa.fa

genomes/hairpin/hairpin_hsa.fa

"Configuration"

"miRDeep2 Installation"

Set the MIRDEEP2_DIR environment variable to the local miRDeep2 installation directory containing the required Perl scripts.

Example:

export MIRDEEP2_DIR="/path/to/miRDeep2"

"Adapter Sequences"

The 5′ and 3′ adapter sequences are library-specific. Before running the pipeline, replace the placeholder values for ADAPTER_5P and ADAPTER_3P near the beginning of the script with the adapter sequences corresponding to your sequencing library.

ADAPTER_5P="REPLACE_WITH_YOUR_5_PRIME_ADAPTER_SEQUENCE"
ADAPTER_3P="REPLACE_WITH_YOUR_3_PRIME_ADAPTER_SEQUENCE"

The placeholder values must be replaced before execution.

"Directory Structure"

.
├── input/

│   └── sample.fastq.gz

├── genomes/

│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa

│   ├── Homo_sapiens.GRCh38.dna.primary_assembly_cleaned.fa

│   ├── INDEX/

│   ├── hairpin/

│   └── mirna/

├── output/

│   ├── fastqc_results/

│   │   ├── raw/

│   │   └── trimmed/

│   ├── trimmed/

│   ├── bowtie/

│   ├── mirdeep2_input/

│   └── mirdeep2_output/

└── mirna_seq_pipeline.sh

"Usage"

1. Clone or download the repository.

2. Place the raw compressed FASTQ files in the input/ directory.

3. Set MIRDEEP2_DIR to the local miRDeep2 installation directory.

4. Replace the adapter placeholders in the script with the appropriate library-specific sequences.

5. Open a terminal and navigate to the project root directory.

6. Run the pipeline:

export MIRDEEP2_DIR="/path/to/miRDeep2"
bash mirna_seq_pipeline.sh

Important: Run the script from the project root directory because the input, output, and genome paths are defined relative to the current working directory.

"Output"

The pipeline generates the following outputs:

FastQC reports: Quality control reports for raw and trimmed reads.

Trimmed FASTQ files: Reads after preprocessing and adapter trimming.

Bowtie-related files: Intermediate alignment files.

miRDeep2 input files: Reformatted FASTA and ARF files.

miRDeep2 output files: Analysis outputs generated separately for each sample.

All outputs are stored in the output/ directory, while reference files and the Bowtie index are stored in genomes/.

"Scope and Limitations"

The pipeline is designed specifically for human miRNA-seq data and uses miRNA-oriented preprocessing and alignment settings.

Each sample is processed individually.

Differential expression analysis is not included in this pipeline and is performed separately.

Adapter sequences must be configured according to the sequencing library protocol.

The pipeline uses the human GRCh38 primary assembly genome from Ensembl release 116.

Reference resources and software versions should be documented to support reproducibility.

The pipeline should be tested with the intended dataset and software environment before being applied to other studies.

"Author"

Developed as a bioinformatics scripting project for human miRNA-seq data analysis.
