# ChiRA - Chimeric Read Analyzer

**Version**: 1.4.14 (Modified with performance optimizations and parallel computing support)

ChiRA is a set of tools to analyze RNA-RNA interactome experimental data such as CLASH, CLEAR-CLIP, PARIS, SPLASH etc. Following are the descriptions of the each tool. Here we provide descriptions about the input and ouptput files. For the detailed description of the other parameters please look at the help texts of tools.

**Note**: fastChiRA is a modified version of ChiRA (based on v1.4.3) with significant performance optimizations and new features. The original code is licensed under GPL v3, and this modified version maintains the same license. All changes are documented in the "Recent Improvements" section below and in [CHANGELOG.md](CHANGELOG.md).

## Version History
- **v1.4.14** (Current, 2026-02-21): Version bump; documentation updates (final table column names, CLI args, future/future.apply removed)
- **v1.4.13** (2026-02-21): Batchtools/IntaRNA fixes: wait for all batches (including last), POSIXct-safe job status counting and manual polling in R scripts, CSV line-ending handling in chira_extract; IntaRNA sequences from loci_seqs.pkl (no seq1/seq2 in CSV); merge_intarna_into_chimeras.py for post-IntaRNA merge when main job times out
- **v1.4.11** (2026-02-17): MPIRE made required dependency for optimal multiprocessing performance, removed fallback code, improved memory efficiency (50-90% reduction) and startup time (2-3x faster)
- **v1.4.10** (2026-02-15): Fixed batchtools submission issues (template path handling, JSON parsing), ensured all paths are absolute for cluster jobs, and refactored scripts for better code organization
- **v1.4.9** (2026-02-15): Added `--parallel_chunks` parameter for configurable chunk parallelism in chira_map.py
- **v1.4.8** (2026-02-15): Improved chunk-based parallelization for very large transcript counts (e.g., human genome with 387K+ transcripts), variable naming consistency improvements, and bug fixes in chira_merge.py
- **v1.4.7** (2026-02-15): New utility scripts (extract_transcripts_from_genome.py), U→T conversion in download_mirbase_mature.py, gffread support, and removal of deprecated scripts
- **v1.4.6** (2026-02-15): Multiprocessing improvements, adaptive I/O buffer sizing, FASTA chunking, I/O bottleneck detection, and code refactoring
- **v1.4.5** (2026-02-02): Parallel computing support, I/O optimizations, automatic memory management, and enhanced performance
- **v1.4.4** (2026-01-26): Performance optimizations, BEDTools compatibility, sample name support, comprehensive bug fixes, and new utility scripts
- **v1.4.3** (Previous): Original version from GitHub (https://github.com/pavanvidem/chira)

See [CHANGELOG.md](CHANGELOG.md) for detailed change history.

## Recent Improvements

### v1.4.11 (2026-02-17) - MPIRE Required & Critical Bug Fixes

**Critical Bug Fixes:**
- Fixed CRL ID validation and EM algorithm return bugs in `chira_quantify.py`
- Fixed alternate alignment filtering for stranded RNA-seq in `chira_map.py`
- Fixed field index mismatches and MPIRE argument passing in `chira_extract.py`
- Collapsed 8 separate `hybridization_genomic_coordinates` columns into a single column

**Multiprocessing:**
- **MPIRE is now required** (moved from optional to `install_requires`)
- All parallel processing uses MPIRE WorkerPool with shared objects
- Added `start_method='fork'` for copy-on-write semantics (50-90% memory reduction)
- Benefits: 50-90% memory reduction, 2-3x faster startup

### v1.4.7-1.4.10 (2026-02-15) - Utility Scripts, Code Refactoring & Batchtools

**New Features:**
- **extract_transcripts_from_genome.py**: New utility script using gffread (replaces `remove_mirna_hairpin_from_fasta.py`)
- **download_mirbase_mature.py**: Added automatic U→T conversion for ChiRA compatibility
- **chira_map.py**: Added `--parallel_chunks` parameter and FASTA chunking support
- **Batchtools support**: Enhanced HPC cluster job submission with absolute path handling

**Code Improvements:**
- All scripts refactored for better code organization (modular functions, consistent structure)
- Improved chunk-based parallelization for very large datasets (e.g., 387K+ transcripts)

### v1.4.6 (2026-02-15) - Multiprocessing Evolution & I/O Optimizations

**Multiprocessing Evolution:**
- ThreadPoolExecutor → ProcessPoolExecutor → MPIRE (bypasses GIL, 2-8x faster)
- Adaptive buffer sizing (8-16MB, 10-50x I/O improvement for large files)
- FASTA chunking support in `chira_map.py` for better memory efficiency
- Progress tracking for large BAM files (reports every 1M reads)

### v1.4.5 (2026-02-02) - Parallel Computing Introduction

**Initial Parallel Computing:**
- Multi-threading support in `chira_quantify.py` and `chira_merge.py` (2-4x faster)
- Enhanced multi-threading for external tools (samtools, pysam) in `chira_map.py`
- GNU sort parallel support (2-4x faster sorting)
- I/O optimizations with 2MB buffers (20-40% faster)

### v1.4.4 (2026-01-26) - Performance Optimizations

**Performance Improvements:**
- **3-10x faster** overall processing, especially in `chira_quantify.py`
- **2-5x faster** FASTQ parsing in `chira_collapse.py` (removed Biopython dependency)
- **2-5x faster** CIGAR string parsing with pre-compiled regex patterns
- **20-40% faster** BAM file processing in `chira_map.py`
- **10-50x faster** EM algorithm with optimized memory operations

**New Features:**
- Added `--sample_name` parameter to `chira_extract.py` for customizable output file names
- Output files now use format: `{sample_name}.chimeras.txt`, `{sample_name}.singletons.txt`, `{sample_name}.interactions.txt`
- Added `--gzip` option to `chira_extract.py` for compressing large output files (saves disk space, faster I/O for large files)
- Added header rows to interactions output file with column descriptions
- Added `mirna_position` column to chimeras output indicating read orientation (miRNA_first or miRNA_last)
- Automatic BEDTools version detection - works with both old and new command formats
- New utility scripts for reference file preparation (see Utility Scripts section)
- Docker support with pre-installed dependencies (see Docker Support section)
- Singularity/Apptainer support - see [SINGULARITY_SETUP.md](SINGULARITY_SETUP.md) for detailed instructions

**Compatibility:**
- BEDTools: Automatically supports both `intersectBed`/`fastaFromBed` (old) and `bedtools intersect`/`bedtools getfasta` (new)
- No code changes needed when switching between BEDTools versions

**Bug Fixes:**
- **chira_map.py**: Fixed first iteration bug that wrote empty string to unmapped FASTA file
- **chira_merge.py**: Added zero-length match checks to prevent division by zero errors
- **chira_extract.py**: Fixed multiple bugs including undefined variables, logic errors in strand/chromosome checks, index out-of-bounds in TPM cutoff, and string slicing issues
- **chira_utilities.py**: Fixed `median()` function bug where input wasn't sorted; improved `get_bedtools_command()` to verify command success
- **chira_quantify.py**: Fixed CRL iteration bug that skipped index 0
- Improved file handling with proper context managers throughout
- Fixed BEDTools command compatibility across versions with automatic detection

For complete details with line-by-line changes, please refer to [CHANGELOG.md](CHANGELOG.md).

## Installation

### Recommended: Container Installation

**The easiest way to use ChiRA is with the provided container images**, which include all dependencies pre-installed:

**Pre-built Docker Image Available:**
A functional Docker image is available at: **`docker.io/nemat1976/chiraplus:v0.0.2`**

This image includes all dependencies and is ready to use. Simply pull and run:

**Docker:**
```bash
# Pull the image (first time only)
docker pull docker.io/nemat1976/chiraplus:v0.0.2

# Run ChiRA commands
docker run --rm -v $(pwd)/data:/app/data -v $(pwd)/output:/app/output \
  docker.io/nemat1976/chiraplus:v0.0.2 chira_collapse.py -i data/input.fastq -o output/collapsed.fasta
```

**Singularity/Apptainer (for HPC systems):**
```bash
# Pull image
singularity pull docker://docker.io/nemat1976/chiraplus:v0.0.2

# Run ChiRA commands (entrypoint script handles environment automatically)
singularity exec -B $(pwd)/data:/app/data -B $(pwd)/output:/app/output \
 chiraplus_v0.0.2.sif chira_collapse.py -i data/input.fastq -o output/collapsed.fasta
```

For detailed Singularity/Apptainer setup and usage instructions, see [SINGULARITY_SETUP.md](SINGULARITY_SETUP.md).

### Docker Image Contents

The Docker image includes:
- **Python packages**: biopython, bcbiogff, pysam, requests, pyliftover, psutil
- **Bioinformatics tools**: bwa, samtools, bedtools, gffread, intarna
- **All ChiRA scripts and utilities**: Pre-installed and executable
- **Environment setup**: Proper PATH and PYTHONPATH configuration

**Note:** Optional tools (blockbuster, clan) are not included in the Docker image by default but can be added if needed for specific use cases.

### Running ChiRA in Docker

**Basic usage:**
```bash
docker run --rm -v $(pwd)/data:/app/data -v $(pwd)/output:/app/output \
  docker.io/nemat1976/chiraplus:v0.0.1  chira_collapse.py -i data/input.fastq -o output/collapsed.fasta
```

**Interactive shell:**
```bash
docker run --rm -it -v $(pwd)/data:/app/data -v $(pwd)/output:/app/output \
  docker.io/nemat1976/chiraplus:v0.0.1 bash
```

**Volume mounts:**
```bash
docker run --rm \
  -v /path/to/data:/app/data \
  -v /path/to/output:/app/output \
  -v /path/to/references:/app/references \
  docker.io/nemat1976/chiraplus:v0.0.1 chira_extract.py [options]
```

### Manual Installation

If you prefer to install dependencies manually:

**Core Python packages (required):**
- biopython
- bcbiogff
- pysam
- **mpire** (required for `chira_quantify.py`, `chira_extract.py`, and `chira_merge.py` parallel processing)
  - Enhanced multiprocessing framework for EM algorithm parallelization, chimera extraction, and transcript processing
  - Benefits: 50-90% memory reduction, 2-3x faster startup, better performance
  - Install with: `pip install mpire` or `conda install -c conda-forge mpire`
  - Required for parallel processing (no fallback)
- requests (for `download_ensembl.py`)

**Optional Python packages:**
- **psutil** (highly recommended for optimal performance)
  - **chira_map.py**: Automatic memory optimization for BAM sorting and I/O bottleneck detection
  - **chira_utilities.py**: Adaptive buffer sizing (8-16MB) for 10-50x I/O performance improvement
  - Install with: `pip install psutil` or `conda install psutil`
  - Falls back to safe defaults if not available (2GB per thread for BAM sorting, 8MB buffer for I/O)
- pyliftover (for `download_mirbase_gff3.py` coordinate liftover)

**R packages (for batchtools HPC cluster support):**
- **batchtools** (required for `--use_batchtools` in `chira_map.py` and for hybridization in `chira_extract.py`)
  - Purpose: Submitting chunk-based batch jobs to HPC cluster schedulers (LSF, SLURM, SGE, etc.)
  - Used by: `submit_chunks_batchtools.R` (mapping), `submit_intarna_batchtools.R` (IntaRNA hybridization)
  - Install with: `conda install -c conda-forge r-batchtools` or `install.packages("batchtools")` in R
  - Benefits: Enables distributing chunk processing across multiple cluster nodes for true parallel computing
- **jsonlite** (required for batchtools JSON configuration parsing)
  - Purpose: Parsing JSON configuration files for batchtools job submission
  - Install with: `conda install -c conda-forge r-jsonlite` or `install.packages("jsonlite")` in R
  - Note: Usually installed automatically as a dependency of `batchtools`, but explicitly installing ensures compatibility
- See [BATCHTOOLS_USAGE.md](BATCHTOOLS_USAGE.md) for detailed usage instructions

**Command-line tools:**
- bwa (recommended for alignment)
- samtools
- bedtools
- gffread (for `extract_transcripts_from_genome.py`)
  - Install with: `conda install -c bioconda gffread`
  - Part of GFF utilities from Johns Hopkins University
  - Documentation: https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread
- intarna (optional, for hybridization prediction)
- blockbuster (optional, for block-based merging in `chira_merge.py`)
- clan (optional, alternative aligner to BWA in `chira_map.py`)
- **GNU coreutils** (for parallel sort support)
  - **Linux**: Usually pre-installed (GNU sort is standard)
  - **macOS**: Install via Homebrew: `brew install coreutils`
  - Required for `--parallel` option in sort (version >= 8.6)
  - Code automatically detects and uses parallel sort when available

**Installation commands:**
```bash
# Core packages (required)
pip install biopython bcbiogff pysam mpire requests

# Optional packages (highly recommended for optimal performance)
pip install psutil  # For automatic memory optimization, I/O bottleneck detection, and adaptive buffer sizing (10-50x I/O improvement)
pip install pyliftover  # For coordinate liftover in download_mirbase_gff3.py

# Command-line tools
conda install -c bioconda bwa samtools bedtools gffread intarna

# Optional tools (for specific use cases)
conda install -c bioconda blockbuster  # For block-based merging in chira_merge.py
conda install -c bioconda clan  # Alternative aligner to BWA in chira_map.py

# GNU coreutils (for parallel sort on macOS)
# Linux: Usually pre-installed
# macOS: brew install coreutils

# R packages (for batchtools HPC cluster support)
conda install -c conda-forge r-batchtools r-jsonlite
# Or in R: install.packages(c("batchtools", "jsonlite"))
```

For a complete list of dependencies, see [DEPENDENCIES.md](DEPENDENCIES.md).

For batchtools HPC cluster usage, see [BATCHTOOLS_USAGE.md](BATCHTOOLS_USAGE.md).

---

## Workflow

The ChiRA pipeline consists of five main steps, from raw FASTQ files to final interaction results. The following workflow shows a typical analysis using split-reference mapping (miRNA and target transcriptomes separately).

### Step 1: Prepare Reference Files

Before starting the analysis, prepare your reference files:

**1.1 Download mature miRNA sequences:**

```bash
# Download species-specific mature miRNAs from miRBase and convert "U" to "T"
download_mirbase_mature.py -s hsa -o mature_mirna_hsa.fasta
```

**1.2 Download Ensembl reference files:**

```bash
# Download GTF, and genome FASTA from Ensembl
download_ensembl.py -s homo_sapiens -g 115 -t 115 -o ./ensembl_files
```

**1.3 Download miRBase GFF3 file:**

```bash
# Download species-specific GFF3 file from miRBase (contains mature miRNA coordinates)
download_mirbase_gff3.py -s hsa -o hsa.gff3
```

**Note:** The GFF3 file from miRBase contains mature miRNA coordinates and can be used directly with ChiRA. You don't need to convert it to GTF format. 

**1.4 Prepare target transcriptome (remove miRNAs):**

```bash
# Remove miRNA entries from Ensembl GTF
remove_mirna_hairpin_from_gtf.py -i ensembl_files/Homo_sapiens.GRCh38.115.gtf \
  -o target_transcriptome.gtf

# Extract transcript sequences from genome FASTA using filtered GTF
# This uses gffread to extract transcripts from the genome FASTA based on the filtered GTF
extract_transcripts_from_genome.py -g target_transcriptome.gtf \
  -f ensembl_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -o target_transcriptome.fasta
```

**1.5 Combine miRNA gff3 and filtered GTF (for annotation):**

```bash
# Concatenate miRNA GTF with target GTF
concatenate_gtf.py -m mature_mirna.gtf -t target_transcriptome.gtf \
  -o combined_annotation.gtf
```

**Result:** You now have:

- `mature_mirna_hsa.fasta`: Mature miRNA sequences (from miRBase) with "U" replaced by "T"
- `target_transcriptome.fasta`: Target transcriptome sequences (without miRNAs)
- `combined_annotation.gtf`: Combined annotation file (if combining miRNA and target annotations)
- `GRCh38.primary.assembly.fasta`: Target genome sequences

---

### Step 2: Collapse Reads

Deduplicate reads from FASTQ file:

```bash
chira_collapse.py -i raw_reads.fastq -o collapsed_reads.fasta -u 6
```

**Input:** Raw FASTQ file (quality and adapter trimmed using cutadapt)  
**Output:** FASTA file with unique sequences and read counts

---

### Step 3: Map Reads

Map collapsed reads to reference transcriptomes:

```bash
# Build indices and map (one command)
# -p 8: Use 8 CPU processes (omit to auto-detect all available cores)
# -s fw: Forward strand mapping (default, recommended for CLASH, CLEAR-CLIP, PARIS, SPLASH)
# -a bwa: Use BWA aligner (default)
chira_map.py -i collapsed_reads.fasta -o mapping_output \
  -f1 target_transcriptome.fasta -f2 mature_mirna_hsa.fasta \
  -a bwa -b -p 8 -s fw

# Or use pre-built indices
chira_map.py -i collapsed_reads.fasta -o mapping_output \
  -x1 index1 -x2 index2 -p 8 -s fw

# split a large collapsed_reads.fasta into multiple chuncks for bwa alignment leveraging batchtools to accelerate alignment (currently support LSF-managed HPC cluster)
chira_map.py -i collapsed_reads.fasta -o mapping_output \
   -f1 target_transcriptome.fasta -f2 mature_mirna_hsa.fasta \
   -b -a bwa -p 8 --chunk_fasta 4 \
  --parallel_chunks 4 --use_batchtools --batchtools_queue long \
  --batchtools_cores 8 --batchtools_memory 64GB \
  --batchtools_walltime 48:00 -s fw \
  --batchtools_max_parallel 4 --batchtools_conda_env ~/miniconda3/envs/chira
```

**Output:** `sorted.bam` BAM files, `sorted.bed` file with all alignments, and `unmapped.fasta` file

---

### Step 4: Merge Alignments

Merge overlapping alignments into loci:

```bash
# Using 8 processes for parallel transcript processing
chira_merge.py -b mapping_output/sorted.bed -o merge_output \
  -g combined_annotation.gtf \
  -f1 target_transcriptome.fasta -f2 mature_mirna_hsa.fasta  \
  -ls 1 -p 8
```

**Input:** `sorted.bed` from Step 3  
**Output:** `genomic_exons.bed`,  `merged.bed`,  `segments.bed`, and `transcriptomic_exons.bed`

---

### Step 5: Quantify CRLs

Build Chimeric Read Loci (CRLs) and quantify:

```bash
# Using 8 threads for parallel EM algorithm
chira_quantify.py -b merge_output/segments.bed \
  -m merge_output/merged.bed -o quantify_output \
  --build_crls_too \
  -cs 0.7 -ls 5 -p 8
```

**Input:** `merged.bed`, and `segments.bed` from Step 4  
**Output:** `loci.counts` with TPM values and CRL assignments

---

### Step 6: Extract Interactions

Extract chimeric reads and summarize interactions:

```bash
# Using 8 processes for parallel extraction and hybridization
chira_extract.py -l quantify_output/loci.counts -o extract_output \
  -f1 target_transcriptome.fasta -f2 mature_mirna_hsa.fasta -n sample1 \
  -g combined_annotation.gtf  -p 8 \
  --accessibility  C --ref GRCh38.primary.assembly.fasta \
  --hybridize --summarize
```
**Input:** `loci.counts` from Step 5  
**Output:** 
- `sample1.chimeras.txt`: Individual chimeric reads
- `sample1.singletons.txt`: Non-chimeric reads
- `sample1.interactions.txt`: Summarized interactions (if `-s` used)

---

### Docker Workflow

If using Docker, the workflow is similar but commands are prefixed with `docker run`. Use the pre-built image using the **Dockerfile**:

```bash
# Pull the pre-built image (first time only)
docker pull docker.io/nemat1976/chiraplus:v0.0.1

# Run each step with volume mounts
docker run --rm -v $(pwd)/data:/app/data -v $(pwd)/output:/app/output \
  docker.io/nemat1976/chiraplus:v0.0.1 chira_collapse.py -i data/input.fastq -o output/collapsed.fasta

# ... continue with remaining steps
```

**Note:** Alternatively, you can build the image from the provided Dockerfile:
```bash
docker build -t chira:latest .
```

---

## Tool Documentation

### chira_collapse.py

**Description:** Deduplicates reads from a FASTQ file by collapsing identical sequences and their UMIs (Unique Molecular Identifiers), if available, then outputs a FASTA file with unique sequences and their read counts. This is typically the first step in the ChiRA pipeline to reduce redundancy before mapping.

**Key Features:**
- Handles UMI-tagged reads (optional UMI trimming from 5' end)
- Counts occurrences of each unique sequence
- Fast raw file parsing (no Biopython dependency)
- Accepts gzipped or uncompressed FASTQ

**Required Arguments:**
- `-i, --fastq`: Input FASTQ file (gzipped or uncompressed)
- `-o, --fasta`: Output FASTA file

**Optional Arguments:**
- `-u, --umi_len`: Length of UMI to trim from 5' end of each read and append to the tag id (default: 0, no UMI)

**Output:** FASTA with unique sequences. Header format: `>sequence_id|UMI|read_count` (or `>sequence_id|read_count` if no UMI). Each unique sequence appears once with its total count.

**Usage Example:**
```bash
chira_collapse.py -i input.fastq -o output.fasta -u 12
```

---

### chira_map.py

**Description:** Maps reads to a reference transcriptome using either BWA-MEM (recommended) or CLAN aligner. Performs two-pass mapping (long then short segments) to capture chimeric alignments. Supports strand-specific mapping and can handle split-reference genomes (e.g., miRNA and target transcriptomes separately).

**Key Features:**
- Two alignment algorithms: BWA-MEM (default, recommended for speed) or CLAN (slower, use only if needed)
- Two-pass mapping strategy: first pass for long segments, second pass for short/unmapped segments
- Strand-specificity options: forward, reverse-complement, or both
- Handles chimeric reads with configurable overlap between segments
- Automatic index building option

**Required Arguments:**
- `-i, --query_fasta`: Path to query FASTA file (reads, typically from chira_collapse.py)
- `-o, --outdir`: Output directory for BAM and BED files
- Either `-x1, --index1` (pre-built index) or `-f1, --ref_fasta1` (reference FASTA to build index)

**Optional Arguments:**
- `-x2, --index2`: Second priority index file (for split-reference)
- `-f2, --ref_fasta2`: Second priority reference FASTA file
- `-b, --build`: Build indices from reference FASTA files

**Split Reference:**
A split reference uses two separate reference FASTA files instead of one combined file. This is useful for:
- **Separating different RNA types**: For example, use `ref_fasta1` for target transcript sequences and `ref_fasta2` for miRNA sequences
- **Better chimeric detection**: Helps identify chimeric reads that span between the two reference types (e.g., miRNA-target interactions)
- **Smaller indices**: Each reference can be indexed separately, potentially reducing memory usage

**How to prepare a split reference:**
1. **Create two separate FASTA files:**

   - `ref1.fasta`: First priority reference (e.g., target transcript sequences)
   - `ref2.fasta`: SEcond priority reference (e.g., mature miRNA sequences from [miRBase](https://www.mirbase.org/)). Notes: convert "U" to "T" is necessary.
  
   **Important:** When preparing `ref2.fasta` (target transcripts):
   - **Remove mature miRNA sequences** from the target transcript reference
   - **Remove miRNA hairpin sequences** from the target transcript reference
   - This prevents false-positive chimeric alignments where miRNA sequences align to both references

2. **Prepare GTF/GFF annotation file:**
   - **Include only mature miRNA annotations** (e.g., `3p_mature_mir`, `5p_mature_mir`, `mature_mir`) download .gff3 from miRBase Downloads: https://www.mirbase.org/download/
   - This ensures proper identification of miRNA vs target loci in downstream analysis

3. **Option A - Build indices automatically:**

   ```bash
   chira_map.py -i reads.fasta -o output_dir -f1 ref1.fasta -f2 ref2.fasta -b -a bwa
   ```

   The `-b` flag will automatically build indices for both references.

4. **Option B - Pre-build indices (recommended for repeated use):**

   ```bash
   # Build index1
   bwa index -p index1 ref1.fasta
   
   # Build index2
   bwa index -p index2 ref2.fasta
   
   # Then use pre-built indices
   chira_map.py -i reads.fasta -o output_dir -x1 index1 -x2 index2 -a bwa
   ```

**Note:** When using split reference, reads are mapped to both references separately, then the results are merged. The tool tracks which reference each alignment came from, which is important for downstream analysis in `chira_merge.py` and `chira_extract.py`.

**Key Parameters:**
- `-a, --aligner`: Alignment program (`bwa` or `clan`, default: `bwa`)
  - **Note:** BWA-MEM is significantly faster than CLAN and is recommended for most use cases. CLAN should only be used if specific alignment characteristics are required.
- `-s, --stranded`: Strand specificity (`fw`, `rc`, or `both`, default: `fw`)
  - **`fw`** (forward/transcript strand): Use for stranded libraries where reads align to the transcript strand (recommended for most protocols like CLASH, CLEAR-CLIP, PARIS, SPLASH)
  - **`rc`** (reverse complement): Use for stranded libraries where reads align to the reverse complement strand (rare, check your protocol)
  - **`both`**: Use only for unstranded libraries where strand information is not preserved (rare for interactome protocols)
  - **Why it matters**: Stranded mapping filters out alignments on the wrong strand, reducing false positives and improving chimeric read detection accuracy
- `-l1, --seed_length1`: Seed length for 1st mapping iteration (default: 12)
- `-l2, --seed_length2`: Seed length for 2nd mapping iteration (default: 16)
- `-s1, --align_score1`: Minimum alignment score for 1st iteration (default: 18)
- `-s2, --align_score2`: Minimum alignment score for 2nd iteration (default: 16)
- `-co, --chimeric_overlap`: Max bases between chimeric segments (default: 2)
- `-p, --processes`: Total number of CPU processes/threads to use (default: auto-detects CPU count) for the main job.
  - **Automatic distribution**: The script automatically divides total processes among parallel BWA jobs
    - **Without chunking**: Total processes divided among parallel BWA jobs (2-4 jobs depending on indices)
    - **With chunking**: Total processes divided among parallel chunks, then each chunk uses its allocated processes specified by `--batchtools_cores`
  - **Multi-threading**: Enables parallel processing for `samtools view`, `pysam.merge`, `pysam.sort`, and `sort` commands
  - **Performance**: 2-4x faster for large BAM files and sorting operations
  - **Memory optimization**: Use `--sort_memory` to specify memory per thread (e.g., "2G", "3G"), or install `psutil` for automatic optimization
  - **Recommendation**: Set to total number of CPU cores for optimal performance
- `--sort_memory`: Memory per thread for BAM sorting (e.g. "2G", "3G"). If not set, auto-calculated from available RAM. Total memory = sort_memory × processes.
- `--chunk_fasta`: Split input FASTA into N chunks for parallel processing (optional, recommended for very large files >1GB)
  - **How it works**: 
    - Creates N chunks from the input FASTA file
    - Processes chunks in batches, with parallel execution controlled by `--parallel_chunks` (default: 2)
    - Each chunk processes all BWA jobs sequentially
    - Remaining chunks are processed in subsequent batches automatically
  - **Example 1**: `--chunk_fasta 10 --batchtools_cores 8 --batchtools_memory 64GB --parallel_chunks 2` (default)
    - Creates 10 chunks from FASTA
    - Runs 2 chunks in parallel (8 processes/8GB each)
    - Processes remaining 8 chunks in subsequent batches of 2
  - **Benefits**: Better I/O performance and memory efficiency for large datasets
  - Each chunk is processed independently through all BWA jobs, then results are merged
- `--parallel_chunks`: Number of chunks to process simultaneously when using `--chunk_fasta` (default: 2)
  - **How to set**: Based on available memory and CPU resources
    - **Default: 2** (recommended for most systems)
    - Small systems (<16 cores, <32GB RAM): `--parallel_chunks 1`
    - Medium systems (16-32 cores, 32-64GB RAM): `--parallel_chunks 2` (default)
    - Large systems (>32 cores, >64GB RAM): `--parallel_chunks 2-4`
    - Very large systems (>64 cores, >128GB RAM): `--parallel_chunks 4-8`
  - **Memory consideration**: Each chunk needs ~2-4GB RAM, so total memory ≈ `parallel_chunks × 4GB`
  - **CPU consideration**: Processes per chunk = `--processes / --parallel_chunks`
    - Each chunk should get at least 4 processes for optimal BWA performance
    - Example: `--batchtools_cores 8 --parallel_chunks 4` → 8 processes per chunk (good)
  - **Note**: Only takes effect when `--chunk_fasta` is specified
- `--use_batchtools`: Enable batchtools for HPC cluster job submission (optional, for LSF/SLURM clusters)
  - **Requirements**: R with `batchtools` and `jsonlite` R packages installed
  - See "R packages" section above for installation instructions
  - **Benefits**: Submit chunk jobs to cluster scheduler for true parallel processing across cluster nodes
  - **Path handling**: All file paths are automatically converted to absolute paths for cluster job execution
  - **Batchtools options** (when `--use_batchtools` is used): `--batchtools_queue` (default: long), `--batchtools_cores` (default: 8), `--batchtools_memory`, `--batchtools_walltime` (default: 240:00), `--batchtools_max_parallel`, `--batchtools_conda_env`, `--batchtools_template` (default: lsf_custom.tmpl or "lsf-simple")
  - See [BATCHTOOLS_USAGE.md](BATCHTOOLS_USAGE.md) for detailed usage and examples

**Outputs:**
- `sorted.bam`: Sorted BAM file
- `sorted.bed`: BED file containing all alignments
- `unmapped.fasta`: FASTA file with unmapped reads (optional)

**Usage Example:**
```bash
# Basic usage with 32 total processes (automatically divided among BWA jobs)
chira_map.py -i reads.fasta -o output_dir \
   -f1 ref1.fasta -f2 ref2.fasta -a bwa -s both -p 32

# With manual memory specification for BAM sorting
chira_map.py -i reads.fasta -o output_dir \
   -f1 ref1.fasta -f2 ref2.fasta -a bwa -s both -p 8 --sort_memory 3G

# For very large FASTA files (>1GB), use chunking for better I/O performance
# Creates 10 chunks, processes them in batches (default: 2 chunks at a time)
chira_map.py -i large_reads.fasta -o output_dir \
   -f1 ref1.fasta -f2 ref2.fasta -a bwa -s both -p 32 --chunk_fasta 10

# Example 1: Custom parallel chunks (4 chunks in parallel)
# - Creates 10 chunks from FASTA
# - Runs 4 chunks in parallel (8 processes each, using all 32 processes)
# - Processes remaining 6 chunks in subsequent batches of 4
chira_map.py -i large_reads.fasta -o output_dir \
   -f1 ref1.fasta -f2 ref2.fasta -a bwa -s both -p 8 --chunk_fasta 10 --parallel_chunks 4

# Example 2: Using batchtools for HPC cluster submission (LSF)
# - Splits FASTA into 20 chunks
# - Submits 20 independent LSF jobs (one per chunk)
# - Each job runs on different cluster node with 8 cores and 8GB memory
# - All paths automatically converted to absolute paths for cluster execution
chira_map.py -i large_reads.fasta -o output_dir \
   -f1 ref1.fasta -f2 ref2.fasta -a bwa -s both \
   --chunk_fasta 20 --use_batchtools \
   --batchtools_queue long \
   --batchtools_cores 8 \
   --batchtools_memory 8GB \
   --batchtools_walltime 240:00
```

**Understanding Chunking and Process Management:**

When using `--chunk_fasta`, the script uses a two-stage approach:

1. **FASTA Splitting**: The input FASTA is split into N chunks (as specified by `--chunk_fasta`)
   - Each chunk contains a subset of reads distributed round-robin
   - Chunks are typically 1-3GB each for optimal I/O performance

2. **Parallel Processing**: Chunks are processed in batches, with parallelism controlled by `--parallel_chunks` (default: 2)
   - Number of chunks running simultaneously is set by `--parallel_chunks` (default: 2)
   - Processes per chunk = `--processes / --parallel_chunks`
   - Each chunk should get at least 4 processes for optimal BWA performance
   - If processes per chunk < 4, the number of parallel chunks is automatically reduced
   - Example: 32 processes with `--parallel_chunks 2` → 2 chunks run in parallel (16 processes each)
   - Example: 32 processes with `--parallel_chunks 4` → 4 chunks run in parallel (8 processes each)
   - Example: 8 processes with `--parallel_chunks 2` → 1 chunk runs at a time (8 processes each, auto-reduced)

3. **Batch Processing**: If you have more chunks than `--parallel_chunks`, remaining chunks are processed in subsequent batches
   - Example: 10 chunks with `--parallel_chunks 2` → batches of 2: (2, 2, 2, 2, 2)
   - Example: 10 chunks with `--parallel_chunks 4` → batches of 4: (4, 4, 2)
   - Example: 10 chunks with `--parallel_chunks 1` → batches of 1: (1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

**Why This Approach:**
- **Prevents oversubscription**: Won't try to run more chunks than you have processes for
- **Efficient resource use**: All available processes are utilized without waste
- **Better I/O**: Chunking still provides I/O benefits even with limited parallelism
- **User control**: Adjust `--parallel_chunks` to match your system's memory and CPU resources
- **Automatic fallback**: If processes per chunk < 4, automatically reduces parallel chunks

**Recommendations:**
- For large files (>1GB): Use `--chunk_fasta 10-20` to create manageable chunks
- Set `--processes` to your total CPU core count for optimal performance
- Set `--parallel_chunks` based on your system:
  - **Small systems** (<16 cores, <32GB RAM): `--parallel_chunks 1`
  - **Medium systems** (16-32 cores, 32-64GB RAM): `--parallel_chunks 2` (default)
  - **Large systems** (>32 cores, >64GB RAM): `--parallel_chunks 2-4`
  - **Very large systems** (>64 cores, >128GB RAM): `--parallel_chunks 4-8`
- Ensure each chunk gets at least 4 processes: `--processes / --parallel_chunks >= 4`

---

### chira_merge.py

**Description:** Merges overlapping aligned positions to define read-concentrated loci (RCLs). If a GTF annotation file is provided, transcriptomic coordinates are converted to genomic coordinates. Segments reads into aligned portions and merges overlapping segments using configurable overlap thresholds.

**Key Features:**
- Converts transcriptomic to genomic coordinates (if GTF provided)
- Segments reads based on alignment parts
- Merges overlapping alignments into loci
- Two merging modes: segment-based or block-based (using Blockbuster)
- Filters by minimum locus size

**Required Arguments:**
- `-b, --bed`: Input BED file with alignments (e.g. from chira_map)
- `-o, --outdir`: Output directory for merged BED and segments

**Optional Arguments:**
- `-g, --gtf`: Annotation GTF file for coordinate conversion
- `-f1, --ref_fasta1`: First priority reference FASTA
- `-f2, --ref_fasta2`: Second priority reference FASTA

**Key Parameters:**
- `-ao, --alignment_overlap`: Minimum fraction overlap among BED entries to merge [0–1.0] (default: 0.7)
- `-so, --segment_overlap`: Merge read positions with greater than this overlap into a segment (default: 0.7)
- `-lt, --length_threshold`: Minimum alignment length as fraction of longest [0.8–1.0] (default: 0.9)
- `-co, --chimeric_overlap`: Max bases between chimeric segments (default: 2)
- `-c, --chimeric_only`: Consider only chimeric reads for merging
- `-ls, --min_locus_size`: Minimum number of alignments per merged locus (default: 1)
- `-bb, --block_based`: Use Blockbuster for block-based merging
- Blockbuster parameters: `-d, --distance` (default: 30), `-mc, --min_cluster_height` (default: 10), `-mb, --min_block_height` (default: 10), `-sc, --scale` (default: 0.1)
- `-p, --processes`: Number of parallel processes for transcript processing (default: auto-detect CPU count)
  - Chunk-based multiprocessing; efficiently handles very large datasets (e.g., human genome with 387K+ transcripts). 4–8x faster for datasets with many transcripts.

**Outputs:**
- `segments.bed`: A BED file with reads categorized into segments
- `merged.bed`: A tabular file with merged alignments
  - Column 4: All alignments merged into that location
  - Column 5: Number of reads supporting the locus

**Usage Example:**
```bash
# Basic usage with 8 processes (auto-detects CPU count if not specified)
chira_merge.py -b mapped.bed -o output_dir -g annotation.gtf -f1 ref1.fasta -ao 0.7 -so 0.7 -p 8

# Auto-detect CPU count (default behavior)
chira_merge.py -b mapped.bed -o output_dir -g annotation.gtf -f1 ref1.fasta -ao 0.7 -so 0.7
```

---

### chira_quantify.py

**Description:** Creates Chimeric Read Loci (CRLs) from merged BED files and quantifies them using an Expectation-Maximization (EM) algorithm to handle multi-mapping reads. Calculates TPM (Transcripts Per Million) values for each CRL.

**Key Features:**
- Builds CRLs by grouping loci that share significant read overlap
- EM algorithm for resolving multi-mapping reads
- TPM normalization for expression quantification
- Configurable CRL building thresholds

**Required Arguments:**
- `-b, --bed`: Input BED file (e.g. segments.bed from chira_merge)
- `-m, --merged_bed`: Input merged BED file (e.g. merged.bed from chira_merge)
- `-o, --outdir`: Output directory (writes loci.counts)

**Key Parameters:**
- `-cs, --crl_share`: Minimum fraction of locus reads that must overlap with all CRL loci to merge (default: 0.7)
- `-ls, --min_locus_size`: Minimum reads per locus to participate in CRL creation (default: 10)
- `-e, --em_threshold`: Max difference in transcript expression between consecutive EM iterations for convergence (default: 0.00001)
- `-crl, --build_crls_too`: Build CRLs in addition to quantification
- `-p, --processes`: Number of parallel processes for EM (default: 0 = use all available cores)
  - Uses MPIRE WorkerPool; 50–90% memory reduction, 2–3x faster startup. Falls back to sequential for very small datasets.
- `--use_sqldb`: Use SQLite backend for very large inputs (50GB+). Disk-backed; lower RAM (~8–16GB) vs in-memory (200–500GB).

**Outputs:**
- `loci.counts`: Tabular file containing reads, their CRLs, and TPM values
  - Each line represents a read-CRL association with TPM quantification
  - Used as input for `chira_extract.py`

**Usage Example:**
```bash
# Basic usage with 8 processes
chira_quantify.py -b segments.bed -m merged.bed -o output_dir -cs 0.7 -ls 10 -p 8

# Use all available CPU cores automatically
chira_quantify.py -b segments.bed -m merged.bed -o output_dir -cs 0.7 -ls 10 -p 0
```

---

### chira_extract.py

**Description:** Extracts the best chimeric alignments for each read and optionally performs RNA-RNA hybridization prediction using IntaRNA. Summarizes interactions at the locus level and identifies miRNA-target pairs. This is the final step that produces the interaction results.

**Key Features:**
- Extracts chimeric reads with best scoring alignments
- Optional RNA-RNA hybridization prediction (IntaRNA)
- Interaction summarization at locus level
- Identifies miRNA vs target loci based on annotation
- Customizable output file names with sample name prefix

**Required Arguments:**
- `-l, --loci`: Input BED file with alignments (e.g. loci.counts from chira_quantify)
- `-o, --out`: Output directory path
- `-n, --sample_name`: Sample name prefix for output files
- `-f1, --ref_fasta1`: First priority reference FASTA file
- `-g, --gtf`: Annotation GTF file (including miRBase GFF3 for mature miRNA)
- `-f2, --ref_fasta2`: Second priority reference FASTA (e.g. miRNA)
- `-f, --ref`: Reference genomic FASTA (for IntaRNA accessibility)

**Key Parameters:**
- `-p, --processes`: Number of processes for extraction and hybridization prep/finish (default: 1)
- `-tc, --tpm_cutoff`: Discard transcripts below this TPM percentile (default: 0)
- `-sc, --score_cutoff`: Discard hybrids below this score [0–1.0] (default: 0.0)
- `-co, --chimeric_overlap`: Max bases between chimeric segments (default: 2)
- `-s, --summarize`: Summarize interactions at locus level
- `-z, --gzip`: Compress output files (chimeras and singletons) with gzip

**IntaRNA / hybridization (requires `-r, --hybridize`):**
- `-r, --hybridize`: Run IntaRNA to hybridize predicted chimeras
- `-ns, --no_seed`: Do not enforce seed interactions
- `-acc, --accessibility`: IntaRNA accessibility: `C` (compute) or `N` (not) (default: `N`)
- `-m, --intarna_mode`: IntaRNA mode: `H` (heuristic), `M` (exact), `S` (seed-only) (default: `H`)
- `-t, --temperature`: IntaRNA temperature in Celsius (default: 37)
- `-sbp, --seed_bp`: IntaRNA --seedBP: base pairs in seed region (default: 5)
- `-smpu, --seed_min_pu`: IntaRNA --seedMinPu: minimal unpaired probability in seed (default: 0)
- `-accw, --acc_width`: IntaRNA --accW: sliding window size for accessibility (default: 150)

**Batchtools (HPC cluster, requires `--hybridize`):**
- `--use_batchtools`: Submit IntaRNA jobs via R batchtools. Requires R with batchtools and IntaRNA on cluster PATH. Prepare writes per-chunk `query.fa`, `target.fa`, and `loci_seqs.pkl`; finish loads `loci_seqs.pkl` and `result.csv` per chunk to write chimeras. See [BATCHTOOLS_USAGE.md](BATCHTOOLS_USAGE.md).
- If the main job times out after submitting IntaRNA jobs, run **merge_intarna_into_chimeras.py** (required: `--outdir`, `--sample-name`, `--n-chunks`) after all IntaRNA jobs complete to merge results and produce chimeras.txt, singletons.txt, and interactions.txt.
- `--remove_intermediate`: Remove intermediate files after success (loci.fa.<n>, loci.bed.<n>, batchtools_work/, *.chimeras.<n>, chimeras-r.<n>, singletons.<n>). Default: keep them.
- `--batchtools_registry`: Registry directory (default: `<outdir>/batchtools_work/registry`).
- `--batchtools_template`: LSF template path or `"lsf-simple"` (default: lsf_custom.tmpl if present).
- `--batchtools_queue`: LSF queue (default: long).
- `--batchtools_cores`: Cores per LSF job (default: 8). Within each job that many IntaRNA pairs run in parallel (faster chunks).
- `--batchtools_memory`: Total memory per job (e.g. 8GB, 64GB). Converted to per-core for LSF.
- `--batchtools_walltime`: Walltime per job (default: 48:00).
- `--batchtools_conda_env`: Conda environment for cluster jobs (optional).
- `--batchtools_max_parallel`: Max concurrent batchtools jobs (default: None = all chunks at once for minimum walltime).

**Command-line settings to finish sooner / shorten walltime:**

- **Skip steps you don't need:** Omit `-r` / `--hybridize` to skip IntaRNA; omit `-s` / `--summarize` if you don't need the interactions file.
- **Submit all IntaRNA jobs at once:** Leave `--batchtools_max_parallel` unset (default) so all chunk jobs are submitted together and run in parallel. Set it only to limit concurrent jobs (e.g. for queue policy).
- **More chunks and parallelism:** Use higher `-p` (e.g. `-p 16` or `-p 32`) so you have more IntaRNA chunk jobs; with default max_parallel they all run concurrently. Use `--batchtools_cores 8` (or higher) so each chunk runs many IntaRNA pairs in parallel within the job; or `--batchtools_cores 1` to run more jobs at once (one core per job).
- **Less work (filter earlier):** Raise `-tc` or `-sc` to drop low-TPM/low-score alignments; fewer chimeras and less IntaRNA.
- **Faster IntaRNA:** Keep `-acc N` and `-m H` (defaults). Merge step uses extra sort threads when available to shorten merge time.

**Example — minimize walltime with hybridization:**
```bash
chira_extract.py -l loci.counts -o out -n mysample -f1 ref.fa -g ref.gtf -f ref_genome.fa \
  -p 32 -tc 0.2 -sc 0.1 -r --use_batchtools \
  --batchtools_cores 8 --batchtools_memory 8GB 
```
(Default: all 32 chunk jobs submitted at once. Use `-tc` / `-sc` only if you accept dropping some low-TPM/low-score chimeras.)

**Example — fastest run without IntaRNA:**
```bash
chira_extract.py -l loci.counts -o out -n mysample -f1 ref.fa -g ref.gtf -p 8
```

**Outputs:**

**Note:** If `--gzip` is specified, output files will have `.gz` extension (e.g., `{sample_name}.chimeras.txt.gz`, `{sample_name}.singletons.txt.gz`). Compression is applied only to final merged files, not intermediate files, for optimal performance.

**1. `{sample_name}.chimeras.txt`** (or `.txt.gz` if `--gzip` is used) - Tabular file with chimeric read information in an extended BED format (tab-separated).

**Header line (35 columns when using hybridization):**
```
read_id	transcript_id_1	transcript_id_2	gene_id_1	gene_id_2	gene_symbol_1	gene_symbol_2	annotation_region_1	annotation_region_2	transcript_start_1	transcript_end_1	transcript_strand_1	transcript_length_1	transcript_start_2	transcript_end_2	transcript_strand_2	transcript_length_2	read_alignment_info	genomic_coordinates_1	genomic_coordinates_2	locus_id_1	locus_id_2	crl_group_id_1	crl_group_id_2	tpm_1	tpm_2	alignment_score_1	alignment_score_2	combined_alignment_score	hybridized_sequences	hybridization_structure	hybridization_positions	hybridization_mfe_kcal_mol	mirna_read_position	hybridized_subsequences
```

**Column descriptions:**
- `read_id`: Read identifier (from collapsed FASTQ)
- `transcript_id_1`, `transcript_id_2`: Transcript IDs for locus 1 and locus 2
- `gene_id_1`, `gene_id_2`: Gene IDs for locus 1 and locus 2
- `gene_symbol_1`, `gene_symbol_2`: Gene symbols for locus 1 and locus 2
- `annotation_region_1`, `annotation_region_2`: Annotation regions (e.g. `3p_mature_mir`, `5p_mature_mir`, `mature_mir` for miRNA; gene/exon types for targets)
- `transcript_start_1`, `transcript_end_1`, `transcript_strand_1`, `transcript_length_1`: Transcriptomic coordinates and alignment length for locus 1
- `transcript_start_2`, `transcript_end_2`, `transcript_strand_2`, `transcript_length_2`: Transcriptomic coordinates and alignment length for locus 2
- `read_alignment_info`: Format `arm1_start,arm1_end,arm2_start,arm2_end,read_length`. Use to determine actual read orientation (compare `arm1_start` vs `arm2_start`: if `arm1_start < arm2_start` then 5' locus1 → locus2 3'; if `arm2_start < arm1_start` then 5' locus2 → locus1 3')
- `genomic_coordinates_1`, `genomic_coordinates_2`: Genomic coordinates (if GTF provided)
- `locus_id_1`, `locus_id_2`: Locus identifiers (`chr:start:end:strand`)
- `crl_group_id_1`, `crl_group_id_2`: CRL group IDs
- `tpm_1`, `tpm_2`: TPM (Transcripts Per Million) for each locus
- `alignment_score_1`, `alignment_score_2`, `combined_alignment_score`: Alignment scores (combined = score1 × score2)
- `hybridized_sequences`: `sequence1&sequence2` or `NA`
- `hybridization_structure`: Dot-bracket RNA-RNA structure or `NA`
- `hybridization_positions`: Hybridization position information or `NA`
- `hybridization_mfe_kcal_mol`: Minimum free energy (kcal/mol) or `NA`
- `mirna_read_position`: `miRNA_first` (5' miRNA → target 3'), `miRNA_last` (5' target → miRNA 3'), or `NA`
- `hybridized_subsequences`: Hybridized subsequence segments from IntaRNA (format: `subseq1&subseq2`, or `NA`)

**Note on chimeric read orientation:**
The tool detects chimeric reads in both orientations:
- **5' miRNA → target 3'**: miRNA at 5' end, target at 3' end
- **5' target → miRNA 3'**: Target at 5' end, miRNA at 3' end

For **split reference** (when `-f2` is provided), the output is standardized so that:
- `locus1` always corresponds to the first reference (typically miRNA from `ref_fasta1`)
- `locus2` always corresponds to the second reference (typically target from `ref_fasta2`)

**To determine the actual read orientation**, check the `read_alignment_info` column:
- It contains: `arm1_start,arm1_end,arm2_start,arm2_end,read_length`
- Compare the start positions:
  - If `arm1_start < arm2_start`: The read is in **5' locus1 → locus2 3'** orientation
  - If `arm2_start < arm1_start`: The read is in **5' locus2 → locus1 3'** orientation (reoriented in output for split reference)

**Example:**
- If `read_alignment_info = "1,20,21,40,50"`: arm1 starts at position 1, arm2 at 21 → Read is 5' locus1 → locus2 3'
- If `read_alignment_info = "21,40,1,20,50"`: arm2 starts at position 1, arm1 at 21 → Read is 5' locus2 → locus1 3' (reoriented in output for split reference)

In the **interactions file**, both orientations are merged into a single entry; the chimeras file column `read_alignment_info` preserves the actual orientation.

**2. `{sample_name}.singletons.txt`** (or `.txt.gz` if `--gzip` is used) - Tabular file with singleton reads (non-chimeric alignments).

Header columns (15 total, tab-separated):
- `read_id`: Read identifier
- `transcript_id`: Transcript ID
- `gene_id`: Gene ID
- `gene_symbol`: Gene symbol
- `annotation_region`: Annotation region type
- `transcript_start`, `transcript_end`, `transcript_strand`: Transcriptomic coordinates
- `transcript_length`: Alignment length
- `read_alignment_info`: Read alignment information
- `genomic_coordinates`: Genomic coordinates (if GTF provided)
- `locus_id`: Locus identifier (format `chr:start:end:strand`)
- `crl_group_id`: CRL group ID
- `tpm`: TPM value
- `alignment_score`: Alignment score

**3. `{sample_name}.interactions.txt`** - Tabular file with detected interactions (if `--summarize` used). This file is always uncompressed for compatibility with downstream analysis tools.

Header columns (28 total, tab-separated):
- `supporting_read_count`: Number of reads supporting this interaction
- `locus_1_chromosome`, `locus_1_start`, `locus_1_end`, `locus_1_strand`: Genomic coordinates for locus 1
- `locus_2_chromosome`, `locus_2_start`, `locus_2_end`, `locus_2_strand`: Genomic coordinates for locus 2
- `locus_1_sequence`, `locus_2_sequence`: Sequences from each locus (or `NA` if not hybridized)
- `hybridization_structure_dotbracket`: RNA-RNA hybridization structure in dot-bracket notation (format: `structure1&structure2`, or `NA`)
- `hybridization_mfe_kcal_mol`: Minimum free energy of hybridization in kcal/mol (or `NA`)
- `hybridized_sequence_segments`: Hybridized sequence segments (format: `seq1\tseq2`, or `NA\tNA`)
- `hybridization_start_positions`: Start positions of hybridization within sequences (format: `pos1&pos2`, or `NA`)
- `hybridization_genomic_coordinates`: Genomic coordinates of hybridization region (format: `refid1:ref_start1-ref_end1:ref_strand1&refid2:ref_start2-ref_end2:ref_strand2`, or `NA`)
- `tpm_locus_1`, `tpm_locus_2`, `tpm_combined`: TPM values per locus and combined
- `alignment_score_locus_1`, `alignment_score_locus_2`, `combined_alignment_score`: Alignment scores per locus and combined (score1 × score2)
- `annotation_region_locus_1`, `annotation_region_locus_2`: Annotation regions (semicolon-separated if multiple; use to identify miRNA vs target)
- `reference_transcript_id_1`, `reference_transcript_id_2`: Reference transcript IDs (semicolon-separated if multiple)
- `gene_name_1`, `gene_name_2`: Gene symbols for locus 1 and locus 2 (semicolon-separated if multiple)

**Note:** The interactions file may include comment lines explaining how to identify miRNA vs target loci. miRNA annotations typically include: `miRNA`, `3p_mature_mir`, `5p_mature_mir`, `mature_mir`.

**Usage Example:**
```bash
# Basic usage
chira_extract.py -l loci.txt -o output_dir -f1 ref1.fasta -n sample1 \
  -g annotation.gtf -r -s -tc 0.1 -sc 0.5

# With gzip compression (recommended for large files)
chira_extract.py -l loci.txt -o output_dir -f1 ref1.fasta -n sample1 \
  -g annotation.gtf -r -s -tc 0.1 -sc 0.5 --gzip
```

---

## Utility Scripts

The following utility scripts are provided to help prepare reference files and annotations for ChiRA analysis. The scripts `process_chunk_batchtools.py` and `process_intarna_chunk_batchtools.py` are invoked internally by the batchtools R workflows and are not intended to be run directly by users.

### merge_intarna_into_chimeras.py

**Description:** Standalone script to merge IntaRNA results into chimeras and produce final sample-level chimeras, singletons, and interactions files. Use this when `chira_extract.py --hybridize --use_batchtools` times out after submitting IntaRNA jobs but before the finish phase (e.g. LSF walltime limit). After all IntaRNA chunk jobs have completed, run this script to: (1) merge each chunk’s `loci_seqs.pkl` and `result.csv` into `*.chimeras-r.<n>`, (2) merge `*.chimeras-r.<n>` → `{sample_name}.chimeras.txt`, (3) merge `*.singletons.<n>` → `{sample_name}.singletons.txt`, (4) generate `{sample_name}.interactions.txt`.

**Required Arguments:** `--outdir`, `--sample-name`, `--n-chunks` (number of chunks, 0..N-1).

**Optional Arguments:** `--chunk-root` (default: `outdir/batchtools_work`), `--remove-input` (remove `*.chimeras.<n>` after writing `*.chimeras-r.<n>`), `--buffer-size`.

**Usage Example:**
```bash
python merge_intarna_into_chimeras.py --outdir /path/to/extract_output --sample-name mysample --n-chunks 32
```

**Dependencies:** Imports from `chira_extract` (parse_intarna_csv, _merge_hybrid_into_chimeras, merge_files, write_interaction_summary). Chunk directories must contain `loci_seqs.pkl` and `result.csv`; `*.chimeras.<n>` and `*.singletons.<n>` must exist under `outdir`.


### download_ensembl.py

**Description:** Downloads GTF and genome FASTA files from Ensembl for a given species and release versions.

**Key Features:**
- Downloads primary assembly genome FASTA (not toplevel)
- Auto-detects assembly names or accepts explicit assembly parameter
- Supports HTTP and FTP with automatic fallback
- Automatically decompresses gzipped files

**Required Arguments:**
- `-s, --species`: Species (e.g., homo_sapiens, mus_musculus, bos_taurus)
- `-g, --genome-version`: Ensembl release version for genome (e.g., 110)
- `-t, --gtf-version`: Ensembl release version for GTF annotation (e.g., 110)
- `-o, --output-dir`: Output directory for downloaded files

**Optional Arguments:**
- `-a, --assembly`: Genome assembly (e.g., GRCh38). Auto-detect if not set.
- `--keep-compressed`: Keep compressed files after decompression (default: remove)
- `--no-decompress`: Do not decompress gzipped files (default: decompress)
- `--timeout`: Download timeout in seconds (default: 60)

**Usage Example:**
```bash
download_ensembl.py -s homo_sapiens -g 110 -t 110 -o ./ensembl_files
```

---

### download_mirbase_mature.py

**Description:** Downloads species-specific mature miRNA sequences from miRBase.

**Key Features:**
- Downloads from specific miRBase version or CURRENT
- Extracts sequences by species code
- Handles gzipped files automatically
- **Automatically converts U (uracil) to T (thymine)** in sequences for ChiRA compatibility (ChiRA expects DNA sequences)

**Required Arguments:**
- `-s, --species`: Species code (e.g., hsa=human, mmu=mouse, bta=bovine, rno=rat)
- `-o, --output`: Output FASTA file for species-specific mature miRNAs

**Optional Arguments:**
- `--mirbase-version`: miRBase version (e.g., "22.1"). Default: CURRENT
- `--keep-full`: Keep full mature.fa after extraction (default: remove)
- `--timeout`: Download timeout in seconds (default: 30)

**Usage Example:**
```bash
download_mirbase_mature.py -s hsa -o mature_mirna_hsa.fasta
```

---

### download_mirbase_gff3.py

**Description:** Downloads species-specific GFF3 file from miRBase containing chromosomal coordinates of microRNAs. Supports coordinate liftover between genome assemblies and chromosome name mapping.

**Key Features:**
- Downloads current or version-specific GFF3 files
- Contains both miRNA_primary_transcript (hairpin precursors) and miRNA (mature sequences)
- **Coordinate liftover**: Convert coordinates between genome assemblies (e.g., hg19 → hg38) using pyliftover
- **Chromosome name mapping**: Rename chromosomes based on a mapping file
- Processing order: Download → Liftover → Rename chromosomes
- Can be used directly with ChiRA (no GTF conversion needed)

**Required Arguments:**
- `-s, --species`: Species code (e.g., hsa=human, mmu=mouse, bta=bovine, rno=rat)
- `-o, --output`: Output GFF3 file for species-specific microRNA annotations

**Optional Arguments:**
- `--mirbase-version`: miRBase version (e.g., "21"). Default: CURRENT
- `--timeout`: Download timeout in seconds (default: 60)
- `-m, --chromosome-mapping`: Tab-separated file: `gff3_chromosome_name<tab>target_chromosome_name`. Applied after liftover (if performed).
- `--source-genome`: Source genome assembly for liftover (e.g., hg19, mm9). Required if `--chain-file` is given.
- `--target-genome`: Target genome assembly for liftover (e.g., hg38, mm10). Required if `--chain-file` is given.
- `--chain-file`: Path to chain file for coordinate liftover (e.g. from UCSC). Requires `--source-genome` and `--target-genome`. Needs `pyliftover` (`pip install pyliftover`).

**Usage Examples:**
```bash
# Download current version
download_mirbase_gff3.py -s hsa -o hsa.gff3

# Download specific version with chromosome mapping
download_mirbase_gff3.py -s hsa -o hsa.gff3 --mirbase-version 21 -m chr_mapping.txt

# Download with coordinate liftover (hg19 to hg38)
download_mirbase_gff3.py -s hsa -o hsa_hg38.gff3 \
  --source-genome hg19 --target-genome hg38 \
  --chain-file hg19ToHg38.over.chain

# Download with both liftover and chromosome mapping
download_mirbase_gff3.py -s hsa -o hsa_processed.gff3 \
  --source-genome hg19 --target-genome hg38 \
  --chain-file hg19ToHg38.over.chain \
  -m chr_mapping.txt
```

**Note:** The GFF3 file from miRBase contains mature miRNA coordinates and can be used directly with ChiRA. You don't need to convert it to GTF format unless you want to combine it with Ensembl annotations.

**Coordinate Liftover:**
- Liftover converts coordinates from one genome assembly to another (e.g., GRCh37/hg19 to GRCh38/hg38)
- Chain files are available from UCSC Genome Browser for common assembly conversions
- Liftover is performed first, then chromosome renaming (if provided)
- Features that cannot be lifted over will retain their original coordinates

---

### remove_mirna_hairpin_from_gtf.py

**Description:** Removes all microRNA entries from an Ensembl GTF file. Used for preparing target-only transcriptome annotations.

**Key Features:**
- Identifies miRNA entries by feature type, biotype, and optional regex pattern
- Flexible pattern matching for custom miRNA identification
- Preserves comment lines (optional removal)

**Required Arguments:**
- `-i, --input`: Input GTF file
- `-o, --output`: Output GTF file (without microRNA entries)

**Optional Arguments:**
- `-p, --pattern`: Regex for matching miRNA in GTF attributes. If not set, only feature type and biotype are used. Example: `'gene_name\s+"[^"]*(?:[Mm][Ii][Rr][_-]|[^"]*[-_][Mm][Ii][Rr][_-])'`
- `--remove-comments`: Remove comment lines from output (default: keep comments)

**Usage Example:**
```bash
remove_mirna_hairpin_from_gtf.py -i annotation.gtf -o annotation_no_mirna.gtf
```

---

### extract_transcripts_from_genome.py

**Description:** Extracts transcript FASTA sequences from a genome FASTA file using gffread based on a filtered GTF file.

This script uses gffread (from [GFF utilities](https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)) to extract transcript sequences directly from the genome FASTA file based on transcript features in a filtered GTF file (e.g., output from `remove_mirna_hairpin_from_gtf.py`). This approach is more accurate than filtering pre-extracted transcript sequences, as it extracts sequences directly from the genome coordinates.

**Key Features:**
- Uses gffread to extract transcript sequences from genome FASTA
- Works with filtered GTF files (e.g., miRNA-removed GTF from `remove_mirna_hairpin_from_gtf.py`)
- Produces transcriptome FASTA without miRNA sequences (if GTF was filtered)
- Automatically validates that gffread is available

**Required Arguments:**
- `-g, --gtf`: Filtered GTF file (e.g., from `remove_mirna_hairpin_from_gtf.py`)
- `-f, --genome-fasta`: Genome FASTA (e.g., primary assembly from Ensembl)
- `-o, --output`: Output transcript FASTA file

**Dependencies:** Requires `gffread` (e.g. `conda install -c bioconda gffread`). gffread is part of the GFF utilities package from Johns Hopkins University.

**Usage Example:**
```bash
# Extract transcripts from genome using filtered GTF (without miRNAs)
extract_transcripts_from_genome.py -g target_transcriptome.gtf \
  -f Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -o target_transcriptome.fasta
```

**Note:** This script replaces the previous `remove_mirna_hairpin_from_fasta.py` approach. Instead of filtering pre-extracted transcript sequences, it extracts sequences directly from the genome FASTA based on the filtered GTF coordinates, ensuring accuracy and completeness.

---

### concatenate_gtf.py

**Description:** Concatenates mature miRNA GTF/GFF3 file with target transcriptome GTF file, removing comment lines from the miRNA GTF.

**Key Features:**
- Combines miRNA GFF3 and target GTF files for split-reference analysis
- Removes comment lines from miRNA GTF
- Optionally removes comment lines from target GTF
- miRBase GFF3 format can be used directly with ChiRA; this script accepts GTF or GFF3.

**Required Arguments:**
- `-m, --mirna-gtf`: Mature miRNA GTF/GFF3 (e.g. from `download_mirbase_gff3.py`). GFF3 works with ChiRA.
- `-t, --target-gtf`: Target transcriptome GTF (e.g. from `remove_mirna_hairpin_from_gtf.py`)
- `-o, --output`: Output combined GTF file

**Optional Arguments:**
- `--remove-target-comments`: Remove comment lines from target GTF (default: keep)

**Usage Example:**
```bash
concatenate_gtf.py -m mature_mirna.gtf -t target_no_mirna.gtf -o combined.gtf
```

---
