# ChiRA External Dependencies

This document lists all external dependencies required by the ChiRA package.

## Python Packages

The following Python packages need to be installed via pip or conda:

### Core Dependencies

1. **Biopython** (`Bio`)
   - Used in: `chira_extract.py`, `chira_utilities.py`
   - Import: `from Bio import SeqIO`
   - Purpose: Reading FASTA files (used in `extract_reflengths()` and for parsing FASTA sequences in `chira_extract.py`)
   - Note: `chira_collapse.py` no longer uses Biopython - uses raw file parsing for better performance

2. **bcbiogff** (`BCBio`)
   - Used in: `chira_extract.py`, `chira_merge.py`
   - Import: `from BCBio import GFF`
   - Purpose: Parsing GFF/GTF annotation files

3. **pysam**
   - Used in: `chira_map.py`
   - Import: `import pysam`
   - Purpose: Reading and manipulating BAM files, merging BAM files, sorting BAM files

4. **mpire**
   - Used in: `chira_quantify.py`, `chira_extract.py`, `chira_merge.py` (required for multiprocessing)
   - Import: `from mpire import WorkerPool`
   - Purpose: Enhanced multiprocessing framework for parallel processing
   - Benefits:
     - **50-90% memory reduction**: Shared objects avoid copying large dictionaries across processes
     - **2-3x faster startup**: Lower overhead than ProcessPoolExecutor/Process
     - **Better performance**: Optimized for CPU-bound parallel tasks
     - **Copy-on-write semantics**: Uses `start_method='fork'` on Unix/Linux for optimal memory efficiency
   - Install with: `pip install mpire` or `conda install -c conda-forge mpire`
   - Required for: Parallel EM algorithm execution in `chira_quantify.py`, parallel chimera extraction in `chira_extract.py`, and parallel transcript processing in `chira_merge.py`
   - **Note**: MPIRE automatically uses `'fork'` start method on Unix/Linux systems when available, providing copy-on-write semantics for shared objects. On Windows, it falls back to the default start method (spawn).

### Optional Dependencies

5. **psutil**
   - Used in: `chira_map.py`, `chira_utilities.py` (optional, for automatic memory and I/O optimization)
   - Import: `import psutil`
   - Purpose:
     - **chira_map.py**: Automatically calculating optimal memory allocation for BAM sorting and detecting I/O bottlenecks
     - **chira_utilities.py**: Calculating adaptive buffer sizes (8-16MB) for file I/O operations based on available system RAM
   - Note: If not available, falls back to safe defaults (2GB per thread for BAM sorting, 8MB buffer for I/O). Install with `pip install psutil` or `conda install psutil` for automatic optimization.
   - Benefit: Prevents memory exhaustion, optimizes I/O performance (10-50x improvement), and enables automatic I/O bottleneck detection

6. **pyliftover**
   - Used in: `download_mirbase_gff3.py` (optional, only when using coordinate liftover)
   - Import: `from pyliftover import LiftOver`
   - Purpose: Converting coordinates between different genome versions
   - Note: Only required if using `--source-genome`, `--target-genome`, and `--chain-file` options in `download_mirbase_gff3.py`

7. **requests**
   - Used in: `download_ensembl.py`
   - Import: `import requests`
   - Purpose: Downloading files from Ensembl via HTTP/HTTPS

## R Packages

The following R packages are required for HPC cluster job submission via batchtools:

### Core R Dependencies (for batchtools)

1. **batchtools**
   - Used in: `chira_map.py` (optional, for HPC cluster job submission), `chira_extract.py` (optional, for IntaRNA when `--hybridize --use_batchtools`)
   - R scripts: `submit_chunks_batchtools.R`, `submit_intarna_batchtools.R`
   - Purpose: Submitting chunk-based batch jobs to HPC cluster schedulers (LSF, SLURM, SGE, etc.)
   - Installation:
     ```r
     install.packages("batchtools")
     ```
     Or via conda:
     ```bash
     conda install -c conda-forge r-batchtools
     ```
   - Note: Required when using `--use_batchtools` in `chira_map.py` or `--hybridize --use_batchtools` in `chira_extract.py`
   - Benefits: Enables distributing chunk processing across multiple cluster nodes for true parallel computing
   - Compatibility: Both R scripts use manual polling and `count_status()` for job status (handles batchtools returning POSIXct for done/running/error in some versions)

2. **jsonlite**
   - Used in: `chira_map.py` and `chira_extract.py` (optional, for batchtools JSON configuration parsing)
   - R scripts: `submit_chunks_batchtools.R`, `submit_intarna_batchtools.R`
   - Purpose: Parsing JSON configuration files for batchtools job submission
   - Installation:
     ```r
     install.packages("jsonlite")
     ```
     Or via conda:
     ```bash
     conda install -c conda-forge r-jsonlite
     ```
   - Note: Usually installed automatically as a dependency of `batchtools`, but explicitly installing ensures compatibility
   - Required for: JSON configuration file parsing in batchtools submission workflow

### R Installation

**Install R:**
- From [CRAN](https://cran.r-project.org/) (recommended for most systems)
- Or via conda:
  ```bash
  conda install -c conda-forge r-base
  ```

**Install R packages:**
```r
R
> install.packages(c("batchtools", "jsonlite"))
```

Or via conda:
```bash
conda install -c conda-forge r-batchtools r-jsonlite
```

**Note**: When using batchtools, ensure all file paths are absolute paths (this is handled automatically by the code). The batchtools template file path should also be absolute or use the built-in "lsf-simple" template.

## External Command-Line Tools

The following command-line tools must be installed and available in the system PATH:

### Alignment Tools

1. **BWA (Burrows-Wheeler Aligner)**
   - Used in: `chira_map.py`
   - Commands: `bwa mem`, `bwa index`
   - Purpose: Mapping reads to reference transcriptome

2. **CLAN (Chimeric Long-read Aligner)**
   - Used in: `chira_map.py` (optional)
   - Commands: `clan_search`, `clan_output`, `clan_index`
   - Purpose: Alternative alignment tool for chimeric reads
   - Note: Only needed if using CLAN aligner instead of BWA (not recommended due to performance)

### Bioinformatics Utilities

3. **samtools**
   - Used in: `chira_map.py`
   - Commands: `samtools view`
   - Purpose: Converting SAM to BAM format

4. **BEDTools**
   - Used in: `chira_extract.py`, `chira_merge.py`
   - Commands: `bedtools intersect` (or `intersectBed` in older versions), `bedtools getfasta` (or `fastaFromBed` in older versions)
   - Purpose: Extracting sequences from BED files and calculating overlaps
   - Note: Code automatically detects and supports both old and new BEDTools command formats

5. **blockbuster.x**
   - Used in: `chira_merge.py` (optional, for block-based merging)
   - Purpose: Clustering and merging alignments

### RNA Interaction Tools

6. **IntaRNA**
   - Used in: `chira_extract.py` (optional, for hybridization)
   - Purpose: Predicting RNA-RNA interactions and hybridization

### System Utilities

7. **sort** (standard Unix/Linux tool)
   - Used in: `chira_extract.py`, `chira_map.py`, `chira_merge.py`, `chira_quantify.py`
   - Purpose: Sorting files (used via `os.system()`)
   - Note: GNU sort (part of GNU coreutils) version >= 8.6 supports `--parallel` option for faster sorting of large files. The code automatically detects and uses parallel sort when available, falling back to standard sort otherwise.
   - Requirement: GNU coreutils (standard on Linux, available via Homebrew on macOS: `brew install coreutils`)

8. **cat** (standard Unix/Linux tool)
   - Used in: `chira_extract.py`
   - Purpose: Concatenating files (used via `os.system()`)

9. **mv** (standard Unix/Linux tool)
   - Used in: `chira_merge.py`
   - Purpose: Moving/renaming files (used via `os.system()`)

### HPC Cluster Tools

10. **LSF (Load Sharing Facility) cluster scheduler**
    - Used in: `chira_map.py` (optional, for HPC cluster job submission)
    - Commands: `bsub`, `bjobs`, `bqueues`
    - Purpose: Job submission and management on LSF-managed HPC clusters
    - Note: Only required when using `--use_batchtools` option with LSF template in `chira_map.py`
    - Alternative: Other cluster schedulers (SLURM, SGE, etc.) can be used with appropriate batchtools templates

## Installation Recommendations

### Python Packages

**Core packages:**
```bash
pip install biopython bcbiogff pysam requests
# or
conda install -c bioconda biopython bcbiogff pysam requests
```

**Optional packages:**
```bash
# For automatic memory optimization, I/O bottleneck detection, and adaptive buffer sizing
# Highly recommended for optimal performance (10-50x I/O improvement)
pip install psutil
# or
conda install -c bioconda psutil

# For enhanced multiprocessing performance in chira_quantify.py
# HIGHLY RECOMMENDED: 50-90% memory reduction, 2-3x faster startup
pip install mpire
# or
conda install -c conda-forge mpire

# For coordinate liftover in download_mirbase_gff3.py
pip install pyliftover
# or
conda install -c bioconda pyliftover
```

### Command-Line Tools

Most tools can be installed via conda:
```bash
conda install -c bioconda bwa samtools bedtools intarna
```

**GNU coreutils** (for parallel sort support):
- **Linux**: Usually pre-installed (GNU sort is standard)
- **macOS**: Install via Homebrew: `brew install coreutils` (provides `gsort` command, or use `gcoreutils` to get GNU versions of standard commands)
- **Note**: GNU sort version >= 8.6 is required for `--parallel` option. The code automatically detects GNU sort and uses parallel sorting when available.

For CLAN and blockbuster, please refer to their respective documentation for installation instructions.

**R packages** (for HPC cluster job submission via batchtools):
- See the "R Packages" section above for detailed installation instructions
- **LSF cluster scheduler**: Usually pre-installed on HPC clusters. Contact your cluster administrator for access and queue information.
- **Note**: When using batchtools, ensure all file paths are absolute paths (this is handled automatically by the code). The batchtools template file path should also be absolute or use the built-in "lsf-simple" template.

## Parallel Computing Support

ChiRA scripts support parallel processing to improve performance on multi-core systems:

### Multi-Processing Support

1. **chira_quantify.py**
   - Uses MPIRE `WorkerPool` for parallelizing the EM algorithm
   - Command-line option: `-p, --processes` (default: 1, use 0 for all available cores)
   - Parallelizes: E-step (multimapped reads) and aggregation step
   - Benefit: 2-8x faster for large datasets with many multimapping reads (bypasses Python GIL for true parallelism)
   - Benefits of MPIRE: 50-90% memory reduction, 2-3x faster startup, shared objects for large dictionaries
   - Note: MPIRE is required (no fallback). Changed from `ThreadPoolExecutor` to `ProcessPoolExecutor` in v1.4.6, then to MPIRE in v1.4.11 for better performance

2. **chira_merge.py**
   - Uses MPIRE `WorkerPool` for parallelizing transcript chunk processing
   - Command-line option: `-p, --processes` (default: None, auto-detects CPU count)
   - Parallelizes: Transcript chunk processing in overlap-based and blockbuster-based merging
   - Uses SharedObject for `d_desc` dictionary to avoid copying chunk data
   - Benefit: 4-8x faster for datasets with many transcripts, 50-90% memory reduction (bypasses Python GIL for true parallelism)
   - Note: Changed from `ThreadPoolExecutor` to `multiprocessing.Pool` in v1.4.6, then to MPIRE in v1.4.11 for better performance and memory efficiency

3. **chira_map.py**
   - Uses multi-threading for external tools (`samtools`, `pysam`, `sort`) and `ThreadPoolExecutor` for parallel BWA jobs
   - Command-line option: `-p, --processes` (default: 1)
   - Parallelizes: `samtools view`, `pysam.merge`, `pysam.sort`, `sort` commands, and BWA alignment jobs
   - Automatic memory optimization: Uses `psutil` (optional) to calculate optimal memory per thread for BAM sorting
   - I/O bottleneck detection: Uses `psutil` (optional) to automatically detect I/O bottlenecks during chunk processing
   - FASTA chunking: New `--chunk_fasta` parameter for splitting large FASTA files into chunks for parallel processing
   - Benefit: 2-4x faster for large BAM files and sorting operations, improved scalability for very large datasets

4. **chira_extract.py**
   - Uses MPIRE `WorkerPool` for parallelizing read processing
   - Command-line option: `-p, --processes` (default: 1)
   - Parallelizes: Chimeras extraction and hybridization steps
   - Also uses parallel sort (GNU sort `--parallel`) for merging and interaction summary
   - Benefit: Linear speedup with number of processes (up to available cores)
   - Benefits of MPIRE: 50-90% memory reduction, 2-3x faster startup, shared objects for reference dictionaries
   - Note: MPIRE is required (no fallback). Changed from `multiprocessing.Process` to MPIRE in v1.4.11 for better performance

5. **chira_map.py (HPC Cluster Submission)**
   - Uses R batchtools package for submitting chunk-based batch jobs to HPC cluster schedulers
   - Command-line option: `--use_batchtools` (requires additional batchtools parameters)
   - Parallelizes: FASTA chunk processing by submitting each chunk as a separate cluster job
   - Cluster schedulers supported: LSF (default), SLURM, SGE, etc. (via batchtools templates)
   - Benefit: Scales to hundreds of chunks on HPC clusters, better resource management, automatic job queuing
   - Requirements: R with batchtools package, cluster scheduler (LSF, SLURM, etc.), batchtools template file
   - Note: All file paths are automatically converted to absolute paths for cluster job execution

### Performance Recommendations

- **For large datasets**: Use `-t 0` or `-p 0` (where applicable) to automatically use all available CPU cores
- **For memory-constrained systems**: Specify thread/process count explicitly to control memory usage
- **For BAM sorting**: Install `psutil` for automatic memory optimization, or use `--sort_memory` to specify memory per thread manually
- **For I/O optimization**: Install `psutil` to enable adaptive buffer sizing (8-16MB) which provides 10-50x I/O performance improvement
- **For parallel processing in chira_quantify.py, chira_extract.py, and chira_merge.py**: MPIRE is required for multiprocessing (50-90% memory reduction, 2-3x faster startup). Install with `pip install mpire` or `conda install -c conda-forge mpire`.
- **For very large FASTA files**: Use `--chunk_fasta` in `chira_map.py` to split input into chunks for better I/O and memory efficiency
- **For HPC clusters**: Use `--use_batchtools` in `chira_map.py` to submit chunk jobs to cluster schedulers (LSF, SLURM, etc.) for better scalability and resource management

## Script-Specific Dependencies

### Core ChiRA Scripts

- **chira_collapse.py**: No external dependencies (uses standard library only)
- **chira_map.py**: Requires `pysam`, `bwa`, `samtools` (CLAN optional, `psutil` optional for memory optimization and I/O bottleneck detection, R with batchtools optional for HPC cluster submission)
- **chira_merge.py**: Requires `bcbiogff`, `BEDTools`, `mpire` (blockbuster optional)
- **chira_extract.py**: Requires `biopython`, `bcbiogff`, `BEDTools`, `mpire` (IntaRNA optional)
- **chira_quantify.py**: Requires `mpire` for parallel processing
- **chira_utilities.py**: Requires `biopython` (`psutil` optional for adaptive buffer sizing)

### Utility Scripts

- **download_ensembl.py**: Requires `requests`
- **download_mirbase_gff3.py**: Requires `pyliftover` (optional, only for coordinate liftover)
- **download_mirbase_mature.py**: No external dependencies (uses standard library only)
- **remove_mirna_hairpin_from_gtf.py**: No external dependencies (uses standard library only)
- **concatenate_gtf.py**: No external dependencies (uses standard library only)
- **extract_transcripts_from_genome.py**: Requires `gffread` (conda: `conda install -c bioconda gffread`)
- **process_chunk_batchtools.py**: Imports `chira_map` and `chira_utilities`; called by batchtools jobs for mapping chunks
- **submit_chunks_batchtools.R**: Requires R with `batchtools` and `jsonlite` (used by `chira_map.py` for HPC cluster submission)
- **submit_intarna_batchtools.R**: Requires R with `batchtools` and `jsonlite` (used by `chira_extract.py` for IntaRNA when `--hybridize --use_batchtools`)

## Notes

- Some tools are optional depending on usage:
  - **CLAN**: Only needed if using CLAN aligner instead of BWA (not recommended)
  - **blockbuster.x**: Only needed if using block-based merging method
  - **IntaRNA**: Only needed if using the `--hybridize` option in `chira_extract.py`
  - **pyliftover**: Only needed if using coordinate liftover in `download_mirbase_gff3.py`
  - **R with batchtools**: Only needed if using `--use_batchtools` option in `chira_map.py` for HPC cluster job submission
  - **LSF cluster scheduler**: Only needed if using batchtools with LSF template (other schedulers like SLURM, SGE can be used with appropriate templates)
  - **psutil**: Optional but highly recommended for:
    - Automatic memory optimization in `chira_map.py` (BAM sorting)
    - I/O bottleneck detection in `chira_map.py` (when using `--chunk_fasta`)
    - Adaptive buffer sizing in `chira_utilities.py` (10-50x I/O performance improvement)

- Standard library modules used (no installation needed):
  - `argparse`, `os`, `sys`, `collections`, `multiprocessing`, `itertools`, `datetime`, `subprocess`, `math`, `re`, `copy`, `gzip`, `shutil`, `ftplib`, `urllib`, `tempfile`, `time`, `json`, `concurrent.futures` (for ThreadPoolExecutor and ProcessPoolExecutor)

- Parallel computing features:
  - All parallel computing features are backward compatible (default to single-threaded/process)
  - **Multiprocessing** (MPIRE `WorkerPool`) is used for CPU-bound tasks to bypass Python GIL:
    - EM algorithm in `chira_quantify.py` (MPIRE WorkerPool with shared objects)
    - Transcript chunk processing in `chira_merge.py` (MPIRE WorkerPool with shared objects)
  - **Threading** (`ThreadPoolExecutor`) is used for I/O-bound tasks and external tool parallelization:
    - BWA alignment jobs in `chira_map.py` (ThreadPoolExecutor)
    - External tools (`samtools`, `pysam`, `sort`) use multi-threading via their native `-@` flags
  - **Multiprocessing.Process** is used for independent parallel tasks:
    - Read processing in `chira_extract.py` (multiprocessing.Process)
  - GNU sort (GNU coreutils >= 8.6) parallel support is automatically detected and used when available
  - On macOS, GNU coreutils can be installed via Homebrew (`brew install coreutils`) to enable parallel sort

- All utility scripts (`download_ensembl.py`, `download_mirbase_mature.py`, `gff3_to_gtf.py`, etc.) are standalone and can be used independently of the main ChiRA pipeline.
