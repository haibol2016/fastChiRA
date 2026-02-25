# Changelog

All notable changes to ChiRA will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.4.14] - 2026-02-25

### Changed
- Version bump to 1.4.14. Documentation updates: final output table column names (chimeras, singletons, interactions) aligned with code; CLI argument docs (chira_extract batchtools: `--remove_intermediate`, removed `--keep_batchtools_work` and `--batchtools_poll_interval`); R packages `future`/`future.apply` no longer required (removed from docs).

## [1.4.13] - 2026-02-22

### Added
- **merge_intarna_into_chimeras.py**: Standalone script to merge IntaRNA results into chimeras and produce final `{sample_name}.chimeras.txt`, `{sample_name}.singletons.txt`, and `{sample_name}.interactions.txt`. Use when `chira_extract.py --hybridize --use_batchtools` times out after submitting IntaRNA jobs: run after all IntaRNA chunk jobs complete, with required `--outdir`, `--sample-name`, `--n-chunks`. Processes all chunks (merge `loci_seqs.pkl` + `result.csv` → `*.chimeras-r.<n>`), then merges chimeras-r and singletons and writes interactions.

### Fixed
- **submit_intarna_batchtools.R**: Wait for all batches (including the last) so the main job does not finish while IntaRNA jobs are still running. Replaced `waitForJobs()` with manual polling using `getJobTable()` and status counting to avoid `'sum' not defined for "POSIXt" objects` when batchtools returns POSIXct for done/running/error.
- **submit_chunks_batchtools.R**: Same POSIXct-safe status counting (`count_status()`) and manual polling instead of `waitForJobs()` when waiting between batches.
- **chira_map.py**: Embedded R in `query_batchtools_status()` now uses `count_status()` for job table columns so status queries work when batchtools returns POSIXct.
- **chira_extract.py**: IntaRNA result CSV is read with `f.read().splitlines()` so files with `\r` or `\r\n` line endings are parsed correctly (avoids false "empty or header-only" warning). Parser accepts both `;` and `,` separators and requires at least 8 fields per row.

### Changed
- **submit_intarna_batchtools.R**, **submit_chunks_batchtools.R**: Added `count_status()` helper: for POSIXct/POSIXt columns count non-NA; otherwise sum logical. All job status counts use this helper. Wait logic uses a repeat loop with `getJobTable()` and `Sys.sleep(30)` instead of `waitForJobs()`.
- **chira_extract.py** (IntaRNA/batchtools): IntaRNA CSV no longer requests or parses seq1/seq2; full sequences are taken from `loci_seqs.pkl` (written in prepare, loaded in finish). Dot-bracket in interactions follows IntaRNA convention (id1 = '(', id2 = ')') without swapping. `_merge_hybrid_into_chimeras` always uses `d_loci_seqs` for seq1/seq2.

## [1.4.12] - 2026-02-20

### Added
- **chira_extract.py**: Added `gene_name_1` and `gene_name_2` columns to the interaction summary table (`{sample_name}.interactions.txt`)
  - Extracts gene symbols from the chimeras file (columns 5 and 6: `gene_symbol_1` and `gene_symbol_2`)
  - Displays all unique gene names involved in each interaction (semicolon-separated if multiple)
  - Columns are added at the end of the interaction summary table, following the same pattern as `reference_transcript_id_1` and `reference_transcript_id_2`
  - Properly handles interaction reversal when `interaction_otherway` is detected (gene names are swapped along with other fields)

## [1.4.11] - 2026-02-17

### Fixed
- **chira_quantify.py**: Fixed CRL ID validation bug (assumed sequential IDs) and EM algorithm return bug (returned initial instead of optimized values)
- **chira_map.py**: Fixed alternate alignment filtering for stranded RNA-seq (quality thresholds now based only on correct-strand alignments)
- **chira_extract.py**: Fixed field index mismatches in `loci.counts` parsing, MPIRE argument passing (dictionaries → tuples), and collapsed 8 separate `hybridization_genomic_coordinates` columns into a single column

### Changed
- **MPIRE is now a required dependency** (moved from optional to `install_requires` in setup.py)
  - All parallel processing now uses MPIRE WorkerPool with shared objects
  - Added `start_method='fork'` for copy-on-write semantics (50-90% memory reduction on Unix/Linux)
  - Removed ProcessPoolExecutor/Process fallback code
  - Benefits: 50-90% memory reduction, 2-3x faster startup, shared memory access

## [1.4.10] - 2026-02-15

### Fixed
- **chira_map.py**: Fixed batchtools path handling (all paths now absolute) and JSON parsing errors
- **submit_chunks_batchtools.R**: Added file existence checks and improved error handling for JSON parsing

### Changed
- **Code refactoring**: All scripts refactored to extract `parse_arguments()` and modularize `main()` functions for better code organization
- **chira_map.py**: All batchtools job paths are now absolute and normalized for cluster execution

## [1.4.9] - 2026-02-15

### Added
- **chira_map.py**: Added `--parallel_chunks` parameter to control chunk parallelism (default: 2)

### Fixed
- **chira_map.py**: Fixed CPU guidance display to use parallel chunks instead of total chunks

## [1.4.8] - 2026-02-15

### Changed
- **chira_merge.py**: Improved chunk-based parallelization for very large datasets (fixed ~1000 transcripts per chunk), updated variable naming (chrom → transcript)

### Fixed
- **chira_merge.py**: Fixed `AttributeError` when processing empty files, corrected variable naming inconsistencies

## [1.4.7] - 2026-02-15

### Added
- **extract_transcripts_from_genome.py**: New utility script to extract transcript FASTA from genome FASTA using gffread
- **Dockerfile**: Added `gffread` package

### Changed
- **download_mirbase_mature.py**: Added automatic U→T conversion in mature miRNA sequences
- **concatenate_gtf.py**: Updated to support miRBase GFF3 format directly

### Removed
- **remove_mirna_hairpin_from_fasta.py**: Replaced by `extract_transcripts_from_genome.py`
- **concatenate_fasta.py**: No longer needed (use standard Unix tools)

## [1.4.6] - 2026-02-15

### Added
- **chira_utilities.py**: Added adaptive buffer sizing (8-16MB) based on available RAM (10-50x I/O improvement)
- **chira_map.py**: Added FASTA chunking (`--chunk_fasta`), I/O bottleneck detection, and progress tracking
- **chira_merge.py**: Re-added parallel processing with `multiprocessing.Pool` (4-8x faster)
- **chira_extract.py**: Refactored into modular functions

### Changed
- **Multiprocessing evolution**: ThreadPoolExecutor → ProcessPoolExecutor (bypasses GIL, 2-8x faster)
- **I/O optimizations**: Fixed 2MB buffers → adaptive 8-16MB buffers (10-50x faster for large files)
- **chira_map.py**: Replaced `os.system()` with `subprocess.run()`, added progress tracking
- **chira_merge.py**: Parameter changed from `-t, --threads` to `-p, --processes`

## [1.4.5] - 2026-02-02

### Added
- **Parallel computing support**: Multi-threading in `chira_quantify.py` and `chira_merge.py` (2-4x faster)
- **I/O optimizations**: Added 2MB buffer sizes for file I/O (20-40% faster)
- **GNU sort parallel support**: Automatic detection and usage (2-4x faster sorting)
- **chira_map.py**: Enhanced multi-threading for external tools (samtools, pysam) with automatic memory optimization

### Changed
- **chira_quantify.py**: Shallow copy optimization (10-50x faster), pre-computed values, optimized dictionary operations
- **chira_map.py**: Automatic memory calculation with `psutil`, parallel sort support

- **chira_collapse.py**:
  - Added buffer sizes for I/O operations
  - Added comprehensive comments explaining optimizations

- **DEPENDENCIES.md**:
  - Added `psutil` as optional dependency for automatic memory optimization
  - Added comprehensive "Parallel Computing Support" section
  - Documented GNU coreutils requirement for parallel sort
  - Updated installation instructions for GNU coreutils on different platforms
  - Added performance recommendations and usage notes

### Fixed
- **chira_map.py**:
  - Fixed division by zero protection in memory calculation (use `max(1, args.processes)`)
  - Fixed parallel sort to check `args.processes > 0` before using `--parallel`
  - Fixed redundant calculation: use cached `ref_start_int` instead of recalculating
  - Removed unused variable initialization (`readseq = None`)

- **chira_merge.py**:
  - Fixed missing `d_desc.clear()` in parallel path for consistency
  - Ensured parallel and sequential paths produce identical results

- **chira_quantify.py**:
  - Verified thread safety: each thread processes different read IDs (no race conditions)
  - Verified algorithm correctness: parallel and sequential paths produce identical results

### Performance Improvements
- **EM Algorithm (chira_quantify.py)**: 2-4x faster with multi-threading for large datasets
- **Chromosome Processing (chira_merge.py)**: 2-4x faster with multi-threading for many chromosomes
- **BAM Operations (chira_map.py)**: 2-4x faster with multi-threaded samtools/pysam operations
- **File Sorting**: 2-4x faster with GNU sort parallel support
- **I/O Operations**: 20-40% faster with optimized buffer sizes
- **Overall Pipeline**: Significant speedup for large datasets on multi-core systems

### Compatibility
- **Backward Compatible**: All parallel features default to single-threaded (num_threads=1, processes=1)
- **GNU sort**: Automatically detects and uses parallel sort when available (version >= 8.6)
- **psutil**: Optional dependency with graceful fallback to safe defaults
- **Algorithm Correctness**: All optimizations preserve original algorithm logic

## [1.4.4] - 2026-01-26

### Added
- **chira_extract.py**:
  - Added `--sample_name` parameter (required) for customizable output file names (line 707)
  - Output files now use format: `{sample_name}.chimeras.txt`, `{sample_name}.singletons.txt`, `{sample_name}.interactions.txt`
  - Added `--gzip` option to compress output files (chimeras and singletons) with gzip (lines 835-836, 844)
    - Only final merged files are compressed (intermediate per-process files remain uncompressed for optimal performance)
    - Compressed files have `.gz` extension (e.g., `{sample_name}.chimeras.txt.gz`)
    - Interactions file is always uncompressed for compatibility with downstream tools
    - Compression improves I/O performance for large files and saves disk space
  - Added header row to interactions output file with all column descriptions (lines 614-622)
  - Added comment lines in interactions file explaining how to identify miRNA vs target loci (lines 625-626)
  - Added `mirna_position` column to `{sample_name}.chimeras.txt` output (lines 268-280, 292, 902)
    - Indicates whether miRNA is at 5' end (`miRNA_first`) or 3' end (`miRNA_last`) of the chimeric read
    - Helps identify read orientation: 5' miRNA-target 3' vs 5' target-miRNA 3'

- **chira_utilities.py**:
  - Added `get_bedtools_command()` function for automatic BEDTools version detection (lines 124-160)
  - Supports both old (individual commands) and new (unified `bedtools` command) formats
  - Module-level regex pattern compilation for CIGAR parsing (lines 8-11)
  - Caching dictionary for BEDTools commands (line 14)

- **DEPENDENCIES.md**:
  - Created comprehensive dependency documentation
  - Documented all Python packages and command-line tools

- **CHANGELOG.md**:
  - Added changelog file to track all modifications

- **download_ensembl.py** (new utility script):
  - Created script to download Ensembl files (cDNA, ncRNA, GTF, genome FASTA)
  - Supports downloading from specific Ensembl release versions
  - Auto-detects assembly names or accepts `--assembly` parameter
  - Downloads primary assembly (not toplevel) for genome FASTA
  - Supports HTTP and FTP protocols with automatic fallback
  - Automatically decompresses gzipped files
  - Parameters: `--species`, `--genome-version`, `--gtf-version`, `--output-dir`

- **download_mirbase_mature.py** (new utility script):
  - Created script to download species-specific mature miRNA sequences from miRBase
  - Supports specific miRBase versions or CURRENT version
  - Extracts sequences by species code (e.g., hsa, mmu, bta)
  - Handles gzipped files and automatic decompression
  - Parameters: `--species`, `--output`, `--mirbase-version`

- **download_mirbase_gff3.py** (new utility script):
  - Created script to download species-specific GFF3 annotation files from miRBase
  - Supports CURRENT version (default) or specific version via `--mirbase-version` parameter
  - Supports chromosome name mapping via `--chromosome_mapping` parameter
  - **Coordinate liftover support**: Convert coordinates between genome assemblies (e.g., hg19 → hg38) using pyliftover
    - Parameters: `--source-genome`, `--target-genome`, `--chain-file` for liftover functionality
    - Processing order: Download → Liftover → Rename chromosomes
    - Handles GFF3 1-based coordinates correctly (converts to 0-based for pyliftover, back to 1-based for output)
    - Updates chromosome names if they change during liftover
    - Features that cannot be lifted over retain their original coordinates
  - Can be used directly with ChiRA (no GTF conversion needed)
  - Parameters: `--species`, `--output`, `--mirbase-version`, `--chromosome_mapping`, `--source-genome`, `--target-genome`, `--chain-file`

- **concatenate_fasta.py** (new utility script):
  - Created script to concatenate multiple FASTA files into a single file
  - Handles various FASTA header formats
  - Useful for combining miRNA and target transcriptome FASTA files
  - Parameters: `--input-files`, `--output`

- **remove_mirna_hairpin_from_gtf.py** (new utility script):
  - Created script to remove microRNA entries from Ensembl GTF files
  - Identifies miRNA entries by feature type, biotype, and optional regex pattern
  - Supports custom regex pattern via `--pattern` parameter for flexible matching
  - Preserves comment lines (optional removal)
  - Used for preparing target-only transcriptome annotations

- **remove_mirna_hairpin_from_fasta.py** (new utility script):
  - Created script to remove miRNA FASTA records from transcriptome FASTA files
  - Uses transcript IDs extracted from GTF file to identify miRNA sequences
  - Supports various FASTA header formats
  - Reuses miRNA detection logic from `remove_mirna_hairpin_from_gtf.py`
  - Used for preparing target-only transcriptome FASTA files

- **concatenate_gtf.py** (new utility script):
  - Created script to concatenate mature miRNA GTF with target transcriptome GTF
  - Removes comment lines from miRNA GTF file
  - Optionally removes comment lines from target GTF
  - Produces combined GTF file for split-reference analysis
  - Parameters: `--mirna-gtf`, `--target-gtf`, `--output`

- **Dockerfile**:
  - Created Docker image definition using micromamba base image
  - Includes all core dependencies (biopython, bcbiogff, pysam, requests, pyliftover)
  - Includes bioinformatics tools (bwa, samtools, bedtools, intarna)
  - Sets up proper environment variables and PATH
  - Makes Python scripts executable
  - Supports containerized execution of ChiRA pipeline

### Changed
- **chira_collapse.py**:
  - Removed Biopython dependency - replaced `SeqIO.parse()` with raw file parsing
  - Implemented direct FASTQ line-by-line parsing (2-5x faster for large files)
  - Optimized string formatting using f-strings
  - Cached UMI length check to avoid repeated conditionals
  - **Code changes**: Removed `from Bio import SeqIO` import, replaced `SeqIO.parse()` loop with raw file reading using `line_num % 4 == 2` to extract sequence lines

- **chira_quantify.py**:
  - Changed `d_crl_reads` from `defaultdict(list)` to `defaultdict(set)` for faster set operations (line 33)
  - Replaced `list.extend()` with `set.update()` for better performance
  - Replaced `copy.deepcopy()` with `dict()` for shallow copy in EM algorithm (10-50x faster, line 194)
  - Pre-computed inverse values to avoid repeated divisions
  - Cached sorted CRL IDs and sorted lengths to avoid repeated sorting
  - Reduced print frequency in EM algorithm (every 10 iterations instead of every iteration)
  - Optimized dictionary operations and loop structures (lines 55-72, 93-101, 155-184)
  - Improved file I/O with context managers (lines 256-270)

- **chira_utilities.py**:
  - Pre-compiled regex patterns as module-level constants for CIGAR string parsing (lines 8-11):
    - `_CIGAR_PATTERN_QUERY` for `query_length()` and `match_positions()`
    - `_CIGAR_PATTERN_ALIGN` for `alignment_length()`
    - `_CIGAR_PATTERN_END` for `alignment_end()`
  - Added caching for BEDTools command detection to avoid repeated subprocess calls
  - Improved `extract_reflengths()` to use context manager for proper file handling
  - Optimized `print_w_time()` to use f-string formatting

- **chira_map.py**:
  - Optimized `write_mapped_bed()` function:
    - Only get sequence when needed (for unmapped reads only)
    - Pre-computed desired strand check

    - Two-pass processing: lightweight first pass to find optimal length, then parse and write
    - Merged condition checks for cleaner code
    - Removed unnecessary caching of single-use variables
  - Optimized `clan_to_bed()` function:
    - Simplified field parsing using tuple unpacking
    - Removed unnecessary empty string checks
    - Streamlined CIGAR string building

- **chira_merge.py**:
  - Updated `transcript_to_genomic_pos()` function to use `chira_utilities.get_bedtools_command('intersect')`
  - Replaced hardcoded `intersectBed` command with version-agnostic function call
  - Added zero-length match checks in `write_segments()` to prevent division by zero
  - Removed unnecessary caching optimizations that didn't provide performance benefit
  - Optimized set operations in `transcript_to_genomic_pos()` for deduplication
  - Enhanced UTR parsing to support Ensembl-specific UTR types (`five_prime_utr`, `three_prime_utr`)
    - Updated `limit_info` to include all UTR types for better Ensembl GTF compatibility
  - Added Biopython deprecation warning suppression (lines 10-21)
    - Suppresses `BiopythonDeprecationWarning` related to `UnknownSeq(length)` from `bcbiogff` package
    - Includes fallback for older Biopython versions that don't have `BiopythonDeprecationWarning`

- **chira_extract.py**:
  - Updated to use `chira_utilities.get_bedtools_command('getfasta')` for BEDTools compatibility (line 803)
  - Updated all function signatures to include `sample_name` parameter:
    - `write_chimeras()` (line 319)
    - `hybridize_and_write()` (line 364)
    - `write_interaction_summary()` (line 540)
  - Updated all file paths to use sample name prefix (lines 322-323, 372, 542, 776-778, 894)
  - Added `--sample_name` argument to argument parser (line 707)
  - Improved output file headers with more descriptive column names (lines 912-945, 947-961, 678-689)
    - Chimeras file: Updated all 34 column names (e.g., `tagid` → `read_id`, `txid1` → `transcript_id_1`, `mfe` → `hybridization_mfe_kcal_mol`)
    - Singletons file: Updated all 14 column names (e.g., `tagid` → `read_id`, `txid` → `transcript_id`)
    - Interactions file: Updated all 24 column names (e.g., `read_count` → `supporting_read_count`, `locus1_chr` → `locus_1_chromosome`)
  - Refactored `merge_files()` function to use Python for header writing and shell commands for sorting (lines 597-644)
    - Header is written directly in Python to avoid shell escaping issues
    - Uses `subprocess.Popen()` for better error handling
    - Improved file path escaping for shell safety
    - Intermediate files are always uncompressed (only final merged files are compressed if `--gzip` is used)

- **README.md**:
  - Added version information and brief summary of improvements
  - Added note about modified version and GPL compliance
  - Streamlined to focus on user-facing documentation, with detailed changes in CHANGELOG.md
  - Updated to reflect new header names in output files (chimeras, singletons, interactions)
  - Added documentation for `--gzip` option in `chira_extract.py`
  - Added Singularity/Apptainer support documentation with reference to `SINGULARITY_SETUP.md`
  - Updated installation section to include both Docker and Singularity examples
  - Fixed typo: `--summerize` → `--summarize` in parameter documentation (lines 577, 651)
  - Updated all column name references to match new descriptive headers

- **DEPENDENCIES.md**:
  - Updated to reflect that `chira_collapse.py` no longer uses Biopython
  - Updated BEDTools commands to reflect automatic version detection
  - Added documentation for new utility scripts and their dependencies
  - Added optional dependencies: pyliftover (for `download_mirbase_gff3.py` coordinate liftover), requests (for `download_ensembl.py`)
  - Documented script-specific dependencies

### Fixed
- **chira_map.py**:
  - Fixed bug where first iteration would write empty string to unmapped FASTA file (line 67-70)
    - Added check: `if prev_readid is not None` before writing previous read
  - Fixed typos: "STRAT" → "START" in all print statements (lines 416, 422, 426, 431)
  - Fixed typo: "Funtion" → "Function" in docstring (line 12)
  - Fixed typo: "alignmnet" → "alignment" in docstring (line 15)
  - Fixed typos: "prioroty" → "priority" in argument help text (lines 249, 255)
  - Fixed inefficiency in `clan_to_bed()` function (line 222)
    - Changed to use pre-split `location_list` instead of re-splitting `mapped_locations`

- **chira_merge.py**:
  - Added zero-length match checks in `write_segments()` to prevent division by zero (lines 134-138, 153-156)
    - Checks `first_match_length` and `current_match_length` before division
    - Checks `last_match_length` and `current_match_length` before division
  - Removed unnecessary caching that didn't improve performance
  - Fixed typos: "alignmnets" → "alignments" (line 33), "alignmnet" → "alignment" (lines 95, 609)
  - Fixed typo: "postion" → "position" (line 33)
  - Fixed typo: "prioroty" → "priority" (line 632)
  - Fixed typos: "blockbuser" → "blockbuster" (line 655), "blockcuster" → "blockbuster" (lines 698, 701)
  - Enhanced UTR parsing to support Ensembl-specific UTR types (lines 436, 459-484)
    - Added support for `five_prime_utr` and `three_prime_utr` feature types in addition to generic `UTR`
    - Updated `limit_info` to include all UTR types for better Ensembl GTF compatibility

- **chira_extract.py**:
  - Fixed undefined variables `d_reflen1`/`d_reflen2` in `extract_and_write()` (lines 249-257)
    - Changed to use correct parameter names `d_ref_lengths1`/`d_ref_lengths2`
  - Fixed logic error in strand/chromosome check in `guess_region()` (line 85)
    - Changed `and` to `or` for correct conditional logic
  - Fixed index out of bounds error in `parse_counts_file()` (lines 539-553)
    - Clamped `tpm_cutoff` to `[0, 1)` range
    - Added check for empty `uniq_tpms` list
  - Fixed logic error in `hybridization_positions()` function (lines 573-582)
    - Corrected iteration logic to find last `'('` in `dotbracket1` and last `')'` in `dotbracket2`
    - Refactored loops for clarity using `range(len(...))`
  - Fixed potential string slicing issue in `hybridize_and_write()` (lines 394-398)
    - Added length check before slicing `record.id[:-3]`
  - Fixed logic bug in `write_chimeras()` for processing last read (line 398)
    - Changed condition from `if read_count >= total_read_count:` to `if l_readlines and chunk_start <= read_count + 1 <= chunk_end:`
    - Ensures last read is processed correctly when it's within the chunk range
  - Fixed typos: "prioroty" → "priority" (line 763), "summerize" → "summarize" (lines 771, 772, 796, 965)
  - Fixed typo: "Outpur direcoty" → "Output directory" (line 782)
  - Enhanced UTR parsing to support Ensembl-specific UTR types (lines 445, 475-484, 53-65, 85-95)
    - Added support for `five_prime_utr` and `three_prime_utr` feature types in addition to generic `UTR`
    - Updated `limit_info` to include all UTR types for better Ensembl GTF compatibility
    - Modified `guess_region()` to use stored `utr_type` for more specific region assignment (5_prime_UTR vs 3_prime_UTR)
    - Position-based UTR determination now only applies to generic 'UTR' types
  - Improved code maintainability with constants for magic indices (lines 19-26, 29)
    - Defined constants: `CHIMERA_IDX_LOCUS1`, `CHIMERA_IDX_LOCUS2`, `CHIMERA_IDX_SEQUENCES`, etc.
    - Defined constant: `MIRNA_REGION_TYPES` for miRNA region type checking
  - Refactored `mirna_position` handling (lines 295-309, 430-433)
    - Pre-populated `chimera` list with "NA" values for hybridization fields
    - Eliminated need for `pop()` and `append()` operations later in code
    - Removed unnecessary fallback code after refactoring

- **chira_utilities.py**:
  - Fixed `median()` function bug where input wasn't sorted (lines 28-35)
    - Added `x_sorted = sorted(x)` before calculating median
    - Changed `int(n/2)` to `n // 2` for integer division
  - Fixed `get_bedtools_command()` to verify command success (line 157)
    - Added `process.returncode == 0` check before caching result
  - Fixed file handling in `extract_reflengths()` to use proper context managers

- **chira_quantify.py**:
  - Fixed CRL iteration range to properly include index 0 in reverse traversal (line 55)
    - Changed from `range(len(d_crl_reads) - 1, 0, -1)` to `range(len(d_crl_reads) - 1, -1, -1)`
  - Fixed typos: "2st" → "2nd" in print statements (lines 111, 128)

- **BEDTools compatibility**:
  - Fixed command compatibility issues across different BEDTools versions
  - Now automatically detects and uses appropriate command format
  - Implemented in `chira_merge.py` and `chira_extract.py`

### Performance Improvements
- **Overall**: 3-10x faster processing for `chira_quantify.py`
- **CIGAR parsing**: 2-5x faster with pre-compiled regex patterns
- **FASTQ parsing**: 2-5x faster with raw file parsing in `chira_collapse.py`
- **BAM processing**: 20-40% faster in `chira_map.py`
- **EM algorithm**: 10-50x faster with shallow copy instead of deepcopy
- **CRL building**: 2-5x faster with set operations instead of list operations

### Compatibility
- **BEDTools**: Now supports both old format (`intersectBed`, `fastaFromBed`) and new format (`bedtools intersect`, `bedtools getfasta`)
- **Automatic detection**: BEDTools version is automatically detected and appropriate commands are used

### New Utility Scripts
- **download_ensembl.py**: Download Ensembl reference files (cDNA, ncRNA, GTF, genome)
  - Downloads primary assembly genome FASTA (not toplevel)
  - Auto-detects assembly names
  - Supports HTTP and FTP with automatic fallback
- **download_mirbase_mature.py**: Download species-specific mature miRNA sequences from miRBase
  - Extracts sequences by species code
  - Supports specific versions or CURRENT
- **download_mirbase_gff3.py**: Download species-specific GFF3 annotation files from miRBase
  - Supports CURRENT version or specific version
  - **Coordinate liftover**: Convert coordinates between genome assemblies (e.g., hg19 → hg38, mm9 → mm10)
    - Uses pyliftover with UCSC chain files
    - Processing order: Download → Liftover → Rename chromosomes
    - Handles coordinate conversion correctly (GFF3 1-based to pyliftover 0-based and back)
    - Updates chromosome names if changed during liftover
  - **Chromosome name mapping**: Rename chromosomes based on a mapping file
  - Can be used directly with ChiRA (no GTF conversion needed)
- **concatenate_fasta.py**: Concatenate multiple FASTA files into a single file
  - Handles various FASTA header formats
  - Useful for combining miRNA and target transcriptome FASTA files
- **remove_mirna_hairpin_from_gtf.py**: Remove miRNA entries from GTF annotation files
  - Flexible regex pattern matching
  - Used for preparing target-only transcriptome annotations
- **remove_mirna_hairpin_from_fasta.py**: Remove miRNA sequences from FASTA files
  - Uses transcript IDs from GTF to identify sequences
  - Used for preparing target-only transcriptome FASTA files
- **concatenate_gtf.py**: Concatenate miRNA and target GTF files for split-reference analysis
  - Removes comment lines from miRNA GTF
  - Produces combined GTF for use as ref1 annotation

### Docker Support
- **Dockerfile**: Complete Docker image with all dependencies pre-installed
  - Uses micromamba for lightweight conda package management
  - Includes all Python packages (biopython, bcbiogff, pysam, requests, pyliftover)
  - Includes bioinformatics tools (bwa, samtools, bedtools, intarna)
  - Sets up proper environment variables and PATH
  - Makes Python scripts executable
  - Includes entrypoint script that automatically sets up conda environment, PATH, and PYTHONPATH
  - Symlinks all conda binaries to `/usr/local/bin` for universal PATH access
  - Ready-to-use containerized environment for ChiRA pipeline

### Singularity/Apptainer Support
- **SINGULARITY_SETUP.md**: Comprehensive guide for using ChiRA with Singularity/Apptainer containers
  - Installation instructions for Linux and macOS
  - Docker image to Singularity conversion methods
  - Environment setup using entrypoint script (recommended approach)
  - Environment isolation flags and best practices
  - Troubleshooting guide for common PATH and environment issues
  - Complete reference for all Singularity options and flags

## [1.4.3] - Previous Version
- Original version before performance optimizations and new features

---

[1.4.4]: https://github.com/original-repo/chira/compare/v1.4.3...v1.4.4
[1.4.3]: https://github.com/original-repo/chira/releases/tag/v1.4.3

