#!/usr/bin/env python
"""Setup script for ChiRA package."""

from setuptools import setup, find_packages
import os

# Read version from chira_utilities.py (single source of truth)
# The version is defined as __version__ = "X.Y.Z" in chira_utilities.py
version = "1.4.14"  # Default fallback version
try:
    # Import chira_utilities to get the version
    import chira_utilities
    version = chira_utilities.__version__
except (ImportError, AttributeError):
    # If import fails, try to extract from file directly
    try:
        with open("chira_utilities.py", "r") as f:
            for line in f:
                if '__version__' in line and '=' in line:
                    import re
                    # Match pattern: __version__ = "X.Y.Z" or __version__ = 'X.Y.Z'
                    match = re.search(r"__version__\s*=\s*['\"]([\d]+\.[\d]+\.[\d]+)['\"]", line)
                    if match:
                        version = match.group(1)
                        break
    except Exception:
        # If extraction fails, use default version
        pass

# Read long description from README
long_description = ""
if os.path.exists("README.md"):
    with open("README.md", "r", encoding="utf-8") as f:
        long_description = f.read()

setup(
    name="chira",
    version=version,
    description="Chimeric Read Analyzer - Tools for analyzing RNA-RNA interactome data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="ChiRA Contributors",
    url="https://github.com/haibol2016/chira",
    py_modules=[
        "chira_map",
        "chira_merge",
        "chira_extract",
        "chira_collapse",
        "chira_quantify",
        "chira_utilities",
    ],
    scripts=[
        "chira_map.py",
        "chira_merge.py",
        "chira_extract.py",
        "chira_collapse.py",
        "chira_quantify.py",
        "concatenate_gtf.py",
        "download_ensembl.py",
        "download_mirbase_gff3.py",
        "download_mirbase_mature.py",
        "extract_transcripts_from_genome.py",
        "remove_mirna_hairpin_from_gtf.py",
        "process_chunk_batchtools.py",
        "process_intarna_chunk_batchtools.py",
        "merge_intarna_into_chimeras.py",
        "submit_chunks_batchtools.R",
        "submit_intarna_batchtools.R",
    ],
    install_requires=[
        # Core Python dependencies (required for basic functionality)
        "biopython",      # For FASTA/sequence parsing (chira_extract.py, chira_utilities.py)
        "bcbio-gff",      # For GFF/GTF annotation file parsing (chira_extract.py, chira_merge.py)
        "pysam",          # For BAM file manipulation (chira_map.py: merging, sorting BAM files)
        "mpire>=2.4.0",  # For enhanced multiprocessing performance (chira_quantify.py, chira_extract.py, chira_merge.py)
                          # - Provides shared objects, lower overhead, and better performance than ProcessPoolExecutor/Process/Pool
                          # - Benefits: 50-90% memory reduction, 2-3x faster startup, better performance
                          # - Required for parallel processing in chira_quantify.py, chira_extract.py, and chira_merge.py (no fallback)
                          # - Version >= 2.4.0 required for shared_objects module support
                          # - Changed from optional to required in v1.4.11
                          # - Uses start_method='fork' on Unix/Linux for copy-on-write semantics (optimal memory efficiency)
    ],
    extras_require={
        "optional": [
            # Optional Python packages (highly recommended for optimal performance)
            "psutil",          # For automatic memory optimization, I/O bottleneck detection, and adaptive buffer sizing
                              # - chira_map.py: Automatic memory allocation for BAM sorting, I/O bottleneck detection
                              # - chira_utilities.py: Adaptive buffer sizing (8-16MB) for 10-50x I/O performance improvement
                              # Install with: pip install psutil
            "requests",        # For downloading files via HTTP/HTTPS
                              # - Used in: download_ensembl.py, download_mirbase_gff3.py, download_mirbase_mature.py
            "pyliftover",      # For download_mirbase_gff3.py coordinate liftover (converting coordinates between genome versions)
        ],
        # Install all optional dependencies: pip install chira[optional]
    },
    # External Command-Line Tools (must be installed separately, not via pip):
    # These tools are required for ChiRA to function and must be available in system PATH.
    # Most can be installed via conda: conda install -c bioconda <tool_name>
    #
    # REQUIRED TOOLS:
    # - BWA (Burrows-Wheeler Aligner)
    #   * Used in: chira_map.py
    #   * Commands: bwa mem, bwa index
    #   * Install: conda install -c bioconda bwa
    #
    # - samtools
    #   * Used in: chira_map.py (samtools view for SAM to BAM conversion)
    #   * Commands: samtools view
    #   * Install: conda install -c bioconda samtools
    #
    # - BEDTools
    #   * Used in: chira_extract.py, chira_merge.py
    #   * Commands: bedtools intersect (or intersectBed), bedtools getfasta (or fastaFromBed)
    #   * Code automatically detects and supports both old and new BEDTools command formats
    #   * Install: conda install -c bioconda bedtools
    #
    # - sort (standard Unix/Linux tool)
    #   * Used in: chira_extract.py, chira_map.py, chira_merge.py, chira_quantify.py
    #   * GNU sort >= 8.6 recommended for parallel sort support (--parallel option)
    #   * Linux: Usually pre-installed (GNU sort is standard)
    #   * macOS: Install via Homebrew: brew install coreutils
    #   * Code automatically detects and uses parallel sort when available
    #
    # OPTIONAL TOOLS:
    # - CLAN (Chimeric Long-read Aligner)
    #   * Used in: chira_map.py (alternative to BWA, not recommended due to performance)
    #   * Commands: clan_search, clan_output, clan_index
    #   * Only needed if using CLAN aligner instead of BWA
    #
    # - blockbuster.x
    #   * Used in: chira_merge.py (optional, for block-based merging method)
    #   * Install: Refer to blockbuster documentation
    #
    # - IntaRNA
    #   * Used in: chira_extract.py (optional, for RNA-RNA interaction prediction/hybridization)
    #   * Only needed if using --hybridize option in chira_extract.py
    #   * Install: conda install -c bioconda intarna
    #
    # - gffread
    #   * Used in: extract_transcripts_from_genome.py (optional utility script)
    #   * Part of GFF utilities from Johns Hopkins University
    #   * Install: conda install -c bioconda gffread
    #
    # CLUSTER COMPUTING (for parallel processing on HPC clusters):
    # - R with batchtools and jsonlite packages (R packages, not Python packages)
    #   * Used in: chira_map.py --use_batchtools, chira_extract.py --use_batchtools (LSF/SLURM/SGE)
    #   * Enables distributing chunk processing across cluster nodes via batchtools
    #   * batchtools: Submitting chunk-based batch jobs to HPC cluster schedulers
    #   * jsonlite: Parsing JSON configuration files for batchtools job submission
    #   * Install R: conda install -c conda-forge r-base
    #   * Install R packages: conda install -c conda-forge r-batchtools r-jsonlite
    #   * Or in R: install.packages(c("batchtools", "jsonlite"))
    #   * jsonlite is usually installed automatically as a dependency of batchtools
    #   * Required only if using --use_batchtools option in chira_map.py
    #   * All file paths are automatically converted to absolute paths for cluster job execution
    #   * See BATCHTOOLS_USAGE.md for detailed usage instructions
    #
    # Quick installation of all required tools:
    #   conda install -c bioconda bwa samtools bedtools
    #
    # For complete dependency information, see DEPENDENCIES.md
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)

