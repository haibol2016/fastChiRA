# ChiRA Installation Guide

## Installing ChiRA in a Conda Environment

### Step 1: Create and Activate Conda Environment

```bash
# Create a new conda environment with Python 3.9 (or 3.10, 3.11)
conda create -n chira python=3.9
conda activate chira
```

### Step 2: Install Required External Tools

ChiRA requires several command-line tools that must be installed separately:

```bash
# Install required tools from bioconda
conda install -c bioconda bwa samtools bedtools

# Optional: Install additional tools if needed
conda install -c bioconda intarna gffread  # Only if using specific features

# Optional: Install R and batchtools for HPC cluster processing
conda install -c conda-forge r-base r-batchtools r-jsonlite  # Only if using --use_batchtools
```

### Step 3: Install ChiRA Package

#### Option A: Install from Source (Recommended for Development)

```bash
# Navigate to ChiRA directory
cd /path/to/chira

# Install in editable mode (allows code modifications)
pip install -e .

# Or install normally
pip install .
```

#### Option B: Install with Optional Dependencies

```bash
# Install with all optional dependencies (recommended for best performance)
pip install -e .[optional]

# This installs:
# - psutil (for memory optimization and I/O performance)
# - requests (for downloading Ensembl files)
# - pyliftover (for coordinate liftover)

# Note: mpire is now a required dependency (installed automatically with pip install chira)
# mpire provides enhanced multiprocessing performance in chira_quantify.py, chira_extract.py, and chira_merge.py
```

### Step 4: Verify Installation

```bash
# Check that ChiRA scripts are available
chira_map.py --version
chira_merge.py --version
chira_extract.py --version

# Check that external tools are available
bwa
samtools --version
bedtools --version
```

## Complete Installation Example

```bash
# 1. Create environment
conda create -n chira python=3.9
conda activate chira

# 2. Install external tools
conda install -c bioconda bwa samtools bedtools

# 3. Navigate to ChiRA directory
cd /path/to/chira

# 4. Install ChiRA with optional dependencies
pip install -e .[optional]

# 5. Verify installation
chira_map.py --version
```

## Installing for Batchtools Support

If you plan to use batchtools for LSF cluster processing (e.g. `chira_map.py --use_batchtools` or `chira_extract.py --hybridize --use_batchtools`):

```bash
# Activate your conda environment
conda activate chira

# Install R and batchtools (recommended: via conda)
conda install -c conda-forge r-base r-batchtools r-jsonlite

# Or install R first, then packages in R:
conda install -c conda-forge r-base
# Then in R:
# install.packages(c("batchtools", "jsonlite"))
```

**Note**: `jsonlite` is required for JSON configuration file parsing. It's usually installed automatically as a dependency of `batchtools`, but explicitly installing it ensures compatibility.

**Important**: All file paths are automatically converted to absolute paths for cluster job execution. No manual path conversion is needed.

**Scripts**: `submit_chunks_batchtools.R` is used for mapping chunks; `submit_intarna_batchtools.R` is used for IntaRNA hybridization when both `--hybridize` and `--use_batchtools` are set in `chira_extract.py`. If the main job times out after submitting IntaRNA jobs, run `merge_intarna_into_chimeras.py` (see [BATCHTOOLS_USAGE.md](BATCHTOOLS_USAGE.md) and [README.md](README.md)) after all IntaRNA jobs complete to produce final chimeras.txt, singletons.txt, and interactions.txt.

## Troubleshooting

### Python Dependencies

If you encounter import errors:

```bash
# Reinstall dependencies
pip install --upgrade biopython bcbio-gff pysam mpire

# Install optional dependencies for better performance
pip install psutil
```

### External Tools Not Found

Ensure tools are in your PATH:

```bash
# Check if tools are accessible
which bwa
which samtools
which bedtools

# If not found, ensure conda environment is activated
conda activate chira
```

### Permission Errors

If you get permission errors during installation:

```bash
# Install to user directory
pip install --user -e .

# Or use conda environment (recommended)
conda activate chira
pip install -e .
```

## Using ChiRA in LSF Job Scripts

After installing ChiRA with `pip install -e .`, the scripts are available as commands in your conda environment.

### Example bsub Job Script

```bash
#!/bin/bash
#BSUB -n 8
#BSUB -R rusage[mem=16000]
#BSUB -W 240:00
#BSUB -J "chira_map[5]"
#BSUB -q long
#BSUB -o logs/out.%J.%I.txt
#BSUB -e logs/err.%J.%I.txt

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate chira

# Load modules (if available)
module load samtools
module load bwa

# Use chira_map.py directly (available after pip install)
chira_map.py --aligner bwa \
   -i input.fasta \
   -o output_dir \
   --chunk_fasta 10 \
   --use_batchtools \
   --batchtools_queue long \
   --batchtools_cores 8 \
   --batchtools_memory 8GB \
   --index1 index1.fa \
   --index2 index2.fa \
   -s fw -p 8 -co 2
```

### Important Notes for Batchtools

1. **R and batchtools must be available**: The conda environment needs R, batchtools, and jsonlite:
   ```bash
   conda activate chira
   conda install -c conda-forge r-base r-batchtools r-jsonlite
   ```
   - `jsonlite` is required for JSON configuration file parsing
   - Usually installed automatically as a dependency of `batchtools`

2. **Conda environment path**: Use `--batchtools_conda_env` to specify the full path (absolute path recommended):
   ```bash
   --batchtools_conda_env /path/to/miniconda3/envs/chira
   ```
   - Relative paths are automatically converted to absolute paths
   - Absolute paths ensure cluster jobs can find the environment

3. **Template file path**: The `--batchtools_template` parameter accepts:
   - Absolute path to custom template file (recommended)
   - Relative path (automatically converted to absolute)
   - Built-in `"lsf-simple"` template (not recommended for most use cases)

4. **Path handling**: All file paths are automatically converted to absolute paths:
   - Registry directory, chunk directory, Python script paths
   - Template file paths (if custom template is used)
   - Chunk FASTA file paths
   - BWA reference index paths
   - This ensures cluster jobs can correctly resolve all paths regardless of working directory

5. **Scripts location**: After `pip install -e .`, scripts are in your PATH:
   - `chira_map.py` - Direct command (recommended)
   - `python -m chira_map` - Alternative method
   - Full path not needed if installed

See `EXAMPLE_bsub_JOB.sh` for a complete working example.

## Uninstalling

```bash
# Deactivate environment
conda deactivate

# Remove environment
conda env remove -n chira

# Or uninstall package only
pip uninstall chira
```

## Dependencies Summary

### Required Python Packages (auto-installed)
- `biopython` - FASTA/sequence parsing
- `bcbio-gff` - GFF/GTF annotation parsing
- `pysam` - BAM file manipulation

### Optional Python Packages (recommended)
- `psutil` - Memory optimization and I/O performance
- `requests` - Downloading Ensembl files
- `pyliftover` - Coordinate liftover

### Required External Tools (install via conda)
- `bwa` - Sequence alignment
- `samtools` - BAM file processing
- `bedtools` - Genomic interval operations

### Optional External Tools
- `intarna` - RNA-RNA interaction prediction
- `gffread` - GFF utilities

### Optional R Packages (for batchtools support)
- `batchtools` - HPC cluster job submission (required for `--use_batchtools`)
- `jsonlite` - JSON file parsing (required for batchtools configuration)

## Notes

- **Python Version**: Requires Python >= 3.6 (Python 3.9+ recommended)
- **Conda Channel**: Use `bioconda` channel for bioinformatics tools
- **Editable Install**: Use `pip install -e .` during development to see code changes immediately
- **Batchtools**: Required only if using `--use_batchtools` in `chira_map.py` or `--hybridize --use_batchtools` in `chira_extract.py`
  - Requires R with `batchtools` and `jsonlite` packages
  - All file paths are automatically converted to absolute paths for cluster job execution
  - See `BATCHTOOLS_USAGE.md` for detailed usage instructions

