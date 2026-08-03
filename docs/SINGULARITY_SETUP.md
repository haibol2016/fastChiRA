# Singularity/Apptainer Setup and Usage Guide for ChiRA

This guide provides comprehensive instructions for using ChiRA with Singularity/Apptainer containers, including installation, best practices, and troubleshooting.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Installation](#installation)
3. [Converting Docker Images to Singularity](#converting-docker-images-to-singularity)
4. [Basic Usage](#basic-usage)
5. [Environment Setup](#environment-setup)
6. [Best Practices](#best-practices)
7. [Environment Isolation](#environment-isolation)
8. [Troubleshooting](#troubleshooting)
9. [Reference](#reference)

---

## Quick Start

**For Linux systems (most common in bioinformatics):**

```bash
# 1. Install Singularity/Apptainer
sudo apt-get install -y apptainer  # Ubuntu/Debian
# OR
sudo yum install -y singularity    # CentOS/RHEL

# 2. Pull ChiRA image
singularity pull docker://docker.io/nemat1976/chiraplus:v0.0.1

# 3. Run ChiRA command (entrypoint script runs automatically)
singularity exec \
  -B /path/to/data:/app/data \
  -B /path/to/output:/app/output \
  chira_latest.sif \
  chira_collapse.py -i data/input.fastq -o output/collapsed.fasta
```

**Recommended command template (entrypoint script handles environment automatically):**

```bash
singularity exec \
  -B /host/data:/app/data:ro \   # Mount input data (read-only)
  -B /host/output:/app/output \  # Mount output directory
  chira_latest.sif \
  chira_collapse.py [options]
```

**Note:** The entrypoint script automatically sets up PATH, PYTHONPATH, and conda environment. You don't need to specify `--no-home --cleanenv` unless you have specific isolation requirements (see [Environment Setup](#environment-setup) and [Environment Isolation](#environment-isolation) sections for details).

---

## Installation

### Linux Systems (Recommended for Production)

Most HPC clusters and Linux servers already have Singularity/Apptainer installed. If not:

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install -y apptainer
```

**CentOS/RHEL:**
```bash
sudo yum install -y epel-release
sudo yum install -y singularity
```

**Verify installation:**
```bash
singularity --version
# OR
apptainer --version
```

### macOS (Development/Testing Only)

Singularity requires Linux kernel features and cannot run natively on macOS. Options:

#### Option 1: Docker Desktop (Easiest for Testing)

```bash
# Install Docker Desktop from: https://www.docker.com/products/docker-desktop
# Then use Docker directly (no Singularity needed for testing)
docker run --rm docker.io/nemat1976/chiraplus:v0.0.1 chira_collapse.py --help
```

#### Option 2: Linux VM (Lima - Recommended for Singularity Testing)

```bash
# Install Lima
brew install lima

# Start Ubuntu VM
limactl start template://ubuntu
limactl shell ubuntu

# Inside VM, install Apptainer
sudo apt-get update
sudo apt-get install -y apptainer
```

---

## Converting Docker Images to Singularity

### Method 1: Pull from Docker Registry (Recommended)

```bash
# Pull directly from Docker Hub
singularity pull docker://docker.io/nemat1976/chiraplus:v0.0.1

# Result: chira_latest.sif (Singularity Image Format)
```

### Method 2: Build from Local Docker Image

```bash
# Save Docker image to tar file
docker save docker.io/nemat1976/chiraplus:v0.0.1 -o chira.tar

# Build Singularity image from tar
singularity build chira.sif docker-archive://chira.tar
```

### Method 3: Build from Dockerfile (Requires Docker Daemon)

```bash
# Build Singularity image directly from Dockerfile (requires Docker daemon)
# First build Docker image: docker build -t docker.io/nemat1976/chiraplus:v0.0.1 .
singularity build chira.sif docker-daemon://docker.io/nemat1976/chiraplus:v0.0.1
```

**Note:** Pull once, reuse the `.sif` file for all subsequent runs. This is much faster than pulling from Docker each time.

---

## Basic Usage

### Running Commands

**Single command execution (entrypoint script runs automatically):**
```bash
singularity exec chira_latest.sif chira_collapse.py --help
```

**With data mounts:**
```bash
singularity exec \
  -B /host/data:/app/data \
  -B /host/output:/app/output \
  chira_latest.sif \
  chira_collapse.py -i data/input.fastq -o output/collapsed.fasta
```

**Note:** The entrypoint script automatically sets up the environment (PATH, PYTHONPATH, conda). You don't need to specify isolation flags unless you have specific requirements (see [Environment Setup](#environment-setup) and [Environment Isolation](#environment-isolation) sections).

**Interactive shell (for debugging):**
```bash
singularity shell \
  -B /host/data:/app/data \
  chira_latest.sif

# Inside container:
# - Check PATH: echo $PATH
# - Test tools: which python, which bwa
# - Test imports: python -c "import Bio"
```

### Verifying Installation

```bash
# Check image information
singularity inspect chira_latest.sif

# Test Python (entrypoint script handles environment automatically)
singularity exec chira_latest.sif python --version

# Test ChiRA scripts (entrypoint script runs automatically)
singularity exec chira_latest.sif chira_collapse.py --help

# Verify tools are accessible (entrypoint script handles PATH)
singularity exec chira_latest.sif which python bwa samtools

# Test Python package imports (entrypoint script handles PYTHONPATH)
singularity exec chira_latest.sif \
  python -c "import Bio; import BCBio; import pysam; print('All packages OK')"
```

---

## Best Practices

### 1. Image Management

**Pull once, reuse many times:**
```bash
# Pull the image once
singularity pull docker://docker.io/nemat1976/chiraplus:v0.0.1

# Reuse the .sif file for all runs (much faster)
singularity exec chira_latest.sif command1
singularity exec chira_latest.sif command2
```

**Set cache directory (useful on shared systems):**
```bash
export SINGULARITY_CACHEDIR=/path/to/cache
singularity pull docker://docker.io/nemat1976/chiraplus:v0.0.1
```

### 2. Data Mounting

**Use explicit bind mounts:**
```bash
# Good: Explicit mounts with clear paths
singularity exec \
  -B /project/data:/app/data:ro \    # Read-only input
  -B /project/output:/app/output \   # Read-write output
  -B /project/tmp:/tmp \             # Custom tmp location
  chira_latest.sif command
```

**Use absolute paths:**
```bash
# Good: Absolute paths
-B /home/user/project/data:/app/data

# Avoid: Relative paths (can be confusing)
-B ./data:/app/data
```

**Mount only what you need:**
```bash
# Good: Specific directories
-B /project/data:/app/data

# Avoid: Entire home directory (security risk)
# NOT: -B /home/user:/home/user
```

### 3. Performance Optimization

**Parallel execution:**
```bash
# Singularity images can run in parallel
singularity exec chira_latest.sif command1 &
singularity exec chira_latest.sif command2 &
wait
```

**Use local .sif files:**
```bash
# Faster: Use local .sif file
singularity exec chira_latest.sif command

# Slower: Pull from Docker each time (don't do this)
singularity exec docker://docker.io/nemat1976/chiraplus:v0.0.1 command
```

### 4. Security

**Run as non-root:**
```bash
# ChiRA Docker image creates a non-root user
# Singularity runs as that user by default
singularity exec chira_latest.sif command  # Safe
```

**Read-only mounts for input:**
```bash
# Mount input data as read-only
-B /data/input:/app/data:ro \
-B /data/output:/app/output \
```

### 5. Complete Workflow Example

```bash
#!/bin/bash
# Production-ready ChiRA workflow script

set -e  # Exit on error

# Configuration
IMAGE="chira_latest.sif"
DATA_DIR="/project/data"
OUTPUT_DIR="/project/output"
TMP_DIR="/project/tmp"

# Verify image exists
if [ ! -f "$IMAGE" ]; then
    echo "Error: Image $IMAGE not found"
    exit 1
fi

# Run ChiRA with best practices
# Entrypoint script runs automatically to set up environment
# Isolation flags are optional - add if you have host environment conflicts
singularity exec \
  --no-home \                    # Optional: Isolate from host $HOME
  --cleanenv \                   # Optional: Isolate from host environment
  -B "$DATA_DIR:/app/data:ro" \  # Read-only input
  -B "$OUTPUT_DIR:/app/output" \ # Read-write output
  -B "$TMP_DIR:/tmp" \           # Custom tmp location
  "$IMAGE" \
  chira_collapse.py \
    -i data/input.fastq \
    -o output/collapsed.fasta

# Check exit status
if [ $? -eq 0 ]; then
    echo "Success: ChiRA completed successfully"
else
    echo "Error: ChiRA failed"
    exit 1
fi
```

---

## Environment Setup

### Entrypoint Script (Recommended - Automatic Environment Setup)

**The entrypoint script is the recommended approach** because it automatically sets up the correct environment without requiring you to know which environment variables need to be activated. It handles all the complexity for you.

**Why the entrypoint script is important:**

The entrypoint script is a **reliable solution for PATH issues** in Singularity containers pulled from Docker images. It **actively sets up the correct environment** inside the container, ensuring all tools and Python packages are accessible regardless of how Singularity handles PATH.

**Key benefits:**
- **No need to know environment variables** - The script handles everything automatically
- **Explicitly activates conda environment** - Ensures conda's bin directory is in PATH
- **Sets PATH correctly** - Guarantees `/usr/local/bin` and `/opt/conda/bin` are accessible
- **Sets PYTHONPATH** - Ensures Python can find ChiRA scripts and packages
- **Solves persistent PATH problems** - Many users report this completely resolves long-standing PATH issues

**How to use it:**

The ChiRA Docker image includes an entrypoint script that sets up the conda environment. The script is automatically executed by Singularity, so you typically don't need to call it explicitly:

```bash
# Simple usage - entrypoint script runs automatically
singularity exec chira_latest.sif \
  chira_collapse.py --help
```

If you need to explicitly call the entrypoint script (e.g., for troubleshooting or when automatic execution doesn't work):

```bash
# Explicitly call the entrypoint script
singularity exec chira_latest.sif \
  /usr/local/bin/_entrypoint.sh \
  chira_collapse.py --help
```

**Entrypoint script contents:**

The script (`/usr/local/bin/docker-entrypoint.sh`) is automatically generated in the Docker image:

```bash
#!/bin/bash
# Docker entrypoint script for ChiRA
# Sets up conda environment and ensures PATH is correct

# Activate conda base environment
if [ -f /opt/conda/etc/profile.d/conda.sh ]; then
    source /opt/conda/etc/profile.d/conda.sh
    conda activate base
fi

# Ensure /usr/local/bin is in PATH (where all tools are symlinked)
export PATH="/usr/local/bin:/opt/conda/bin:${PATH}"

# Set PYTHONPATH for ChiRA scripts
export PYTHONPATH="/app:${PYTHONPATH}"

# Execute command or start bash
if [ $# -eq 0 ]; then
    exec /bin/bash
else
    exec "$@"
fi
```

**What the entrypoint script does:**
- Activates conda base environment
- Sets PATH to include `/usr/local/bin` and `/opt/conda/bin`
- Sets PYTHONPATH to `/app`
- Executes your command

**Why it works:**
The entrypoint script **proactively sets up the environment** by explicitly:
1. Sources conda's activation script
2. Activates the base conda environment
3. Sets PATH to include all necessary directories
4. Sets PYTHONPATH for Python package discovery

This approach ensures the environment is correct **regardless of how Singularity initializes PATH**, making it a robust solution for PATH-related problems.

---

## Environment Isolation

### Understanding Isolation Flags

Isolation flags control how the host environment interacts with the container. These are useful when you need to prevent host environment variables or filesystem mounts from interfering with the container.

| Flag | What It Does | Use When |
|------|--------------|----------|
| `--no-home` | Prevents mounting host `$HOME` | Host has conflicting software in `~/.local/` or `~/.conda/` |
| `--cleanenv` | Prevents host environment variables from entering container | Host `PATH`, `PYTHONPATH`, or `LD_LIBRARY_PATH` conflicts |
| `--containall` | Maximum isolation (no `$HOME`, `/tmp`, or current directory mounts) | Need complete filesystem isolation |

**For most cases, use `--no-home --cleanenv`:**

```bash
singularity exec \
  --no-home \      # Prevent host $HOME from shadowing container tools
  --cleanenv \     # Prevent host PATH/PYTHONPATH from interfering
  chira_latest.sif \
  chira_collapse.py [options]
```

**Why this works:**
- `--no-home`: Prevents host `$HOME` (which may contain `~/.local/bin/python`, `~/.conda/`, etc.) from shadowing container tools
- `--cleanenv`: Prevents host environment variables from interfering with container's PATH and PYTHONPATH

**Combining with entrypoint script:**

You can combine the entrypoint script with isolation flags for maximum reliability:

```bash
singularity exec \
  --no-home --cleanenv \
  chira_latest.sif \
  /usr/local/bin/_entrypoint.sh \
  chira_collapse.py [options]
```

### Decision Guide

| Situation | Recommended Solution |
|-----------|---------------------|
| **Default / Don't know what to use** | **Use entrypoint script (automatic, handles everything)** |
| **Persistent PATH issues / Tools not found** | **Use entrypoint script (often definitive solution)** |
| **Want simplest solution** | **Use entrypoint script (no need to know environment variables)** |
| Host has Python/conda in `$HOME` | Use `--no-home` (or combine with entrypoint script) |
| Host `PATH`/`PYTHONPATH` conflicts | Use `--cleanenv` (or combine with entrypoint script) |
| Both issues | Use `--no-home --cleanenv` (or combine with entrypoint script) |
| Need maximum isolation | Use `--containall --cleanenv` (or combine with entrypoint script) |
| Want maximum reliability | Combine both approaches (entrypoint script + isolation flags) |

---

## Troubleshooting

### Problem: Tools Not Found

**Symptoms:** Commands like `python`, `bwa`, or `samtools` are not found even though they're installed in the container.

**Solutions (in order of recommendation):**

1. **Use entrypoint script (Recommended - primary solution):**
   ```bash
   # The entrypoint script explicitly sets up PATH and conda environment
   # This is the recommended approach - no need to know which environment variables to set
   # Usually works automatically:
   singularity exec chira_latest.sif command
   
   # Or explicitly call it if needed:
   singularity exec chira_latest.sif /usr/local/bin/_entrypoint.sh command
   ```

2. **Use isolation flags (Alternative approach):**
   ```bash
   # Use isolation flags if you prefer or need additional isolation
   singularity exec --no-home --cleanenv chira_latest.sif command
   ```

3. **Combine both approaches (Maximum reliability):**
   ```bash
   # Combine entrypoint script with isolation flags for maximum reliability
   singularity exec --no-home --cleanenv chira_latest.sif /usr/local/bin/_entrypoint.sh command
   ```

4. **Check PATH inside container:**
   ```bash
   singularity shell chira_latest.sif
   # Inside: echo $PATH
   # Inside: ls /usr/local/bin | grep -E "(python|bwa|samtools)"
   ```

5. **Override PATH:**
   ```bash
   export SINGULARITYENV_PATH="/usr/local/bin:/opt/conda/bin:$PATH"
   singularity exec chira_latest.sif command
   ```

6. **Use absolute paths:**
   ```bash
   singularity exec chira_latest.sif /usr/local/bin/python /usr/local/bin/chira_collapse.py --help
   ```

### Problem: Python Packages Not Found

**Symptoms:** `ImportError` when trying to import `Bio`, `BCBio`, or `pysam`.

**Solutions:**

1. **Use entrypoint script (Recommended - ensures Python packages are found):**
   ```bash
   # The entrypoint script sets PYTHONPATH and activates conda environment
   # This ensures Python can find all installed packages
   singularity exec chira_latest.sif /usr/local/bin/_entrypoint.sh python -c "import Bio"
   ```

2. **Use `--no-home --cleanenv`:**
   ```bash
   singularity exec --no-home --cleanenv chira_latest.sif python -c "import Bio"
   ```

3. **Check Python package location:**
   ```bash
   singularity exec chira_latest.sif python -c "import Bio; print(Bio.__file__)"
   ```

### Problem: Permission Denied

**Symptoms:** Cannot write to mounted output directory.

**Solutions:**

1. **Check directory permissions:**
   ```bash
   ls -ld /path/to/output
   chmod 755 /path/to/output  # If needed
   ```

2. **Use `--fakeroot` (if available):**
   ```bash
   singularity exec --fakeroot chira_latest.sif command
   ```

3. **Run as appropriate user:**
   ```bash
   # Container runs as non-root user by default
   # Ensure output directory is writable by that user
   ```

### Problem: Slow Performance

**Symptoms:** Container startup or execution is slow.

**Solutions:**

1. **Use local .sif file instead of pulling:**
   ```bash
   # Pull once
   singularity pull docker://docker.io/nemat1976/chiraplus:v0.0.1
   
   # Reuse .sif file
   singularity exec chira_latest.sif command
   ```

2. **Set cache directory:**
   ```bash
   export SINGULARITY_CACHEDIR=/fast/storage/cache
   ```

3. **Use faster storage for image:**
   ```bash
   # Store .sif file on fast storage (SSD, local disk)
   # Not on network filesystem
   ```

### Debugging Checklist

If something doesn't work, systematically check:

1. **PATH inside container (entrypoint script sets this automatically):**
   ```bash
   # With entrypoint script (recommended)
   singularity exec chira_latest.sif /usr/local/bin/_entrypoint.sh bash -c 'echo $PATH'
   
   # Or check without entrypoint script
   singularity exec chira_latest.sif echo $PATH
   ```

2. **Tools accessible (entrypoint script ensures PATH is correct):**
   ```bash
   # With entrypoint script (recommended)
   singularity exec chira_latest.sif which python bwa samtools
   
   # Or explicitly call entrypoint script
   singularity exec chira_latest.sif /usr/local/bin/_entrypoint.sh which python bwa samtools
   ```

3. **Python packages (entrypoint script sets PYTHONPATH):**
   ```bash
   # With entrypoint script (recommended)
   singularity exec chira_latest.sif \
     python -c "import Bio; import BCBio; import pysam; print('OK')"
   
   # Or explicitly call entrypoint script
   singularity exec chira_latest.sif /usr/local/bin/_entrypoint.sh \
     python -c "import Bio; import BCBio; import pysam; print('OK')"
   ```

4. **Mounted directories:**
   ```bash
   singularity shell -B /data:/app/data chira_latest.sif
   # Inside: ls /app/data
   ```

5. **Environment variables (entrypoint script sets these):**
   ```bash
   # With entrypoint script (recommended)
   singularity exec chira_latest.sif env | grep -E "(PATH|PYTHONPATH)"
   
   # Or explicitly call entrypoint script
   singularity exec chira_latest.sif /usr/local/bin/_entrypoint.sh env | grep -E "(PATH|PYTHONPATH)"
   ```

---

## Reference

### Command Templates

**Template 1: Entrypoint script (Recommended - automatic environment setup):**
```bash
singularity exec \
  -B /host/data:/app/data:ro \
  -B /host/output:/app/output \
  chira_latest.sif \
  chira_collapse.py [options]
```

**Template 2: Isolation flags (Alternative - when you need host isolation):**
```bash
singularity exec \
  -B /host/data:/app/data:ro \
  -B /host/output:/app/output \
  chira_latest.sif \
  /usr/local/bin/_entrypoint.sh \
  chira_collapse.py [options]
```

**Template 3: Maximum reliability (entrypoint script + isolation flags):**
```bash
singularity exec \
  --no-home \
  --cleanenv \
  -B /host/data:/app/data:ro \
  -B /host/output:/app/output \
  chira_latest.sif \
  /usr/local/bin/_entrypoint.sh \
  chira_collapse.py [options]
```

### Common Commands

| Task | Command |
|------|---------|
| Pull image | `singularity pull docker://docker.io/nemat1976/chiraplus:v0.0.1` |
| Inspect image | `singularity inspect chira_latest.sif` |
| Run command | `singularity exec chira_latest.sif command` |
| Interactive shell | `singularity shell chira_latest.sif` |
| Check version | `singularity --version` |

### Isolation Flags Summary

| Flag | Purpose | When to Use |
|------|---------|-------------|
| `--no-home` | Prevent `$HOME` mount | Host has conflicting software in home directory |
| `--cleanenv` | Prevent env var inheritance | Host PATH/PYTHONPATH conflicts |
| `--containall` | Maximum filesystem isolation | Need complete isolation |
| `--fakeroot` | Run as fake root | Need root privileges (if available) |

### Additional Running Options

These options are useful for specialized use cases but not required for most ChiRA workflows:

**GPU Support (`--nv`):**
```bash
# Enable NVIDIA GPU support (if using GPU-accelerated tools)
singularity exec --nv chira_latest.sif command
```

**Writable Temporary Filesystem (`--writable-tmpfs`):**
```bash
# Create writable temporary filesystem in memory (for tools that need to write)
singularity exec --writable-tmpfs chira_latest.sif command
```

**Set Working Directory (`--pwd`):**
```bash
# Set working directory inside container
singularity exec --pwd /app/data chira_latest.sif command
```

**Set Environment Variables (`--env`):**
```bash
# Set environment variables (alternative to SINGULARITYENV_ prefix)
singularity exec --env VAR=value chira_latest.sif command
```

**Scratch Directories (`--scratch`):**
```bash
# Create temporary scratch directories (useful for HPC systems)
singularity exec --scratch /tmp/scratch1 --scratch /tmp/scratch2 chira_latest.sif command
```

**Overlay Filesystem (`--overlay`):**
```bash
# Use overlay filesystem for persistent writable layer
singularity exec --overlay overlay.img chira_latest.sif command
```

**Network Options (`--network`):**
```bash
# Control network access (useful for security or network isolation)
singularity exec --network none chira_latest.sif command  # No network
singularity exec --network bridge chira_latest.sif command  # Bridge network
```

---

## Additional Resources

- **Singularity Documentation:** https://docs.sylabs.io/
- **Apptainer Documentation:** https://apptainer.org/docs/
- **ChiRA Documentation:** See `README.md` for ChiRA-specific usage

---

**Last Updated:** January 2026
