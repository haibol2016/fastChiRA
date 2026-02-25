# Using batchtools for LSF Cluster Parallel Processing

## Overview

ChiRA supports R batchtools as a backend for submitting chunk processing jobs to LSF clusters. This enables **true parallel computing** by distributing chunks across multiple cluster nodes, each running independently. This is simpler and more reliable than Dask for LSF environments.

**Key Benefits:**
- **Parallel Processing**: Each chunk runs as an independent LSF job on different cluster nodes
- **Scalability**: Process hundreds of chunks simultaneously across the cluster
- **Resource Control**: Specify cores, memory, and walltime per job
- **Visibility**: All jobs visible in `bjobs` for easy monitoring
- **Reliability**: Each job is independent - failures don't affect other jobs

## Prerequisites

1. **R installed** (from [CRAN](https://cran.r-project.org/) or via conda):
   ```bash
   conda install -c conda-forge r-base
   ```

2. **R packages** - `batchtools` and `jsonlite`:
   ```r
   install.packages(c("batchtools", "jsonlite"))
   ```
   Or via conda:
   ```bash
   conda install -c conda-forge r-batchtools r-jsonlite
   ```
   - **batchtools**: Submitting chunk-based batch jobs to HPC cluster schedulers (LSF, SLURM, SGE, etc.)
   - **jsonlite**: Parsing JSON configuration files for batchtools job submission
   - **Note**: `jsonlite` is usually installed automatically as a dependency of `batchtools`.

3. **LSF scheduler** available on your cluster

4. **Shared filesystem** accessible to all compute nodes

5. **Conda environment** (optional but recommended) with ChiRA dependencies

6. **Absolute paths**: All file paths are automatically converted to absolute paths for cluster job execution
   - This ensures jobs running on different cluster nodes can correctly resolve all file paths
   - No manual path conversion needed - handled automatically by the code

## Quick Start: Running with Batchtools

### Basic Example

```bash
python chira_map.py \
  -i input.fasta \
  -o output_dir \
  --index1 index1.fa \
  --index2 index2.fa \
  --chunk_fasta 20 \
  --use_batchtools \
  --batchtools_queue long \
  --batchtools_cores 8 \
  --batchtools_memory 8GB
```

**What happens:**
1. FASTA file is split into 20 chunks
2. 20 independent LSF jobs are submitted (one per chunk)
3. Each job runs on a different cluster node with 8 cores and 8GB memory
4. Jobs run in parallel across the cluster
5. Results are automatically collected and merged

### All batchtools Options

- `--use_batchtools`: **REQUIRED** - Enable batchtools mode for LSF cluster parallel processing
- `--batchtools_queue`: LSF queue name (default: `long`)
- `--batchtools_cores`: **Cores per LSF job** (default: 8)
  - Each chunk job runs with this many cores on a cluster node
  - Example: `--batchtools_cores 16` → each job uses 16 cores
  - **Independent of `--processes`** (which applies to main job only)
- `--batchtools_memory`: **Total memory per LSF job** (automatically converted to per-core for LSF), e.g., `16GB` (default: auto-calculated)
  - **IMPORTANT**: You specify TOTAL memory, but LSF `rusage[mem=...]` is PER CORE
  - Code automatically converts: `total_memory ÷ cores = per_core_memory`
  - Example: `--batchtools_cores 8 --batchtools_memory 16GB` → LSF gets `rusage[mem=2GB]` (16GB ÷ 8 = 2GB per core)
  - Auto-calculation: `cores × 0.5GB` total (e.g., 8 cores → 4GB total → 0.5GB per core)
- `--batchtools_walltime`: Walltime limit per job, e.g., `240:00` (default: `240:00`)
- `--batchtools_max_parallel`: **Max concurrent running jobs** (default: unlimited, all submitted at once)
  - Limits how many jobs run simultaneously
  - Example: `--batchtools_max_parallel 2` → only 2 jobs run at a time, others wait
  - Useful if cluster has job limits per user
- `--batchtools_conda_env`: Conda environment path (default: auto-detect from `CONDA_DEFAULT_ENV`)
- `--batchtools_template`: LSF template file (default: `lsf_custom.tmpl` in ChiRA directory)
  - Can be an absolute path to a custom template file
  - Can be the built-in `"lsf-simple"` template (not recommended for most use cases)
  - Relative paths are automatically resolved to absolute paths
  - Uses proven template based on InPAS implementation
  - Only change if you need custom LSF directives

**Minimizing walltime (IntaRNA / chira_extract):** Leave `--batchtools_max_parallel` unset so all chunk jobs are submitted at once; use high `-p` and `--batchtools_cores` (or `--batchtools_cores 1` to maximize concurrent jobs).

### Example: Full Command with Parallel Processing

```bash
python chira_map.py \
  -i reads.fasta \
  -o results \
  --index1 index1.fa \
  --index2 index2.fa \
  --chunk_fasta 20 \
  --use_batchtools \
  --batchtools_queue long \
  --batchtools_cores 8 \
  --batchtools_memory 8GB \
  --batchtools_walltime 240:00 \
  --batchtools_max_parallel 10 \
  --batchtools_conda_env ~/miniconda3/envs/chira
```

**This command:**
- Splits `reads.fasta` into 20 chunks
- Submits 20 LSF jobs to the `long` queue
- Each job uses 8 cores and 8GB total memory (1GB per core)
- Limits to 10 concurrent jobs (if cluster has limits)
- All 20 chunks process in parallel across cluster nodes

### Example: Limiting Concurrent Jobs

If your cluster limits concurrent jobs per user, use `--batchtools_max_parallel`:

```bash
python chira_map.py \
  -i large.fasta \
  -o results \
  --index1 index1.fa \
  --index2 index2.fa \
  --chunk_fasta 50 \
  --use_batchtools \
  --batchtools_queue long \
  --batchtools_cores 16 \
  --batchtools_memory 16GB \
  --batchtools_max_parallel 2
```

**This ensures:**
- Only 2 jobs run at a time (respects cluster limits)
- All 50 chunks will be processed sequentially in batches of 2
- Each batch waits for completion before submitting the next

## How Parallel Processing Works

1. **Chunk Preparation**: Python splits FASTA file into N chunks (specified by `--chunk_fasta`)
2. **Path Normalization**: All file paths are automatically converted to absolute paths:
   - Registry directory, chunk directory, Python script path
   - Template file path (if custom template is used)
   - Individual chunk FASTA file paths
   - BWA reference index paths
   - This ensures cluster jobs can correctly resolve all paths regardless of working directory
3. **Configuration Files**: Python creates JSON configuration files with all absolute paths:
   - `config.json`: Contains all job parameters and absolute paths
   - `chunks.json`: Contains list of all chunks with absolute paths (one file for all chunks - this is the intended batchtools design)
4. **Job Submission**: Python calls R script (`submit_chunks_batchtools.R`) that uses batchtools to submit N independent LSF jobs
5. **Parallel Execution**: Each LSF job runs on a **different cluster node**:
   - Job 1 processes chunk 1 on node A
   - Job 2 processes chunk 2 on node B
   - Job 3 processes chunk 3 on node C
   - ... and so on
6. **Resource Allocation**: Each job gets:
   - `--batchtools_cores` cores
   - `--batchtools_memory` total memory
   - Runs independently with its own resources
7. **Result Collection**: Python monitors completion and collects BAM files from each chunk
8. **Merging**: Python merges all chunk BAMs into final BAM files

**Key Point**: This is **true parallel processing** - chunks run simultaneously on different nodes, not sequentially on one node.

**Important**: The `chunks.json` file contains information for **all chunks** in a single file. This is the intended design for batchtools - the `batchMap()` function reads this single file and distributes tasks to individual jobs. Each job receives its chunk information via batchtools' task system.

## Monitoring Jobs

### View Submitted Jobs

```bash
# View all your jobs
bjobs -u $USER

# View jobs with specific prefix
bjobs -J chira_bt_*
```

### Check Batchtools Registry

The registry directory is created in your output directory:
```
output_dir/batchtools_registry_<timestamp>/
```

You can check job status in R:
```r
library(batchtools)
reg <- loadRegistry("output_dir/batchtools_registry_<timestamp>")
getStatus(reg)
```

### Check Configuration Files

The configuration files are created in the registry directory:
- `config.json`: Contains all job parameters and absolute paths
- `chunks.json`: Contains list of all chunks with absolute paths (one file for all chunks)

You can inspect these files to verify paths are correct:
```bash
cat output_dir/batchtools_registry_<timestamp>/config.json
cat output_dir/batchtools_registry_<timestamp>/chunks.json
```

**Note**: All paths in these files are absolute paths to ensure cluster jobs can correctly resolve files regardless of working directory.

## Advantages Over Dask

1. **Simpler**: Direct LSF job submission, no worker management
2. **More Reliable**: Each chunk is an independent LSF job
3. **Better Visibility**: All jobs visible in `bjobs`
4. **Easier Debugging**: Each job has its own log file
5. **No Network Overhead**: No task serialization/deserialization

## Troubleshooting

### R/batchtools Not Found

Ensure R is in your PATH:
```bash
which Rscript
```

Install R if needed:
```bash
conda install -c conda-forge r-base
```

Install batchtools and jsonlite R packages if needed:
```r
install.packages(c("batchtools", "jsonlite"))
```
Or via conda:
```bash
conda install -c conda-forge r-batchtools r-jsonlite
```

**Note**: These are R packages (not Python packages). They must be installed in your R environment.

### JSON Parsing Errors

If you encounter JSON parsing errors during batchtools submission:

1. **Check file encoding**: All JSON files are written with UTF-8 encoding
2. **Check file paths**: All paths should be absolute (handled automatically)
3. **Check file content**: Inspect the JSON files in the registry directory:
   ```bash
   cat output_dir/batchtools_registry_<timestamp>/config.json
   cat output_dir/batchtools_registry_<timestamp>/chunks.json
   ```
4. **Check R error messages**: The R script provides detailed error messages with file content previews if JSON parsing fails

Common causes:
- Invalid characters in file paths (should be handled automatically by path normalization)
- File system issues (permissions, disk space)
- Corrupted JSON files (check file integrity)

### Jobs Not Starting

Check LSF queue status:
```bash
bqueues
```

Check job requirements match your LSF configuration:
- `--batchtools_cores` should match `#BSUB -n` in your job script
- `--batchtools_memory` should match `#BSUB -R rusage[mem=...]`

### Conda Environment Not Found

Specify conda environment explicitly (use absolute path):
```bash
--batchtools_conda_env /path/to/miniconda3/envs/chira
```

Or ensure `CONDA_DEFAULT_ENV` is set:
```bash
export CONDA_DEFAULT_ENV=chira_env
```

**Note**: The conda environment path is automatically converted to an absolute path if a relative path is provided.

### nodename Parameter

The `nodename = "localhost"` parameter in `submit_chunks_batchtools.R` (line 51) specifies the hostname where the LSF scheduler runs.

**For most clusters: Keep as `"localhost"`**
- You submit jobs from the login node where LSF commands are available
- LSF scheduler runs on the same node
- This is the default and works for 99% of LSF clusters

**Only change if:**
- Your cluster documentation specifies a different submission host
- LSF scheduler runs on a different machine
- You get connection errors when submitting jobs

**To auto-detect your hostname:**
```r
nodename = Sys.info()["nodename"]
```

**To manually specify:**
```r
nodename = "your-lsf-hostname"
```

## Template Files

**Default Template**: ChiRA uses `lsf_custom.tmpl` by default, which is based on the proven InPAS implementation. This template is included with ChiRA and works for most LSF clusters.

**The default template (`lsf_custom.tmpl`) handles:**
- Queue specification (`#BSUB -q`)
- Resource requirements (`#BSUB -n`, `#BSUB -R rusage[mem=...]`)
- Walltime (`#BSUB -W`)
- Job name (`#BSUB -J`)
- Log files (`#BSUB -o`, `#BSUB -e`)
- Host spanning (`#BSUB -R "span[hosts=1]"`)
- Threading environment variables (OMP_NUM_THREADS, etc.)

**Template Path Handling:**
- **Relative paths**: Automatically resolved to absolute paths
- **Absolute paths**: Used as-is
- **Built-in template**: Use `"lsf-simple"` to use batchtools' built-in simple LSF template (not recommended for most use cases)

**To use a different template:**
```bash
# Absolute path (recommended)
--batchtools_template /path/to/your/custom.tmpl

# Relative path (automatically converted to absolute)
--batchtools_template ./my_custom.tmpl
```

**To use built-in "lsf-simple" (not recommended):**
```bash
--batchtools_template lsf-simple
```

**Only create a custom template if:**
- You need to load specific modules (e.g., `module load python/3.9`)
- You need additional custom LSF directives
- Your cluster has non-standard LSF configuration

**Important**: Template file paths are automatically converted to absolute paths to ensure cluster jobs can find the template file regardless of working directory.

See `lsf_custom.tmpl` for the default template structure.

## chira_extract.py Hybridization with Batchtools

When using `--hybridize` and `--use_batchtools` in `chira_extract.py`, IntaRNA jobs are submitted via a separate R script:

- **submit_intarna_batchtools.R**: Submits one IntaRNA job per chunk; each job runs IntaRNA once per locus pair (only real chimeric pairs; all-vs-all is not used). Reads `config.json` and `jobs.json` from `batchtools_work/`, creates a registry, and submits LSF jobs. Uses the same LSF template and resource options as chunk mapping (`--batchtools_*`). The script **waits for all jobs to complete** before returning (manual polling with `getJobTable` and status counting; no `waitForJobs` to avoid POSIXct compatibility issues in some batchtools versions).

**Flow:** Phase 1 prepare (per-chunk FASTA, manifests, and `loci_seqs.pkl`) → Phase 2 R script (submit and wait) → Phase 3 finish: load `loci_seqs.pkl` and `result.csv` per chunk, merge IntaRNA results into chimeras, write `*.chimeras-r.<n>` and merge/summarize to final outputs.

**Per-chunk files:** Prepare writes `batchtools_work/<n>/query.fa`, `target.fa`, and `loci_seqs.pkl` (locus ID → full sequence). IntaRNA is not asked for seq1/seq2 in the CSV; full sequences come from `loci_seqs.pkl` when merging. Each job reads `query.fa` and `target.fa` and writes `batchtools_work/<n>/result.csv`.

**Important:** The output directory (`-o`) must be on a **shared filesystem** visible from all compute nodes. If `result.csv` is missing or contains only a header, the job likely failed because the path was not visible or writable on the compute node—use an absolute path for `-o` on the shared filesystem and check LSF job logs.

**If the main job times out after submitting IntaRNA jobs:** Run **merge_intarna_into_chimeras.py** after all IntaRNA jobs have completed. It processes all chunks (merge `loci_seqs.pkl` + `result.csv` into `*.chimeras-r.<n>`), then merges `*.chimeras-r.<n>` → `{sample_name}.chimeras.txt`, `*.singletons.<n>` → `{sample_name}.singletons.txt`, and generates `{sample_name}.interactions.txt`. Required args: `--outdir`, `--sample-name`, `--n-chunks`. Optional: `--chunk-root` (default `outdir/batchtools_work`), `--remove-input`, `--buffer-size`. See script docstring and [README.md](README.md) Utility Scripts.

**Options:** Same batchtools options as `chira_map.py` (`--batchtools_queue`, `--batchtools_cores`, `--batchtools_memory`, `--batchtools_walltime`, `--batchtools_template`, `--batchtools_conda_env`, `--batchtools_max_parallel`, `--batchtools_registry`). IntaRNA always runs once per locus pair per chunk (only real chimeric pairs; all-vs-all is not supported).

## R Scripts: Job Waiting and POSIXct Compatibility

Both `submit_chunks_batchtools.R` and `submit_intarna_batchtools.R` use **manual polling** (repeat: `getJobTable` → count completed → `Sys.sleep(POLL_SLEEP_SECONDS)`; IntaRNA script reads sleep from config `poll_sleep_seconds`, default 120) instead of `waitForJobs()`. This avoids errors in environments where batchtools returns `done`/`running`/`error` as POSIXct timestamps (e.g. `'sum' not defined for "POSIXt" objects`). Status counts use a small helper `count_status(x)` that works for both logical and POSIXct columns.

## Files Created

- `process_chunk_batchtools.py`: Standalone script for processing one chunk (called by each batchtools job for mapping)
- `process_intarna_chunk_batchtools.py`: Worker script invoked by IntaRNA batchtools jobs (one run per chunk)
- `submit_chunks_batchtools.R`: R script that submits mapping chunk jobs via batchtools
- `submit_intarna_batchtools.R`: R script that submits IntaRNA jobs via batchtools (used by `chira_extract.py --hybridize --use_batchtools`)
- `merge_intarna_into_chimeras.py`: Standalone script to merge IntaRNA results into chimeras and produce final chimeras.txt, singletons.txt, and interactions.txt (use when the main `chira_extract.py` job times out after submitting IntaRNA jobs)
- `lsf_custom.tmpl`: Default LSF template file (included with ChiRA)
- `output_dir/batchtools_registry_<timestamp>/`: Batchtools registry directory containing:
  - `config.json`: Job configuration with all absolute paths
  - `chunks.json`: List of all chunks with absolute paths (one file for all chunks - intended design)
  - Other batchtools registry files (jobs, logs, etc.)

**Path Handling**: All file paths in configuration files are absolute paths to ensure cluster jobs can correctly resolve files regardless of working directory.

