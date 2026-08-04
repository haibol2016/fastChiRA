# Batchtools (LSF)

Submit chunk jobs to an LSF cluster via R **batchtools** + **jsonlite**. Packaged assets (under `chira/share/`):

- `submit_chunks_batchtools.R` — map chunks → `process_map_chunk_batchtools.py`
- `submit_intarna_batchtools.R` — IntaRNA chunks → `process_intarna_chunk_batchtools.py`
- `lsf_custom.tmpl` — default LSF template (also used when cwd is the share dir)

Install: [install.md](install.md). CLI flags: [cli.md](cli.md).

## Prerequisites

```bash
conda install -c conda-forge r-base r-batchtools r-jsonlite
```

Also: LSF (`bsub`), shared filesystem, absolute paths (fastChiRA converts paths automatically).

## Map (`chira_map.py`)

Requires `--chunk_fasta N` and `--use_batchtools` (BWA only).

```bash
chira_map.py -i reads.fa -o map_out -x1 idx1 -x2 idx2 \
  --chunk_fasta 20 --use_batchtools \
  --batchtools_queue long \
  --batchtools_cores 8 \
  --batchtools_memory 16GB \
  --batchtools_walltime 240:00 \
  --batchtools_conda_env ~/miniconda3/envs/chira \
  -s fw -p 8
```

| Option | Default | Notes |
|--------|---------|-------|
| `--batchtools_queue` | `long` | LSF queue |
| `--batchtools_cores` | **8 if unset** | Cores **per LSF job** (independent of `-p`) |
| `--batchtools_memory` | auto | Total memory per job; converted to per-core for LSF |
| `--batchtools_walltime` | `240:00` | Per job |
| `--batchtools_max_parallel` | unlimited | Cap concurrent jobs |
| `--batchtools_template` | packaged `lsf_custom.tmpl` | Or path / `lsf-simple` |
| `--batchtools_conda_env` | `$CONDA_DEFAULT_ENV` | Absolute path recommended |

Monitor: `bjobs -u $USER`, `bjobs -J 'chira_bt_*'`.

## Extract IntaRNA (`chira_extract.py`)

Requires `-r/--hybridize` and `--use_batchtools`. One IntaRNA pairwise run per chunk.

```bash
chira_extract.py -l quant_out/loci.counts -o extract_out -n sample1 \
  -f1 mirna.fa -f2 target.fa -g annotation.gtf -f genome.fa \
  -r -s -p 8 \
  --use_batchtools \
  --batchtools_cores 1 \
  --batchtools_memory 8GB \
  --batchtools_walltime 48:00
```

| Option | Default | Notes |
|--------|---------|-------|
| `--batchtools_cores` | **1** | IntaRNA is single-threaded; `1` maximizes concurrency |
| `--batchtools_walltime` | `48:00` | Per job |
| `--batchtools_registry` | `<outdir>/batchtools_work/registry` | Registry dir |
| `--remove_intermediate` | off | Clean chunk intermediates after success |

Work dirs: `<outdir>/batchtools_work/<n>/{query.fa,target.fa,result.csv,loci_seqs.pkl}`.

### Recovery if the parent job times out

After all IntaRNA LSF jobs finish:

```bash
merge_intarna_into_chimeras.py \
  --outdir extract_out --sample-name sample1 --n-chunks N
```

## LSF job wrapper example

Submit a head job that itself calls `chira_map.py --use_batchtools` (array index optional):

```bash
#!/bin/bash
#BSUB -n 8
#BSUB -R "rusage[mem=16000] span[hosts=1]"
#BSUB -W 240:00
#BSUB -J chira_map
#BSUB -q long
#BSUB -o logs/out.%J.txt
#BSUB -e logs/err.%J.txt

source ~/miniconda3/etc/profile.d/conda.sh
conda activate chira
mkdir -p logs

chira_map.py --aligner bwa \
  -i reads.fa -o map_out \
  --index1 /path/to/index1 --index2 /path/to/index2 \
  --chunk_fasta 10 --use_batchtools \
  --batchtools_queue long \
  --batchtools_cores 8 \
  --batchtools_memory 8GB \
  --batchtools_walltime 240:00 \
  --batchtools_conda_env ~/miniconda3/envs/chira \
  -s fw -p 8 -co 2
```

## LSF templates

**Packaged default:** `chira/share/lsf_custom.tmpl` (conda activate, modules, resource lines).

**Built-in `lsf-simple`** (batchtools package) roughly:

```bash
#!/bin/bash
#BSUB <%= resources %>
#BSUB -J <%= job.name %>
#BSUB -o <%= log.file %>
#BSUB -e <%= log.file %>
<%= rscript %>
```

Prefer the packaged `lsf_custom.tmpl` unless you need a custom site template (`--batchtools_template /path/to/file.tmpl`).

## Notes

- Both R submit scripts poll with `getJobTable()` + `Sys.sleep` (POSIXct-safe `count_status()`); they do not rely on `waitForJobs()`.
- Leave `--batchtools_max_parallel` unset to submit all chunks at once when the queue allows.
