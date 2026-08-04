# fastChiRA — Chimeric Read Analyzer

**Version:** 1.4.15 · **License:** GPL-3.0-or-later

fastChiRA analyzes RNA–RNA interactome data (CLASH, CLEAR-CLIP, PARIS, SPLASH, and similar). It is a performance-oriented fork of [ChiRA v1.4.3](https://github.com/pavanvidem/chira) with parallel processing, FASTA chunking, and HPC batchtools support. Full history: [CHANGELOG.md](CHANGELOG.md).

Pipeline:

```
FASTQ → collapse → map → merge → quantify → extract → chimeras / singletons [/ interactions]
```

---

## Improvements vs ChiRA ([pavanvidem/chira](https://github.com/pavanvidem/chira/))

Relative to [ChiRA v1.4.3](https://github.com/pavanvidem/chira). Upstream ships only the five core scripts; fastChiRA keeps those CLIs and adds parallelism, HPC support, packaging, and reference-prep utilities.

### Parallel computing (MPIRE)

- **Required** `mpire>=2.4.0` for multi-process work in `chira_merge.py`, `chira_quantify.py`, and `chira_extract.py`
- `WorkerPool` with **shared objects** (large RAM savings and faster startup vs pickling big dicts to every worker)
- `-p/--processes` on merge / quantify / extract (quantify: `0` = all cores)

### HPC batchtools (LSF)

Not in upstream ChiRA. Packaged under `chira/share/` (`submit_*.R`, `lsf_custom.tmpl`):

- **Map:** `--chunk_fasta` + `--use_batchtools` → one LSF job per FASTA chunk (`process_map_chunk_batchtools.py`)
- **Extract IntaRNA:** `-r` + `--use_batchtools` → one IntaRNA job per chunk (`process_intarna_chunk_batchtools.py`)
- Recovery if the parent job times out: `merge_intarna_into_chimeras.py`
- Absolute paths, conda env passthrough, queue/cores/memory/walltime — [docs/batchtools.md](docs/batchtools.md)

### Utility scripts (new)

| Script | Purpose |
|--------|---------|
| `download_ensembl.py` | Ensembl GTF + genome FASTA |
| `download_mirbase_mature.py` | Mature miRNA FASTA (U→T) |
| `download_mirbase_gff3.py` | miRBase GFF3 (+ optional liftover) |
| `remove_mirna_hairpin_from_gtf.py` | Strip miRNAs from Ensembl GTF |
| `extract_transcripts_from_genome.py` | Transcript FASTA via gffread |
| `concatenate_gtf.py` | Combine miRNA GFF/GTF + target GTF |

### Shared utilities & I/O

- **File I/O buffers** (larger than Python defaults to cut syscall overhead on big FASTA/BED/tables):
  - **Adaptive** via `chira.utilities.get_adaptive_buffer_size()` (≈8–16 MB when `psutil` is installed, else a safe default; scales down with concurrent open files):
    - **`chira_map.py`** — splitting query FASTA into chunks; writing `mapped.bed` / `unmapped.fasta`; reading CLAN `out.map`
    - **`chira_merge.py`** — reading input BED / writing `segments.temp.bed`; loading loci for overlap or blockbuster merge; writing `merged.bed`; GTF→exon BED parsing; genomic coordinate conversion
  - **Fixed 2 MB** (`BUFFER_SIZE`) in **`chira_collapse.py`** (FASTQ→FASTA) and **`chira_extract.py`** / **`merge_intarna_into_chimeras.py`** (GTF, chimera/singleton chunks, IntaRNA pair FASTAs, merge, interactions)
  - **Fixed 1 MB** in **`chira_quantify.py`** for CRL / `loci.counts` text I/O (large inputs can also use mmap)
- BEDTools old/new command detection
- GNU sort `--parallel` when available
- Unbuffered logging suited to batch jobs

### Per-tool (vs upstream scripts)

#### [`chira_collapse.py`](https://github.com/pavanvidem/chira/blob/master/chira_collapse.py)

- Raw FASTQ parsing (no Biopython `SeqIO` on the hot path); gzipped FASTQ (`.gz`)

#### [`chira_map.py`](https://github.com/pavanvidem/chira/blob/master/chira_map.py)

- Parallel BWA (`-p`); FASTA chunking (`--chunk_fasta`, `--parallel_chunks`)
- **Batchtools (LSF)** when `--chunk_fasta N --use_batchtools`: one job per FASTA chunk via `submit_chunks_batchtools.R` → `process_map_chunk_batchtools.py` (not in upstream ChiRA)
- Adaptive BAM sort memory; stranded-filter fixes; cluster-safe absolute paths

#### [`chira_merge.py`](https://github.com/pavanvidem/chira/blob/master/chira_merge.py)

- Parallel locus merging (MPIRE; upstream is single-process)
- Chunked transcript workload; cached ref-id / regex / batched writes / adaptive buffers
- Empty-file handling; parallel ≡ sequential results

#### [`chira_quantify.py`](https://github.com/pavanvidem/chira/blob/master/chira_quantify.py)

- Parallel EM via MPIRE (`-p`; default all cores)
- Optional SQLite backend for huge inputs (`--use_sqldb`)
- mmap-friendly I/O; EM/CRL correctness fixes

#### [`chira_extract.py`](https://github.com/pavanvidem/chira/blob/master/chira_extract.py)

- Parallel extraction / merge (MPIRE); `-n/--sample_name`, `-z/--gzip`
- **Batchtools (LSF)** for IntaRNA when `-r/--hybridize --use_batchtools`: one job per chunk via `submit_intarna_batchtools.R` → `process_intarna_chunk_batchtools.py` (not in upstream ChiRA)
- If the parent extract job times out after submitting those jobs: `merge_intarna_into_chimeras.py` finishes chimeras / singletons / interactions
- Interactions: `gene_name_1` / `gene_name_2`; single `hybridization_genomic_coordinates` column

### Packaging & containers

- Installable package (`src/chira/`, `pip install -e .`) with console entry points
- Docker / Singularity path — [docs/containers.md](docs/containers.md)

---

## Installation

Package layout is `src/chira/`. CLI names (`chira_map.py`, etc.) are unchanged after install.

```bash
conda create -n chira python=3.9 && conda activate chira
conda install -c bioconda bwa samtools bedtools
# optional: intarna gffread; for HPC: conda install -c conda-forge r-base r-batchtools r-jsonlite
pip install -e '.[optional]'          # from repo root
chira_map.py --version
```

### Dependencies

| Kind | Packages / tools |
|------|------------------|
| Python (required) | `biopython`, `bcbio-gff`, `pysam`, `mpire>=2.4.0` |
| Python (`.[optional]`) | `psutil`, `requests`, `pyliftover` |
| Required CLI | `bwa`, `samtools`, `bedtools` |
| Optional CLI | IntaRNA (`--hybridize`), gffread (transcript extract), blockbuster (`-bb`), CLAN (`-a clan`), GNU sort ≥ 8.6 |
| HPC | R + `batchtools` + `jsonlite`; LSF (`bsub`) for default template |

### Container

Prebuilt image: `docker.io/nemat1976/chiraplus:v0.0.2` (or `docker build -t chiraplus .`).

```bash
docker run --rm -v "$PWD/data:/app/data" -v "$PWD/output:/app/output" \
  docker.io/nemat1976/chiraplus:v0.0.2 \
  chira_collapse.py -i data/input.fastq -o output/collapsed.fasta
```

Details: [docs/install.md](docs/install.md) · [docs/cli.md](docs/cli.md) · [docs/batchtools.md](docs/batchtools.md) · [docs/containers.md](docs/containers.md)

---

## Quick workflow

Typical split-reference analysis (miRNA + target transcriptome):

```bash
# Reference preparation
download_mirbase_mature.py -s hsa -o mature_mirna.fa
download_ensembl.py -s homo_sapiens -g 115 -t 115 -o ensembl_files
download_mirbase_gff3.py -s hsa -o hsa.gff3
remove_mirna_hairpin_from_gtf.py \
  -i ensembl_files/Homo_sapiens.GRCh38.115.gtf -o target.gtf
extract_transcripts_from_genome.py -g target.gtf \
  -f ensembl_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa -o target.fa
concatenate_gtf.py -m hsa.gff3 -t target.gtf -o annotation.gtf

# Pipeline
chira_collapse.py -i raw.fastq.gz -o reads.fa -u 6

chira_map.py -i reads.fa -o map_out \
  -f1 mature_mirna.fa -f2 target.fa -b -a bwa -s fw -p 8
# or with prebuilt indices: -x1 idx1 -x2 idx2
# large FASTA / HPC: --chunk_fasta 10 [--use_batchtools ...]

chira_merge.py -b map_out/sorted.bed -o merge_out -g annotation.gtf \
  -f1 mature_mirna.fa -f2 target.fa -ao 0.7 -so 0.7 -p 8

chira_quantify.py -b merge_out/segments.bed -m merge_out/merged.bed \
  -o quant_out -cs 0.7 -ls 10 -p 0

chira_extract.py -l quant_out/loci.counts -o extract_out -n sample1 \
  -f1 mature_mirna.fa -f2 target.fa -g annotation.gtf -s -p 8
# with IntaRNA: add -r -f ensembl_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa
# IntaRNA on LSF: add -r --use_batchtools [--batchtools_cores 1 ...]
```

---

## Tool documentation

Full option lists, defaults, and outputs: [docs/cli.md](docs/cli.md)

Also: `TOOL.py --help` after install.

### `chira_collapse.py`

Collapse FASTQ (plain or gzip) to unique FASTA; optional 5′ UMI trim.

| | |
|--|--|
| Required | `-i` FASTQ, `-o` FASTA |
| Optional | `-u/--umi_len` (default `0`) |
| Output | FASTA headers `{id}\|count` or `{id}\|umi\|count` |

### `chira_map.py`

Map reads with BWA (default) or CLAN; write sorted BAM/BED.

| | |
|--|--|
| Required | `-i` query FASTA, `-o` outdir; and either `-x1` **or** (`-f1` + `-b`) |
| Refs | `-x2` / `-f2` for second-priority index/FASTA |
| Key options | `-a {bwa,clan}` (default `bwa`); `-s {fw,rc,both}` (default `fw`); `-p` processes (default: all CPUs); `-co` chimeric overlap (default `2`) |
| Chunking (BWA only) | `--chunk_fasta N`; `--parallel_chunks` (default `2`, single-node) |
| Batchtools | `--use_batchtools` (requires `--chunk_fasta`); queue `long`; cores **8 if unset**; walltime `240:00`; template packaged `lsf_custom.tmpl` |
| Outputs | `sorted.bam`, `sorted.bed`, `unmapped.fasta` |

### `chira_merge.py`

Merge overlapping alignments into segments/loci; optional GTF genomic coordinates.

| | |
|--|--|
| Required | `-b` BED (e.g. `sorted.bed`), `-o` outdir |
| Optional | `-g` GTF; `-f1`/`-f2`; `-ao`/`-so` 0.7; `-lt` 0.9; `-co` 2; `-ls` 1; `-c` chimeric-only; `-bb` blockbuster; `-p` processes (default auto) |
| Outputs | `segments.bed`, `merged.bed` |

### `chira_quantify.py`

Build chimeric read loci (CRLs) and EM-quantify multimappers → TPM / read–CRL fractions.

| | |
|--|--|
| Required | `-b` segments.bed, `-m` merged.bed, `-o` outdir |
| Optional | `-cs` 0.7; `-ls` 10; `-e` 1e-5; `-crl`; `-p` (default `0` = all cores); `--use_sqldb` for huge inputs |
| Output | `loci.counts` (**no header**; 13 columns — see [cli.md](docs/cli.md)) |

### `chira_extract.py`

Extract chimeric and singleton reads from `loci.counts`; optional IntaRNA and locus summary.

| | |
|--|--|
| Required | `-l` loci.counts, `-o` outdir, `-n` sample name, `-f1` ref FASTA |
| Optional refs | `-g` GTF; `-f2` second FASTA; `-f` genomic FASTA (**required if** `-r` **and** `-g`) |
| Key options | `-p` (default **8**); `-tc`/`-sc` cutoffs; `-co` 2; `-s` summarize; `-z` gzip chimeras/singletons; `-r` hybridize |
| IntaRNA | `-ns`, `-acc {C,N}` (default `N`), `-m {H,M,S}` (default `H`), `-t` 37, `-sbp` 5, … |
| Batchtools | `--use_batchtools` (with `-r`); cores default **1**; walltime `48:00`; queue `long`; `--remove_intermediate` |
| Outputs | `{sample}.chimeras.txt`[.gz], `{sample}.singletons.txt`[.gz], and with `-s`: `{sample}.interactions.txt` (column lists in [cli.md](docs/cli.md)) |

---

## Utility scripts

| Script | Role |
|--------|------|
| `download_mirbase_mature.py` | Mature miRNA FASTA by species (`-s`, `-o`); U→T |
| `download_mirbase_gff3.py` | miRBase GFF3; optional liftover (`pyliftover`) |
| `download_ensembl.py` | Ensembl GTF + genome FASTA (`requests`) |
| `remove_mirna_hairpin_from_gtf.py` | Strip miRNA entries from Ensembl GTF |
| `extract_transcripts_from_genome.py` | Transcript FASTA via **gffread** (`-g`, `-f`, `-o`) |
| `concatenate_gtf.py` | Concatenate miRNA GFF/GTF + target GTF |
| `merge_intarna_into_chimeras.py` | Recover finals if extract IntaRNA batchtools job times out |

Full options: [docs/cli.md](docs/cli.md#utility-scripts).

---

## HPC batchtools

Requires R with `batchtools` + `jsonlite`. Packaged under `chira/share/`: `submit_chunks_batchtools.R`, `submit_intarna_batchtools.R`, `lsf_custom.tmpl`.

**Map:** `--chunk_fasta N --use_batchtools` → one LSF job per chunk via `process_map_chunk_batchtools.py`. Cores **8 if unset**; walltime `240:00`.

**Extract IntaRNA:** `-r --use_batchtools` → one IntaRNA pairwise run per chunk via `process_intarna_chunk_batchtools.py`. Cores default **1**; walltime `48:00`; work dir `<outdir>/batchtools_work/`.

If the extract job times out after IntaRNA jobs finish:

```bash
merge_intarna_into_chimeras.py --outdir extract_out --sample-name sample1 --n-chunks N
```

Guide: [docs/batchtools.md](docs/batchtools.md)

---

## Package layout

```
src/chira/
  collapse.py  map.py  merge.py  quantify.py  extract.py  utilities.py
  batchtools/   # cluster workers + merge_intarna helper
  utils/        # reference-prep CLIs
  share/        # *.R, lsf_custom.tmpl
docs/           # install, cli, batchtools, containers
```

Importable as `import chira`; version: `chira.__version__`.
