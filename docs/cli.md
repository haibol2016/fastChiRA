# Command-line reference

Pipeline and utility CLIs for fastChiRA. Install: [install.md](install.md). HPC: [batchtools.md](batchtools.md). Containers: [containers.md](containers.md).

Pipeline order: **collapse → map → merge → quantify → extract**.

**Reference FASTAs (`-f1` / `-f2`):** Use the **same** `-f1, --ref_fasta1` and (when used) `-f2, --ref_fasta2` files at every step that accepts them — typically **map → merge → extract**. Do not swap paths, subsets, or different builds of the same transcriptome between steps; mismatched references break length/context lookups and can corrupt loci, coordinates, and hybridization sequences. If you mapped with only `-f1`, keep using only that file downstream; if you used a split reference (`-f1` targets + `-f2` miRNAs), pass **both** the same way through merge and extract.

---

## Pipeline tools

## chira_collapse.py

Deduplicate FASTQ reads by collapsing identical sequences (and optional UMIs) into a FASTA of unique sequences with counts. Typically the first pipeline step.

### Required

| Option | Description |
|--------|-------------|
| `-i, --fastq` | Input FASTQ (gzipped or uncompressed) |
| `-o, --fasta` | Output FASTA |

### Optional

| Option | Default | Description |
|--------|---------|-------------|
| `-u, --umi_len` | `0` | UMI length; trim from 5′ end and append to the tag id (`0` = no UMI) |
| `-v, --version` | — | Print version and exit |

### Output

FASTA with one entry per unique `(sequence, UMI)`:

- With UMI: `>id|UMI|read_count`
- Without UMI: `>id|read_count`

### Examples

```bash
chira_collapse.py -i input.fastq.gz -o collapsed.fasta
chira_collapse.py -i input.fastq -o collapsed.fasta -u 12
```

---

## chira_map.py

Map collapsed reads to one or two reference transcriptomes with **BWA-MEM** (default, recommended) or **CLAN**. Uses a two-pass strategy (longer then shorter segments) to recover chimeric alignments. Supports strand-specific mapping and optional split references (e.g. targets + mature miRNAs).

**Either** provide a prebuilt index (`-x1`) **or** reference FASTA plus `-b` to build indices. Validation requires `-x1` or `-f1`, and either `-b` or `-x1` ( `-b` and `-x*` are mutually exclusive).

**Chunking** (`--chunk_fasta` / `--use_batchtools`) applies to **BWA only**. CLAN runs without FASTA chunking.

### Required

| Option | Description |
|--------|-------------|
| `-i, --query_fasta` | Query FASTA (usually from `chira_collapse.py`) |
| `-o, --outdir` | Output directory |

Plus index/reference as above: `-x1` / `-x2` **or** `-f1` / `-f2` with `-b`.

### Optional — aligner, strand, processes

| Option | Default | Description |
|--------|---------|-------------|
| `-a, --aligner` | `bwa` | `bwa` or `clan` |
| `-s, --stranded` | `fw` | `fw` (transcript strand; typical for CLASH/CLEAR-CLIP/PARIS/SPLASH), `rc` (antisense), `both` (unstranded) |
| `-p, --processes` | CPU count | Total threads for the main job (BWA jobs share this; with batchtools, LSF jobs use `--batchtools_cores` instead) |
| `--sort_memory` | auto | Memory per thread for BAM sort (e.g. `2G`); total ≈ value × processes |
| `-b, --build` | off | Build indices from `-f1`/`-f2` |
| `-x1, --index1` / `-x2, --index2` | — | Prebuilt index paths |
| `-f1, --ref_fasta1` / `-f2, --ref_fasta2` | — | Reference FASTAs (split reference: e.g. `-f1` non-miRNAs.fa and `-f2` miRNAs.fa; convert miRNA U→T). Reuse these **same** paths in merge and extract |
| `-co, --chimeric_overlap` | `2` | Max bases allowed between chimeric segments |
| `-v, --version` | — | Print version |

### Optional — seed / score (both aligners)

| Option | Default | Description |
|--------|---------|-------------|
| `-l1, --seed_length1` | `12` | 1st-pass seed (`bwa-mem -k`) |
| `-l2, --seed_length2` | `16` | 2nd-pass seed |
| `-s1, --align_score1` | `18` | Min score 1st pass (`bwa-mem -T` / `clan_search -l`) |
| `-s2, --align_score2` | `16` | Min score 2nd pass (must be ≤ 1st) |
| `-h1, --nhits1` | `50` | Max multi-hits per read (1st) |
| `-h2, --nhits2` | `100` | Max multi-hits (2nd) |

#### BWA / CLAN scoring parameters

| Option | Default | Option | Default |
|--------|---------|--------|---------|
| `-ma1/--match1` | `1` | `-ma2/--match2` | `1` |
| `-mm1/--mismatch1` | `4` | `-mm2/--mismatch2` | `6` |
| `-go1/--gapopen1` | `6` | `-go2/--gapopen2` | `100` |
| `-ge1/--gapext1` | `1` | `-ge2/--gapext2` | `100` |

### Optional — chunking (BWA only)

| Option | Default | Description |
|--------|---------|-------------|
| `--chunk_fasta` | off | Split query FASTA into N chunks (useful for large files) |
| `--parallel_chunks` | `2` | Chunks run at once on a single node (ignored when `--use_batchtools`) |

Single-node: processes per chunk ≈ `--processes / --parallel_chunks`. Aim for ≥4 threads per chunk.

### Optional — batchtools (BWA + `--chunk_fasta`)

Requires R with `batchtools` / `jsonlite`. See [batchtools.md](batchtools.md).

| Option | Default | Description |
|--------|---------|-------------|
| `--use_batchtools` | off | Submit one LSF job per chunk |
| `--batchtools_queue` | `long` | LSF queue |
| `--batchtools_cores` | `8` if unset | Cores **per LSF job** (independent of `-p`) |
| `--batchtools_memory` | auto | Total memory per job (e.g. `16GB`); converted to per-core for LSF |
| `--batchtools_walltime` | `240:00` | Walltime per job |
| `--batchtools_max_parallel` | unlimited | Cap concurrent LSF jobs |
| `--batchtools_template` | packaged `lsf_custom.tmpl` | Path to LSF template, or `lsf-simple`. Packaged default: `chira/share/lsf_custom.tmpl` |
| `--batchtools_conda_env` | `CONDA_DEFAULT_ENV` | Conda env for workers |

### Outputs (in `--outdir`)

| File | Description |
|------|-------------|
| `sorted.bam` | Name-sorted BAM of mapped alignments |
| `sorted.bed` | Deduplicated BED of alignments (input to merge) |
| `unmapped.fasta` | Unmapped reads |

### Examples

```bash
# Build indices and map with split reference
chira_map.py -i collapsed.fasta -o map_out -f1 targets.fa -f2 mature_mirna.fa -b -a bwa -s fw -p 16

# Prebuilt indices + large FASTA chunking on one node
chira_map.py -i large.fasta -o map_out -x1 idx1 -x2 idx2 -p 32 --chunk_fasta 10 --parallel_chunks 2

# Cluster: one LSF job per chunk
chira_map.py -i large.fasta -o map_out -f1 targets.fa -f2 mirna.fa -b \
  --chunk_fasta 20 --use_batchtools \
  --batchtools_cores 8 --batchtools_memory 16GB --batchtools_walltime 240:00
```

---

## chira_merge.py

Merge overlapping alignments into read-concentrated loci. With a GTF, convert transcriptomic to genomic coordinates. Writes segmented reads and merged loci for quantification.

### Required

| Option | Description |
|--------|-------------|
| `-b, --bed` | Alignments BED (e.g. `sorted.bed` from map) |
| `-o, --outdir` | Output directory |

### Optional

| Option | Default | Description |
|--------|---------|-------------|
| `-g, --gtf` | — | Annotation GTF/GFF for coordinate conversion |
| `-f1, --ref_fasta1` / `-f2, --ref_fasta2` | — | Reference FASTAs; must match those used in **map** (and later in **extract**) |
| `-ao, --alignment_overlap` | `0.7` | Min fraction overlap to merge BED intervals `[0–1]` |
| `-so, --segment_overlap` | `0.7` | Overlap fraction to merge positions into a segment |
| `-lt, --length_threshold` | `0.9` | Keep alignments ≥ this fraction of the longest `[0.8–1]` |
| `-co, --chimeric_overlap` | `2` | Max bases between chimeric segments |
| `-c, --chimeric_only` | off | Merge chimeric reads only |
| `-ls, --min_locus_size` | `1` | Min alignments per merged locus |
| `-bb, --block_based` | off | Use Blockbuster instead of overlap merging |
| `-d, --distance` | `30` | Blockbuster distance |
| `-mc, --min_cluster_height` | `10` | Blockbuster minClusterHeight |
| `-mb, --min_block_height` | `10` | Blockbuster minBlockHeight |
| `-sc, --scale` | `0.1` | Blockbuster scale |
| `-p, --processes` | auto | Parallel transcript processing (`1` = sequential) |
| `-v, --version` | — | Print version |

### Outputs

| File | Description |
|------|-------------|
| `segments.bed` | Reads split into aligned segments |
| `merged.bed` | Merged loci: `transcript`, `start`, `end`, `strand`, semicolon-joined alignment IDs |

### Examples  

```bash
chira_merge.py -b map_out/sorted.bed -o merge_out -g combined.gtf -f1 targets.fa -p 8
chira_merge.py -b map_out/sorted.bed -o merge_out -g combined.gtf -f1 targets.fa -bb -d 30
```

---

## chira_quantify.py

Build chimeric read loci (CRLs) from merge outputs and quantify multi-mapping reads with an EM algorithm (TPM per CRL).

### Required

| Option | Description |
|--------|-------------|
| `-b, --bed` | Segments BED (`segments.bed` from merge) |
| `-m, --merged_bed` | Merged loci BED (`merged.bed` from merge) |
| `-o, --outdir` | Output directory (writes `loci.counts`) |

### Optional

| Option | Default | Description |
|--------|---------|-------------|
| `-cs, --crl_share` | `0.7` | Min Jaccard-like share of locus reads to join a CRL |
| `-ls, --min_locus_size` | `10` | Min reads for a locus to enter CRL building |
| `-e, --em_threshold` | `1e-5` | EM convergence threshold on transcript expression |
| `-crl, --build_crls_too` | off | Also emit intermediate CRL artifacts during build |
| `-p, --processes` | `0` | EM workers; `0` = all cores; `1` = sequential (MPIRE) |
| `--use_sqldb` | off | SQLite backend for very large inputs (lower RAM, slower) |
| `-v, --version` | — | Print version |

### Output: `loci.counts`

Tab-separated, **no header**, **13 columns**:

| Col | Field | Description |
|-----|-------|-------------|
| 0 | `segment_id` | Read/segment identifier |
| 1 | `transcript_id` | Transcript ID |
| 2 | `locus_id` | Locus identifier |
| 3 | `crl_id` | CRL group ID |
| 4 | `transcript_start` | Transcript start |
| 5 | `transcript_end` | Transcript end |
| 6 | `transcript_strand` | Strand |
| 7 | `cigar` | CIGAR |
| 8 | `genomic_pos` | Genomic position (if GTF used upstream) |
| 9 | `locus_pos` | `chr:start:end:strand` |
| 10 | `locus_share` | Fraction of CRL reads at this locus |
| 11 | `read_crl_fraction` | EM probability that the read belongs to this CRL |
| 12 | `crl_tpm` | TPM for the CRL |

### Examples

```bash
chira_quantify.py -b merge_out/segments.bed -m merge_out/merged.bed -o quant_out -p 8
chira_quantify.py -b merge_out/segments.bed -m merge_out/merged.bed -o quant_out -cs 0.7 -ls 10 --use_sqldb
```

---

## chira_extract.py

Extract best chimeric (and singleton) alignments per read; optionally run IntaRNA hybridization and summarize locus-level interactions. Final reporting step.

### Required

| Option | Description |
|--------|-------------|
| `-l, --loci` | Quantify output (e.g. `loci.counts`) |
| `-o, --out` | Output directory |
| `-n, --sample_name` | Prefix for output files |
| `-f1, --ref_fasta1` | First-priority reference FASTA (same file as in **map** / **merge**) |

**Conditional:** with `--hybridize` **and** `-g/--gtf`, **`-f/--ref`** (genomic FASTA) is required for IntaRNA accessibility.

### Optional — references and filtering

| Option | Default | Description |
|--------|---------|-------------|
| `-g, --gtf` | — | Annotation GTF/GFF (miRBase GFF3 OK for mature miRNAs) |
| `-f2, --ref_fasta2` | — | Second-priority FASTA (e.g. miRNA); same file as in **map** / **merge** when a split reference was used |
| `-f, --ref` | — | Genomic FASTA (IntaRNA accessibility when hybridizing with GTF) |
| `-p, --processes` | `8` | Processes for extraction, hybridization prep/finish, merge |
| `-tc, --tpm_cutoff` | `0` | Drop transcripts below this TPM percentile `[0, 1)` |
| `-sc, --score_cutoff` | `0.0` | Drop hybrids below this score `[0–1)` |
| `-co, --chimeric_overlap` | `2` | Max bases between chimeric segments |
| `-s, --summarize` | off | Write locus-level `interactions.txt` |
| `-z, --gzip` | off | Gzip final chimeras/singletons |
| `--rebuild_index` | off | Rebuild `loci.counts.idx` even if present |
| `--remove_index` | off | Delete `.idx` after run |
| `--remove_intermediate` | off | Delete chunk intermediates after success |
| `-v, --version` | — | Print version |

### Optional — IntaRNA (`-r/--hybridize`)

| Option | Default | Description |
|--------|---------|-------------|
| `-r, --hybridize` | off | Run IntaRNA on predicted chimeras |
| `-ns, --no_seed` | off | Do not enforce seed interactions |
| `-acc, --accessibility` | `N` | `C` = compute accessibility, `N` = not |
| `-m, --intarna_mode` | `H` | `H` heuristic, `M` exact, `S` seed-only |
| `-t, --temperature` | `37` | °C (0–100) |
| `-sbp, --seed_bp` | `5` | IntaRNA `--seedBP` (2–19) |
| `-smpu, --seed_min_pu` | `0` | IntaRNA `--seedMinPu` |
| `-accw, --acc_width` | `150` | IntaRNA `--accW` |

### Optional — batchtools (requires `--hybridize`)

See [batchtools.md](batchtools.md). If the parent job times out after submitting IntaRNA chunks, finish with `merge_intarna_into_chimeras.py` ([cli.md](cli.md#utility-scripts)).

| Option | Default | Description |
|--------|---------|-------------|
| `--use_batchtools` | off | Submit IntaRNA chunk jobs via R batchtools |
| `--batchtools_registry` | `<outdir>/batchtools_work/registry` | Registry directory |
| `--batchtools_template` | packaged `lsf_custom.tmpl` | Template path or `lsf-simple` (`chira/share/lsf_custom.tmpl`) |
| `--batchtools_queue` | `long` | LSF queue |
| `--batchtools_cores` | `1` | Cores per LSF job (IntaRNA is single-threaded; `1` maximizes concurrent jobs) |
| `--batchtools_memory` | — | Total memory per job (converted to per-core) |
| `--batchtools_walltime` | `48:00` | Walltime per job |
| `--batchtools_conda_env` | — | Conda env for workers |
| `--batchtools_max_parallel` | all chunks | Cap concurrent jobs |

### Outputs

Files are prefixed with `{sample_name}`. With `-z`, chimeras/singletons gain `.gz`. Interactions stay uncompressed.

#### `{sample_name}.chimeras.txt` (35 columns; header present)

`read_id`, `transcript_id_1/2`, `gene_id_1/2`, `gene_symbol_1/2`, `annotation_region_1/2`, transcript coords/lengths for both arms, `read_alignment_info` (`arm1_start,arm1_end,arm2_start,arm2_end,read_length`), genomic coords, `locus_id_*`, `crl_group_id_*`, `tpm_*`, alignment scores, then hybridization fields (`hybridized_sequences`, `hybridization_structure`, `hybridization_positions`, `hybridization_mfe_kcal_mol`, `mirna_read_position`, `hybridized_subsequences`; `NA` if not hybridized).

Orientation: compare arm starts in `read_alignment_info` (`arm1_start < arm2_start` → 5′ locus1 → locus2 3′).

#### `{sample_name}.singletons.txt` (15 columns; header)

`read_id`, `transcript_id`, `gene_id`, `gene_symbol`, `annotation_region`, transcript coords/length, `read_alignment_info`, `genomic_coordinates`, `locus_id`, `crl_group_id`, `tpm`, `alignment_score`.

#### `{sample_name}.interactions.txt` (with `-s`; 28 columns; header)

Includes `supporting_read_count`, locus genomic coords, sequences / hybridization fields, TPMs and scores, annotation regions, transcript IDs, and **`gene_name_1`**, **`gene_name_2`**. Comment lines note how to spot miRNA vs target via region annotations (`miRNA`, `3p_mature_mir`, `5p_mature_mir`, `mature_mir`).

### Examples

```bash
# No hybridization
chira_extract.py -l quant_out/loci.counts -o extract_out -n sample1 \
  -f1 targets.fa -g combined.gtf -p 8 -s

# Hybridize + summarize (genomic FASTA required with GTF)
chira_extract.py -l quant_out/loci.counts -o extract_out -n sample1 \
  -f1 targets.fa -f2 mirna.fa -g combined.gtf -f genome.fa \
  -r -s -tc 0.1 -sc 0.1 -p 8 --gzip

# IntaRNA via batchtools
chira_extract.py -l quant_out/loci.counts -o extract_out -n sample1 \
  -f1 targets.fa -g combined.gtf -f genome.fa -r -s -p 32 \
  --use_batchtools --batchtools_cores 1 --batchtools_memory 8GB \
  --batchtools_walltime 48:00
```

---

## Utility scripts

Helper CLIs for reference/annotation prep and IntaRNA recovery after batchtools timeouts.

Typical reference prep for split-reference miRNA–target analysis:

1. `download_ensembl.py` — genome + GTF
2. `download_mirbase_mature.py` / `download_mirbase_gff3.py` — mature miRNAs
3. `remove_mirna_hairpin_from_gtf.py` — strip miRNAs from Ensembl GTF
4. `extract_transcripts_from_genome.py` — target transcriptome FASTA
5. `concatenate_gtf.py` — combine miRNA GFF3 + filtered GTF

**Internal batchtools workers** (not for direct use): `process_map_chunk_batchtools.py` and `process_intarna_chunk_batchtools.py`. Prefer the main pipeline CLIs and [batchtools.md](batchtools.md).

Default LSF template: packaged `chira/share/lsf_custom.tmpl` (via `chira.utilities.share_path`).

## download_mirbase_mature.py

Download species-specific mature miRNA FASTA from miRBase. Sequences are converted **U→T** for DNA-compatible indexing.

### Required

| Option | Description |
|--------|-------------|
| `-s, --species` | Species code (`hsa`, `mmu`, `bta`, `rno`, …) |
| `-o, --output` | Output FASTA |

### Optional

| Option | Default | Description |
|--------|---------|-------------|
| `--mirbase-version` | `CURRENT` | miRBase version (e.g. `22.1`) |
| `--keep-full` | off | Keep temporary full `mature.fa` |
| `--timeout` | `30` | Download timeout (seconds) |
| `-v, --version` | — | Print version |

### Example

```bash
download_mirbase_mature.py -s hsa -o mature_hsa.fasta
```

---

## download_mirbase_gff3.py

Download species-specific miRBase GFF3 (primary transcripts + mature miRNAs). Optional UCSC chain liftover and chromosome renaming. Order: download → liftover → rename. GFF3 can be used directly with fastChiRA (no GTF conversion required).

### Required

| Option | Description |
|--------|-------------|
| `-s, --species` | Species code |
| `-o, --output` | Output GFF3 |

### Optional

| Option | Default | Description |
|--------|---------|-------------|
| `--mirbase-version` | `CURRENT` | miRBase version (e.g. `21`) |
| `--timeout` | `60` | Download timeout (seconds) |
| `-m, --chromosome-mapping` | — | TSV: `gff3_chr<TAB>target_chr` (after liftover) |
| `--source-genome` | — | Source assembly for liftover (required with `--chain-file`) |
| `--target-genome` | — | Target assembly (required with `--chain-file`) |
| `--chain-file` | — | UCSC chain file; needs `pyliftover` |
| `-v, --version` | — | Print version |

### Examples

```bash
download_mirbase_gff3.py -s hsa -o hsa.gff3

download_mirbase_gff3.py -s hsa -o hsa_hg38.gff3 \
  --source-genome hg19 --target-genome hg38 \
  --chain-file hg19ToHg38.over.chain
```

---

## download_ensembl.py

Download Ensembl primary-assembly genome FASTA and GTF for a species/release. Auto-detects assembly unless `-a` is set; HTTP/FTP with fallback; decompresses by default.

### Required

| Option | Description |
|--------|-------------|
| `-s, --species` | e.g. `homo_sapiens`, `mus_musculus` |
| `-g, --genome-version` | Ensembl release for genome (e.g. `110`) |
| `-t, --gtf-version` | Ensembl release for GTF |
| `-o, --output-dir` | Output directory |

### Optional

| Option | Default | Description |
|--------|---------|-------------|
| `-a, --assembly` | auto | Assembly name (e.g. `GRCh38`) |
| `--keep-compressed` | off | Keep `.gz` after decompress |
| `--no-decompress` | off | Leave files gzipped |
| `--timeout` | `60` | Download timeout (seconds) |
| `-v, --version` | — | Print version |

### Example

```bash
download_ensembl.py -s homo_sapiens -g 110 -t 110 -o ./ensembl
```

---

## remove_mirna_hairpin_from_gtf.py

Remove microRNA entries from an Ensembl GTF (by feature type / biotype, optional attribute regex) so target transcriptomes exclude miRNA hairpins.

### Required

| Option | Description |
|--------|-------------|
| `-i, --input` | Input GTF |
| `-o, --output` | Output GTF without miRNA entries |

### Optional

| Option | Default | Description |
|--------|---------|-------------|
| `-p, --pattern` | feature/biotype only | Regex on GTF attributes for extra miRNA matching |
| `--remove-comments` | off | Drop comment lines (default: keep) |
| `-v, --version` | — | Print version |

### Example

```bash
remove_mirna_hairpin_from_gtf.py -i Homo_sapiens.GRCh38.110.gtf -o targets_no_mirna.gtf
```

---

## extract_transcripts_from_genome.py

Extract transcript sequences from a genome FASTA with **gffread**, using a (usually miRNA-filtered) GTF. Prefer this over filtering a pre-built transcriptome FASTA. Requires `gffread` on `PATH` (e.g. `conda install -c bioconda gffread`).

### Required

| Option | Description |
|--------|-------------|
| `-g, --gtf` | Filtered GTF (e.g. from `remove_mirna_hairpin_from_gtf.py`) |
| `-f, --genome-fasta` | Genome FASTA (e.g. Ensembl primary assembly) |
| `-o, --output` | Output transcript FASTA |

### Example

```bash
extract_transcripts_from_genome.py \
  -g targets_no_mirna.gtf \
  -f Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -o target_transcriptome.fasta
```

---

## concatenate_gtf.py

Concatenate mature miRNA GTF/GFF3 with a target transcriptome GTF (comments stripped from the miRNA file). Suitable for split-reference annotation used by merge/extract.

### Required

| Option | Description |
|--------|-------------|
| `-m, --mirna-gtf` | Mature miRNA GTF/GFF3 (e.g. from `download_mirbase_gff3.py`) |
| `-t, --target-gtf` | Target GTF (e.g. from `remove_mirna_hairpin_from_gtf.py`) |
| `-o, --output` | Combined GTF |

### Optional

| Option | Default | Description |
|--------|---------|-------------|
| `--remove-target-comments` | off | Also strip comments from the target GTF |
| `-v, --version` | — | Print version |

### Example

```bash
concatenate_gtf.py -m hsa.gff3 -t targets_no_mirna.gtf -o combined.gtf
```

---

## merge_intarna_into_chimeras.py

Standalone recovery when `chira_extract.py --hybridize --use_batchtools` submits IntaRNA jobs but the parent job exits before the finish phase (e.g. LSF walltime). After all chunk jobs complete:

1. For each chunk `i`: merge `loci_seqs.pkl` + `result.csv` into `{sample}.chimeras-r.{i}`
2. Merge `*.chimeras-r.*` → `{sample}.chimeras.txt`
3. Merge `*.singletons.*` → `{sample}.singletons.txt`
4. Write `{sample}.interactions.txt`

Chunk dirs must contain `loci_seqs.pkl` and `result.csv`; `{sample}.chimeras.{i}` and `{sample}.singletons.{i}` must exist under `--outdir`.

### Required

| Option | Description |
|--------|-------------|
| `--outdir` | Extract output directory |
| `--sample-name` | Sample prefix |
| `--n-chunks` | Number of chunks (`0 .. N-1`) |

### Optional

| Option | Default | Description |
|--------|---------|-------------|
| `--chunk-root` | `<outdir>/batchtools_work` | Parent of chunk directories `0..N-1` |
| `--remove-input` | off | Remove `*.chimeras.<n>` after writing `*.chimeras-r.<n>` |
| `--buffer-size` | package default | I/O buffer size in bytes |

### Example

```bash
merge_intarna_into_chimeras.py \
  --outdir /path/to/extract_out \
  --sample-name sample1 \
  --n-chunks 32
```
