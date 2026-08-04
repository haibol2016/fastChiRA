# Installation

## Quick install

```bash
conda create -n chira python=3.9 && conda activate chira
conda install -c bioconda bwa samtools bedtools
pip install -e '.[optional]'   # from repository root (pyproject.toml)
chira_map.py --version
```

Optional:

```bash
conda install -c bioconda intarna gffread          # hybridization / transcript extract
conda install -c conda-forge r-base r-batchtools r-jsonlite   # --use_batchtools
```

Package layout is `src/chira/`. After install, CLIs such as `chira_map.py` are on `PATH` (also `python -m chira.map`).

See also: [cli.md](cli.md) · [batchtools.md](batchtools.md) · [containers.md](containers.md)

## Dependencies

| Kind | Packages / tools |
|------|------------------|
| Python (required) | `biopython`, `bcbio-gff`, `pysam`, `mpire>=2.4.0` |
| Python (`pip install '.[optional]'`) | `psutil`, `requests`, `pyliftover` |
| Required CLI | `bwa`, `samtools`, `bedtools` |
| Optional CLI | IntaRNA (`--hybridize`), gffread, blockbuster (`-bb`), CLAN (`-a clan`), GNU sort ≥ 8.6 |
| HPC | R + `batchtools` + `jsonlite`; LSF for the default template |

### Where each is used (short)

| Component | Used by |
|-----------|---------|
| biopython | extract, utilities (FASTA lengths) |
| bcbio-gff | merge, extract (GTF/GFF) |
| pysam | map (BAM merge/sort) |
| mpire | quantify, extract, merge (parallel) |
| psutil | map / utilities (memory & I/O sizing; optional) |
| requests | `download_ensembl.py`, miRBase downloaders |
| pyliftover | `download_mirbase_gff3.py` liftover |
| bwa, samtools | map |
| bedtools | merge, extract |
| IntaRNA | extract `--hybridize` |
| gffread | `extract_transcripts_from_genome.py` |
| R batchtools | map/extract `--use_batchtools` |

## Verify

```bash
chira_map.py --version
chira_merge.py --version
chira_extract.py --version
bwa; samtools --version; bedtools --version
```

## Uninstall

```bash
pip uninstall chira
# optional: conda env remove -n chira
```
