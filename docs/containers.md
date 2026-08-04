# Containers (Docker & Singularity/Apptainer)

Install from source: [install.md](install.md). CLI reference: [cli.md](cli.md).

## Docker

Prebuilt: `docker.io/nemat1976/chiraplus:v0.0.2`  
Or build: `docker build -t chiraplus .` (deps from pinned [`environment.yml`](../environment.yml))

```bash
docker pull docker.io/nemat1976/chiraplus:v0.0.2

docker run --rm \
  -v "$PWD/data:/app/data" \
  -v "$PWD/output:/app/output" \
  docker.io/nemat1976/chiraplus:v0.0.2 \
  chira_collapse.py -i data/input.fastq -o output/collapsed.fasta
```

Interactive shell (no command → bash via entrypoint):

```bash
docker run --rm -it \
  -v "$PWD/data:/app/data" -v "$PWD/output:/app/output" \
  docker.io/nemat1976/chiraplus:v0.0.2
```

Image includes Python deps, `bwa`/`samtools`/`bedtools`/`gffread`/`intarna`, R batchtools, and installed fastChiRA CLIs. Optional tools (blockbuster, CLAN) are not bundled.

## Singularity / Apptainer

```bash
# Pull (or convert a local Docker image)
singularity pull docker://docker.io/nemat1976/chiraplus:v0.0.2

singularity exec \
  -B "$PWD/data:/app/data" \
  -B "$PWD/output:/app/output" \
  chiraplus_v0.0.2.sif \
  chira_collapse.py -i data/input.fastq -o output/collapsed.fasta
```

The image entrypoint (`/usr/local/bin/docker-entrypoint.sh`, also `_entrypoint.sh`) activates conda and sets `PATH`. Call it explicitly only if your site strips entrypoints:

```bash
singularity exec chiraplus_v0.0.2.sif /usr/local/bin/_entrypoint.sh chira_map.py --help
```

### Tips

| Topic | Guidance |
|-------|----------|
| Bind mounts | Map host data/output (and references) into `/app/...` or another fixed path; use the same paths in CLI args |
| Isolation | Use `--cleanenv` if host `PATH`/`PYTHONPATH` interfere |
| Non-root user | Image user `chira` (UID 9001); avoid writing only to root-owned paths |
| Build from Docker | `singularity build chira.sif docker-daemon://chiraplus:latest` after `docker build` |

### Quick checks

```bash
singularity exec chiraplus_v0.0.2.sif chira_map.py --version
singularity exec chiraplus_v0.0.2.sif which bwa samtools bedtools
singularity exec chiraplus_v0.0.2.sif python -c "import chira; print(chira.__version__)"
```
