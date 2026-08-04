# fastChiRA Dockerfile
# micromamba base; Docker + Singularity/Apptainer friendly
# Build: docker build -t chiraplus .
# Run:   docker run --rm -v "$PWD/data:/home/chira/data" -v "$PWD/output:/home/chira/output" chiraplus chira_map.py --help
# Deps:  environment.yml (pinned versions)

FROM mambaorg/micromamba:1.5.8-jammy

USER root

# High UID for Singularity host compatibility
RUN groupadd -r -g 9001 chira && \
    useradd -r -u 9001 -g chira -m -d /home/chira chira && \
    mkdir -p /app /home/chira/{data,output,scratch} && \
    chown -R chira:chira /app /home/chira

WORKDIR /app

# System utils + pinned conda/pip stack from environment.yml
COPY --chown=chira:chira environment.yml ./
RUN apt-get update && \
    apt-get install -y --no-install-recommends curl wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    micromamba install -y -n base -f environment.yml && \
    micromamba clean -afy && \
    # Symlink conda binaries into /usr/local/bin (helps when Singularity overrides PATH)
    for bin in /opt/conda/bin/*; do \
        if [ -f "$bin" ] && [ -x "$bin" ]; then \
            ln -sf "$bin" "/usr/local/bin/$(basename "$bin")" 2>/dev/null || true; \
        fi; \
    done && \
    chown -R chira:chira /opt/conda

ENV CONDA_PREFIX=/opt/conda \
    PATH=/usr/local/bin:/opt/conda/bin:/usr/local/sbin:/usr/sbin:/usr/bin:/sbin:/bin \
    HOME=/home/chira \
    TMPDIR=/tmp \
    PYTHONUNBUFFERED=1

# Install fastChiRA from src/ layout (console scripts → PATH; share/ R+tmpl as package data)
# --no-deps: Python deps already provided by environment.yml
COPY --chown=chira:chira pyproject.toml setup.py MANIFEST.in README.md LICENSE ./
COPY --chown=chira:chira src ./src
RUN /opt/conda/bin/pip install --no-cache-dir --no-deps . && \
    chown -R chira:chira /opt/conda /app && \
    chira_map.py --version && \
    python -c "from BCBio import GFF; from chira import utilities; import os; p=utilities.share_path('lsf_custom.tmpl'); assert os.path.exists(p), p"

# Entrypoint: activate conda, then run command (or bash)
RUN cat > /usr/local/bin/docker-entrypoint.sh << 'EOF'
#!/bin/bash
if [ -f /opt/conda/etc/profile.d/conda.sh ]; then
    # shellcheck source=/dev/null
    source /opt/conda/etc/profile.d/conda.sh
    conda activate base
fi
export PATH="/usr/local/bin:/opt/conda/bin:${PATH}"
if [ $# -eq 0 ]; then
    exec /bin/bash
else
    exec "$@"
fi
EOF

# Make entrypoint executable and symlink for Singularity compatibility
RUN chmod +x /usr/local/bin/docker-entrypoint.sh && \
    ln -sf /usr/local/bin/docker-entrypoint.sh /usr/local/bin/_entrypoint.sh

USER chira

LABEL org.opencontainers.image.title="fastChiRA" \
      org.opencontainers.image.description="fastChiRA - Chimeric Read Analyzer" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
