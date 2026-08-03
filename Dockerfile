# ChiRA Dockerfile
# Uses micromamba for lightweight conda package management
# Designed for compatibility with both Docker and Singularity/Apptainer

FROM mambaorg/micromamba:1.5.8-jammy

# Create non-root user for Singularity compatibility
# Use high-range UID to avoid conflicts with host users
USER root
RUN groupadd -r -g 9001 chira && \
    useradd -r -u 9001 -g chira -m -d /home/chira chira && \
    mkdir -p /app /home/chira/{data,output,scratch} && \
    chown -R chira:chira /app /home/chira

WORKDIR /app

# Install system dependencies and set up conda in single layer
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        curl \
        wget \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    # Install conda packages
    micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        biopython \
        bcbiogff \
        pysam \
        requests \
        pyliftover \
        psutil \
        bwa \
        samtools \
        bedtools \
        gffread \
        intarna \
        r-base \
        r-batchtools \
        r-jsonlite \
        r-future \
        r-future.apply \
        && \
    # Install MPIRE with version requirement (>=2.4.0 for shared_objects support)
    # Try conda-forge first, fallback to pip if version is too old
    (micromamba install -y -n base -c conda-forge "mpire>=2.4.0" 2>/dev/null || \
     /opt/conda/bin/pip install --no-cache-dir "mpire>=2.4.0") && \
    micromamba clean -afy && \
    # Create symlinks for ALL conda binaries in /usr/local/bin so they're always in PATH
    # This ensures tools are available even if PATH is overridden in Singularity
    # This includes Python, which needs to be accessible for the scripts to work
    # Using symlinks (not copies) to save space and avoid duplication
    for bin in /opt/conda/bin/*; do \
        if [ -f "$bin" ] && [ -x "$bin" ]; then \
            ln -sf "$bin" /usr/local/bin/$(basename "$bin") 2>/dev/null || true; \
        fi; \
    done && \
    # Fix permissions for conda directory
    chown -R chira:chira /opt/conda

# Set environment variables for both Docker and Singularity
# Note: ENV variables are preserved when converting Docker to Singularity
ENV CONDA_PREFIX=/opt/conda
# Include /opt/conda/bin in PATH so Python and Python packages are accessible
# /usr/local/bin has symlinks to tools, but Python needs to find its packages
ENV PATH=/usr/local/bin:/opt/conda/bin:/usr/local/sbin:/usr/sbin:/usr/bin:/sbin:/bin
ENV HOME=/home/chira
ENV TMPDIR=/tmp
# Unbuffered Python stdout/stderr so logs appear in real time (e.g. docker logs, batch jobs)
ENV PYTHONUNBUFFERED=1

# Copy packaging metadata and source tree, then install the package
COPY --chown=chira:chira pyproject.toml setup.py MANIFEST.in README.md LICENSE ./
COPY --chown=chira:chira src ./src
RUN /opt/conda/bin/pip install --no-cache-dir '.[optional]' && \
    chown -R chira:chira /opt/conda /app

# Generate entrypoint script inline (makes Dockerfile self-contained)
RUN cat > /usr/local/bin/docker-entrypoint.sh << 'EOF'
#!/bin/bash
# Docker entrypoint script for ChiRA
# Sets up conda environment and ensures PATH is correct
# This script is automatically executed in Docker and can be called manually in Singularity

# Activate conda base environment
# This ensures conda's bin directory is in PATH
if [ -f /opt/conda/etc/profile.d/conda.sh ]; then
    source /opt/conda/etc/profile.d/conda.sh
    conda activate base
fi

# Ensure /usr/local/bin is in PATH (where we symlinked all binaries)
export PATH="/usr/local/bin:/opt/conda/bin:${PATH}"

# If arguments are provided, execute them; otherwise, start bash
if [ $# -eq 0 ]; then
    exec /bin/bash
else
    exec "$@"
fi
EOF

RUN chmod +x /usr/local/bin/docker-entrypoint.sh && \
    # Also create a symlink with a shorter name for easier Singularity usage
    ln -sf /usr/local/bin/docker-entrypoint.sh /usr/local/bin/_entrypoint.sh && \
    chown -R chira:chira /app

# Switch to non-root user
USER chira

# Set entrypoint - this will be executed in Docker automatically
# In Singularity, users can call it manually: singularity exec image.sif /usr/local/bin/_entrypoint.sh command
# The entrypoint script handles default command (bash) when no arguments are provided
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
