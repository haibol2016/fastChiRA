#!/bin/bash
# Example script for running chira_map.py with chunking and batchtools
# This script demonstrates how to split large FASTA files and process chunks via LSF cluster

# Basic example: Split into 20 chunks, submit to LSF cluster
python /Users/haiboliu/chira/chira_map.py \
  -i input_reads.fasta \
  -o output_directory \
  -x1 /path/to/index1.fa \
  -x2 /path/to/index2.fa \
  --chunk_fasta 20 \
  --use_batchtools \
  --batchtools_queue long \
  --batchtools_cores 8 \
  --batchtools_memory 8GB

# Full example with all options:
# python /Users/haiboliu/chira/chira_map.py \
#   -i input_reads.fasta \
#   -o output_directory \
#   -x1 /path/to/index1.fa \
#   -x2 /path/to/index2.fa \
#   --chunk_fasta 20 \
#   --use_batchtools \
#   --batchtools_queue long \
#   --batchtools_cores 16 \
#   --batchtools_memory 16GB \
#   --batchtools_walltime 240:00 \
#   --batchtools_conda_env chira_env \
#   --batchtools_max_parallel 10

# If you need to build indices first:
# python /Users/haiboliu/chira/chira_map.py \
#   -i input_reads.fasta \
#   -o output_directory \
#   -f1 /path/to/reference1.fa \
#   -f2 /path/to/reference2.fa \
#   --build \
#   --chunk_fasta 20 \
#   --use_batchtools \
#   --batchtools_queue long \
#   --batchtools_cores 8 \
#   --batchtools_memory 8GB

