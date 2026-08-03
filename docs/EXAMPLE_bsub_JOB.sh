#!/bin/bash
# Example LSF job script for running ChiRA with batchtools
# After installing ChiRA with: pip install -e . in conda environment

#BSUB -n 8  # minimal numbers of processors required for a parallel job
#BSUB -R rusage[mem=16000] # ask for memory 16GB
#BSUB -W 240:00 # limit the job to be finished in 240 hours
#BSUB -J "chira_map[5]"  # Array job with 5 tasks
#BSUB -q long  # which queue we want to run in
#BSUB -o logs/out.%J.%I.txt # log
#BSUB -e logs/err.%J.%I.txt # error
#BSUB -R "span[hosts=1]" # All hosts on the same chassis

i=$(($LSB_JOBINDEX - 1))
mkdir -p logs

# Activate conda environment (where ChiRA is installed)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate chira

# Load required modules (if available on your cluster)
module load samtools
module load bwa

set -euo pipefail

# Define input/output paths
non_mir_fa=v47/003.Homo_sapiens.GRCh38.v113.cdna.and.ncrna.exclude.mirna.hairpins_gffread.fa
mir_fa=v47/004.miRBase22.1.human.mature.mir.U2T.fa
gtf=v47/005.human.Gencode.v47.primary_assembly.annotation-exclude.mirna.hairpins.add.mirbase22.mature.mir.gtf
genome_fasta=v47/GRCh38.primary_assembly.genome.fa

in=results/012.collapse.fasta.out/*
fasta=($in/*.fa)
name=(`ls $in/*.fa | perl -p -e s'{.+/(.+?).fa}{$1}'`)
dir=(`ls $in/*.fa | perl -p -e s'{.+/(.+?)/.+?.fa}{$1}'`)
out=results/010.chira.map.out.batchtools/${dir[$i]}/${name[$i]}
mkdir -p $out

index1=v47/index1
index2=v47/index2

# Method 1: Use chira_map.py directly (recommended after pip install -e .)
# The script is now in your PATH after installation
chira_map.py --aligner bwa \
   -i ${fasta[$i]} \
   -o $out \
   --chunk_fasta 10 \
   --use_batchtools \
   --batchtools_queue long \
   --batchtools_cores 8 \
   --batchtools_max_parallel 10 \
   --batchtools_memory 8GB \
   --batchtools_walltime 240:00 \
   --batchtools_conda_env ~/miniconda3/envs/chira \
   --index1 $index1 \
   --index2 $index2 \
   -s fw -p 8 -co 2

# Alternative Method 2: Use python -m (if scripts not in PATH)
# python -m chira.map --aligner bwa \
#    -i ${fasta[$i]} \
#    -o $out \
#    --chunk_fasta 10 \
#    --use_batchtools \
#    --batchtools_queue long \
#    --batchtools_cores 8 \
#    --batchtools_max_parallel 10 \
#    --batchtools_memory 8GB \
#    --batchtools_walltime 240:00 \
#    --batchtools_conda_env ~/miniconda3/envs/chira \
#    --index1 $index1 \
#    --index2 $index2 \
#    -s fw -p 8 -co 2

# Alternative Method 3: Use full path to script (if not installed)
# python /path/to/chira/chira_map.py --aligner bwa \
#    -i ${fasta[$i]} \
#    -o $out \
#    ... (same parameters as above)

