#!/usr/bin/env python
"""
Standalone script to process a single chunk for batchtools job submission.
This script is called by batchtools jobs to process one chunk.
"""
import sys
import os
import json
import chira_utilities

# Import the alignment function from chira_map
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chira_map import align_with_bwa


def process_chunk_standalone(chunk_file, chunk_idx, chunk_dir, 
                             alignment_job_types_json, 
                             per_chunk_processes):
    """
    Process a single chunk - standalone version for batchtools.
    
    Args:
        chunk_file: Path to chunk FASTA file
        chunk_idx: Chunk index number
        chunk_dir: Directory containing chunks
        alignment_job_types_json: JSON string of alignment job types
        per_chunk_processes: Number of CPU processes for BWA
    """
    if not os.path.isfile(chunk_file):
        raise FileNotFoundError(f"Chunk FASTA not found: {chunk_file}")
    if not os.path.isdir(chunk_dir):
        raise NotADirectoryError(f"Chunk directory not found: {chunk_dir}")

    # Parse alignment job types from JSON
    alignment_job_types = json.loads(alignment_job_types_json)

    chunk_outdir = os.path.join(chunk_dir, f"chunk_{chunk_idx:03d}_out")
    os.makedirs(chunk_outdir, exist_ok=True)
    
    chunk_bams = {}
    for job_type in alignment_job_types:
        (align_type, index_type, refindex, seed_length, align_score,
         match_score, mismatch_score, gap_o, gap_e, n_aligns) = job_type
        
        job_name = f"{index_type}.{align_type}"
        try:
            align_with_bwa(align_type, index_type, chunk_file, refindex, chunk_outdir,
                          seed_length, align_score, match_score, mismatch_score,
                          gap_o, gap_e, n_aligns, per_chunk_processes)
            chunk_bam = os.path.join(chunk_outdir, index_type + "." + align_type + ".bam")
            chunk_bams[job_name] = chunk_bam
            print(f"INFO: Completed {job_name} for chunk {chunk_idx}", file=sys.stderr)
        except Exception as e:
            print(f"ERROR: Chunk {chunk_idx} {job_name} failed: {str(e)}", file=sys.stderr)
            raise
    
    # Write completion marker
    completion_file = os.path.join(chunk_outdir, "chunk_complete.txt")
    with open(completion_file, 'w') as f:
        f.write(f"Chunk {chunk_idx} completed successfully\n")
    
    return chunk_bams


def parse_arguments():
    """
    Parse command-line arguments for chunk processing.
    
    Returns:
        tuple: (chunk_file, chunk_idx, chunk_dir, alignment_job_types_json, per_chunk_processes)
    """
    if len(sys.argv) != 6:
        print("Usage: process_chunk_batchtools.py <chunk_file> <chunk_idx> <chunk_dir> <alignment_job_types_json> <per_chunk_processes>", file=sys.stderr)
        sys.exit(1)
    
    chunk_file = sys.argv[1]
    chunk_idx = int(sys.argv[2])
    chunk_dir = sys.argv[3]
    alignment_job_types_json = sys.argv[4]
    per_chunk_processes = int(sys.argv[5])
    
    return chunk_file, chunk_idx, chunk_dir, alignment_job_types_json, per_chunk_processes


def main():
    """Main function to process a single chunk."""
    chunk_file, chunk_idx, chunk_dir, alignment_job_types_json, per_chunk_processes = parse_arguments()
    
    try:
        process_chunk_standalone(chunk_file, chunk_idx, chunk_dir, alignment_job_types_json, per_chunk_processes)
        print(f"SUCCESS: Chunk {chunk_idx} processed", file=sys.stderr)
        sys.exit(0)
    except Exception as e:
        print(f"FAILED: Chunk {chunk_idx} failed: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

