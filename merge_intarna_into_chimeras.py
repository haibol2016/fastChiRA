#!/usr/bin/env python
"""
Standalone helper: process all chunks 0..N-1 (IntaRNA → chimeras-r), then merge and summarize.
This script can be used if the main job times out after submitting batchtools jobs to predict hybridization for each chunk using IntaRNA.

Script flow:
  1. For each chunk i in 0..N-1: merge loci_seqs.pkl + result.csv into outdir/sample_name.chimeras.i
     → write outdir/sample_name.chimeras-r.i (chunk dirs from --chunk-root/i).
  2. Merge *.chimeras-r.<n> → outdir/sample_name.chimeras.txt.
  3. Merge *.singletons.<n> → outdir/sample_name.singletons.txt.
  4. Generate outdir/sample_name.interactions.txt from the merged chimeras.

Usage as script:
    python merge_intarna_into_chimeras.py --outdir <dir> --sample-name <name> --n-chunks N [--chunk-root <dir>] [--remove-input] [--buffer-size <bytes>]

As module: merge_intarna_into_chimeras, merge_chimeras_r_to_txt, merge_singletons_to_txt, write_interactions_txt.
"""
import os
import sys
import argparse
import pickle

# Allow importing from same package when run as script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from chira_extract import (
    parse_intarna_csv,
    _merge_hybrid_into_chimeras,
    merge_files,
    write_interaction_summary,
    BUFFER_SIZE,
)


def merge_intarna_into_chimeras(
    chunk_dir,
    file_chimeras,
    output_file,
    remove_input=False,
    buffer_size=None,
):
    """
    Incorporate IntaRNA result into a chimera chunk to produce the -r output.

    Loads d_loci_seqs from chunk_dir/loci_seqs.pkl and d_hybrids from
    chunk_dir/result.csv, then merges hybrid columns into file_chimeras
    and writes output_file (e.g. *.chimeras-r.<n>).

    Parameters
    ----------
    chunk_dir : str
        Directory containing loci_seqs.pkl and result.csv.
    file_chimeras : str
        Path to the input chimera chunk (e.g. sample.chimeras.<n>).
    output_file : str
        Path to the output file (e.g. sample.chimeras-r.<n>).
    remove_input : bool, optional
        If True, remove file_chimeras after writing output. Default False.
    buffer_size : int, optional
        I/O buffer size. Defaults to chira_extract.BUFFER_SIZE.
    """
    if buffer_size is None:
        buffer_size = BUFFER_SIZE

    loci_pkl = os.path.join(chunk_dir, "loci_seqs.pkl")
    result_csv = os.path.join(chunk_dir, "result.csv")

    d_loci_seqs = {}
    if os.path.exists(loci_pkl):
        with open(loci_pkl, "rb") as f:
            d_loci_seqs = pickle.load(f)
    else:
        sys.stderr.write(f"Warning: {loci_pkl} not found; full sequences may be NA for some rows.\n")

    d_hybrids = parse_intarna_csv(result_csv) if os.path.exists(result_csv) else {}
    _merge_hybrid_into_chimeras(
        file_chimeras,
        output_file,
        d_loci_seqs or None,
        d_hybrids,
        buffer_size,
        remove_input=remove_input,
    )


# Header rows for merged chimera and singleton files (must match chira_extract column order)
HEADER_CHIMERAS = "\t".join([
    "read_id", "transcript_id_1", "transcript_id_2", "gene_id_1", "gene_id_2",
    "gene_symbol_1", "gene_symbol_2", "annotation_region_1", "annotation_region_2",
    "transcript_start_1", "transcript_end_1", "transcript_strand_1", "transcript_length_1",
    "transcript_start_2", "transcript_end_2", "transcript_strand_2", "transcript_length_2",
    "read_alignment_info", "genomic_coordinates_1", "genomic_coordinates_2",
    "locus_id_1", "locus_id_2", "crl_group_id_1", "crl_group_id_2",
    "tpm_1", "tpm_2", "alignment_score_1", "alignment_score_2", "combined_alignment_score",
    "hybridized_sequences", "hybridization_structure", "hybridization_positions",
    "hybridization_mfe_kcal_mol", "mirna_read_position", "hybridized_subsequences",
])
HEADER_SINGLETONS = "\t".join([
    "read_id", "transcript_id", "gene_id", "gene_symbol", "annotation_region",
    "transcript_start", "transcript_end", "transcript_strand", "transcript_length",
    "read_alignment_info", "genomic_coordinates", "locus_id", "crl_group_id",
    "tpm", "alignment_score",
])


def merge_chimeras_r_to_txt(outdir, sample_name, n_processes, compress=False, remove_intermediate=False):
    """
    Merge *.chimeras-r.<n> files into a single sample_name.chimeras.txt (or .gz if compress).

    Looks for sample_name.chimeras-r.0, .1, ... .(n_processes-1) under outdir and merges
    them with a header into outdir/sample_name.chimeras.txt.
    """
    chimeras_prefix = os.path.join(outdir, sample_name + ".chimeras-r")
    chimeras_file = os.path.join(outdir, sample_name + ".chimeras.txt")
    if compress:
        chimeras_file += ".gz"
    merge_files(chimeras_prefix, chimeras_file, HEADER_CHIMERAS, n_processes, compress, remove_intermediate)


def merge_singletons_to_txt(outdir, sample_name, n_processes, compress=False, remove_intermediate=False):
    """
    Merge *.singletons.<n> files into a single sample_name.singletons.txt (or .gz if compress).

    Looks for sample_name.singletons.0, .1, ... .(n_processes-1) under outdir and merges
    them with a header into outdir/sample_name.singletons.txt.
    """
    singletons_prefix = os.path.join(outdir, sample_name + ".singletons")
    singletons_file = os.path.join(outdir, sample_name + ".singletons.txt")
    if compress:
        singletons_file += ".gz"
    merge_files(singletons_prefix, singletons_file, HEADER_SINGLETONS, n_processes, compress, remove_intermediate)


def write_interactions_txt(outdir, sample_name, compress=False, num_threads=4):
    """
    Generate sample_name.interactions.txt from the merged sample_name.chimeras.txt.

    Reads outdir/sample_name.chimeras.txt (or .gz if compress), aggregates by locus pair,
    and writes outdir/sample_name.interactions.txt with supporting read counts and metadata.
    """
    write_interaction_summary(outdir, sample_name, compress=compress, num_threads=num_threads)


def main():
    parser = argparse.ArgumentParser(
        description="Process all chunks 0..N-1 (IntaRNA → chimeras-r), then merge and summarize to chimeras.txt, singletons.txt, interactions.txt."
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--sample-name", required=True, help="Sample name prefix for output files")
    parser.add_argument(
        "--n-chunks",
        type=int,
        required=True,
        metavar="N",
        help="Number of chunks (0 to N-1). For each i: merge outdir/sample_name.chimeras.i -> chimeras-r.i using chunk_root/i.",
    )
    parser.add_argument(
        "--chunk-root",
        help="Parent of chunk dirs (default: outdir/batchtools_work). Each chunk_root/i must contain loci_seqs.pkl and result.csv",
    )
    parser.add_argument(
        "--remove-input",
        action="store_true",
        help="Remove each *.chimeras.<n> after writing *.chimeras-r.<n>",
    )
    parser.add_argument(
        "--buffer-size",
        type=int,
        default=None,
        help=f"I/O buffer size in bytes (default: {BUFFER_SIZE})",
    )
    args = parser.parse_args()

    outdir = os.path.abspath(args.outdir)
    sample_name = args.sample_name
    chunk_root = os.path.abspath(args.chunk_root or os.path.join(outdir, "batchtools_work"))
    buffer_size = args.buffer_size or BUFFER_SIZE

    for i in range(args.n_chunks):
        chunk_dir = os.path.join(chunk_root, str(i))
        file_chimeras = os.path.join(outdir, f"{sample_name}.chimeras.{i}")
        output_file = os.path.join(outdir, f"{sample_name}.chimeras-r.{i}")
        if not os.path.exists(file_chimeras):
            sys.stderr.write(f"Skip chunk {i}: {file_chimeras} not found\n")
            continue
        if not os.path.isdir(chunk_dir):
            sys.stderr.write(f"Skip chunk {i}: chunk dir not found: {chunk_dir}\n")
            continue
        merge_intarna_into_chimeras(
            chunk_dir=chunk_dir,
            file_chimeras=file_chimeras,
            output_file=output_file,
            remove_input=args.remove_input,
            buffer_size=buffer_size,
        )

    merge_chimeras_r_to_txt(outdir, sample_name, args.n_chunks, False, False)
    merge_singletons_to_txt(outdir, sample_name, args.n_chunks, False, False)
    write_interactions_txt(outdir, sample_name, False, 4)


if __name__ == "__main__":
    main()
