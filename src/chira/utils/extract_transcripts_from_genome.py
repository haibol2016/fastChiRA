#!/usr/bin/env python
"""
Extract transcript FASTA sequences from a genome FASTA file using gffread.

This script uses gffread (https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread) to extract transcript sequences from a genome
FASTA file based on a filtered GTF file (e.g., output from remove_mirna_hairpin_from_gtf.py).

gffread extracts sequences based on transcript features in the GTF file, producing
a transcriptome FASTA file without miRNA sequences (if the GTF was filtered to remove miRNAs).
"""

import argparse
import sys
import subprocess
import os

# Real-time output when stdout is not a TTY (e.g. batch jobs)
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)


def check_gffread():
    """
    Check if gffread is available in the system PATH.
    
    Returns:
        True if gffread is available, False otherwise
    """
    try:
        result = subprocess.run(['gffread', '-h'], 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              timeout=5)
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def extract_transcripts_with_gffread(gtf_file, genome_fasta, output_fasta):
    """
    Extract transcript sequences from genome FASTA using gffread.
    
    Args:
        gtf_file: Path to filtered GTF file (e.g., from remove_mirna_hairpin_from_gtf.py)
        genome_fasta: Path to genome FASTA file
        output_fasta: Path to output transcript FASTA file
    
    Returns:
        True if successful, False otherwise
    """
    # Check if gffread is available
    if not check_gffread():
        print("Error: gffread not found in PATH", file=sys.stderr)
        print("Please install gffread (GFF utilities from Johns Hopkins University):", file=sys.stderr)
        print("  conda install -c bioconda gffread", file=sys.stderr)
        print("  or", file=sys.stderr)
        print("  Download from: https://github.com/gpertea/gffread", file=sys.stderr)
        print("  Documentation: https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread", file=sys.stderr)
        return False
    
    # Check if input files exist
    if not os.path.exists(gtf_file):
        print(f"Error: GTF file '{gtf_file}' not found", file=sys.stderr)
        return False
    
    if not os.path.exists(genome_fasta):
        print(f"Error: Genome FASTA file '{genome_fasta}' not found", file=sys.stderr)
        return False
    
    # Build gffread command
    # -w: write transcript sequences to FASTA file
    # -g: genome FASTA file
    # -W: use transcript features (default)
    # Last argument: input GTF file
    cmd = ['gffread', '-w', output_fasta, '-g', genome_fasta, gtf_file]
    
    try:
        print(f"Running: {' '.join(cmd)}")
        print(f"Extracting transcripts from genome FASTA using GTF file...")
        
        result = subprocess.run(cmd, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              text=True,
                              check=False)
        
        if result.returncode != 0:
            print(f"Error: gffread failed with return code {result.returncode}", file=sys.stderr)
            if result.stderr:
                print(f"Error message: {result.stderr}", file=sys.stderr)
            return False
        
        # Check if output file was created
        if not os.path.exists(output_fasta):
            print(f"Error: Output file '{output_fasta}' was not created", file=sys.stderr)
            return False
        
        # Count sequences in output file
        sequence_count = 0
        with open(output_fasta, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    sequence_count += 1
        
        print(f"Successfully extracted {sequence_count} transcript sequences")
        print(f"Output written to: {output_fasta}")
        return True
        
    except Exception as e:
        print(f"Error running gffread: {e}", file=sys.stderr)
        return False


def parse_arguments():
    """Parse command-line arguments. Order: required GTF, genome, output; version."""
    parser = argparse.ArgumentParser(
        description='Extract transcript FASTA sequences from a genome FASTA file using gffread',
        usage='%(prog)s [-h] [-v,--version]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-g', '--gtf', action='store', dest='gtf_file', required=True,
                        metavar='', help='Filtered GTF (e.g. from remove_mirna_hairpin_from_gtf.py)')
    parser.add_argument('-f', '--genome-fasta', action='store', dest='genome_fasta', required=True,
                        metavar='', help='Genome FASTA (e.g. primary assembly from Ensembl)')
    parser.add_argument('-o', '--output', action='store', dest='output_fasta', required=True,
                        metavar='', help='Output transcript FASTA file')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    
    return parser.parse_args()


def main():
    """Main function to orchestrate the transcript extraction workflow."""
    args = parse_arguments()
    
    success = extract_transcripts_with_gffread(
        args.gtf_file,
        args.genome_fasta,
        args.output_fasta
    )
    
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()

