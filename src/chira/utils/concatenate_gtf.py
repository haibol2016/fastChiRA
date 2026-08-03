#!/usr/bin/env python
"""
Concatenate mature miRNA GTF/GFF3 file with target transcriptome GTF file.

This script combines:
1. Mature miRNA GTF/GFF3 file (e.g., from miRBase via download_mirbase_gff3.py) - with comment lines removed
   Note: miRBase provides GFF3 format, which ChiRA can handle directly. This script accepts GTF format.
2. Target transcriptome GTF file (output from remove_mirna_hairpin_from_gtf.py) - with miRNA hairpins removed

The output is a combined GTF file suitable for use as annotation file in ChiRA analysis.
"""

import argparse
import sys
import re

# Real-time output when stdout is not a TTY (e.g. batch jobs)
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)


def concatenate_gtf_files(mirna_gtf, target_gtf, output_gtf, keep_target_comments=True):
    """
    Concatenate mature miRNA GTF/GFF3 with target transcriptome GTF.
    
    Args:
        mirna_gtf: Path to mature miRNA GTF/GFF3 file (e.g., from miRBase via download_mirbase_gff3.py)
                   Note: miRBase provides GFF3 format, which ChiRA can handle directly. This script accepts GTF format.
        target_gtf: Path to target transcriptome GTF file (with miRNA hairpins removed)
        output_gtf: Path to output combined GTF file
        keep_target_comments: If True, keep comment lines from target GTF (default: True)
    """
    mirna_count = 0
    target_count = 0
    
    version_pattern = re.compile(r'transcript_id\s+"([^".]+)\.\d+"')
    
    try:
        with open(output_gtf, 'w') as fh_out:
            # First, write target transcriptome GTF (with optional comments)
            print(f"Reading target transcriptome GTF: {target_gtf}")
            with open(target_gtf, 'r') as fh_target:
                for line in fh_target:
                    if line.startswith('#'):
                        if keep_target_comments:
                            fh_out.write(line)
                    else:
                        # remove version number from transcript_id
                        line = version_pattern.sub(r'transcript_id "\1"', line)
                        fh_out.write(line)
                        target_count += 1
            
            # Then, append mature miRNA GTF (without comment lines)
            print(f"Reading mature miRNA GTF: {mirna_gtf}")
            with open(mirna_gtf, 'r') as fh_mirna:
                for line in fh_mirna:
                    # Skip comment lines from miRNA GTF
                    if line.startswith('#'):
                        continue
                    # remove version number from transcript_id
                    line = version_pattern.sub(r'transcript_id "\1"', line)
                    fh_out.write(line)
                    mirna_count += 1
        
        print(f"\nCombined GTF file created:")
        print(f"  - {target_count} lines from target transcriptome")
        print(f"  - {mirna_count} lines from mature miRNAs")
        print(f"  - Total: {target_count + mirna_count} lines")
        print(f"Output written to: {output_gtf}")
        
    except FileNotFoundError as e:
        print(f"Error: File not found: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing GTF files: {e}", file=sys.stderr)
        sys.exit(1)


def parse_arguments():
    """Parse command-line arguments. Order: required I/O, optional flags, version."""
    parser = argparse.ArgumentParser(
        description='Concatenate mature miRNA GTF/GFF3 file with target transcriptome GTF file',
        usage='%(prog)s [-h] [-v,--version]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required I/O (order: miRNA GTF, target GTF, output)
    parser.add_argument('-m', '--mirna-gtf', action='store', dest='mirna_gtf', required=True,
                        metavar='', help='Mature miRNA GTF/GFF3 (e.g. from download_mirbase_gff3.py). GFF3 works with ChiRA.')
    parser.add_argument('-t', '--target-gtf', action='store', dest='target_gtf', required=True,
                        metavar='', help='Target transcriptome GTF (e.g. from remove_mirna_hairpin_from_gtf.py)')
    parser.add_argument('-o', '--output', action='store', dest='output_gtf', required=True,
                        metavar='', help='Output combined GTF file')
    parser.add_argument('--remove-target-comments', action='store_true', dest='remove_target_comments',
                        help='Remove comment lines from target GTF (default: keep)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    
    return parser.parse_args()


def main():
    """Main function to orchestrate the GTF concatenation workflow."""
    args = parse_arguments()
    
    concatenate_gtf_files(
        args.mirna_gtf,
        args.target_gtf,
        args.output_gtf,
        keep_target_comments=not args.remove_target_comments
    )


if __name__ == "__main__":
    main()

