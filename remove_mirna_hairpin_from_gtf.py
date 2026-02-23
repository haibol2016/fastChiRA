#!/usr/bin/env python
"""
Remove all included microRNA entries from an Ensembl GTF file.

This script filters out all lines in a GTF file that are related to microRNAs,
including miRNA genes, miRNA transcripts, and miRNA features.

make sure to remove version number from transcript_id
"""

import argparse
import sys
import re

# Real-time output when stdout is not a TTY (e.g. batch jobs)
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)


def is_mirna_line(line, mirna_pattern=None):
    """
    Check if a GTF line is related to microRNA.
    
    Args:
        line: A line from a GTF file
        mirna_pattern: Optional compiled regex pattern for matching miRNA in attributes.
                      If None, uses default checks (feature type and biotype only).
    
    Returns:
        True if the line is related to microRNA, False otherwise
    """
    # Skip comment lines
    if line.startswith('#'):
        return False
    
    # Split the line into fields
    fields = line.strip().split('\t')
    if len(fields) < 9:
        return False
    
    # Check feature type
    source = fields[1]
    attributes = fields[8]
    return (source == 'miRBase' or mirna_pattern.search(attributes))


def remove_mirna_from_gtf(input_gtf, output_gtf, keep_comments=True, mirna_pattern=None):
    """
    Remove all microRNA entries from a GTF file.
    
    Args:
        input_gtf: Path to input GTF file
        output_gtf: Path to output GTF file
        keep_comments: If True, keep comment lines in output (default: True)
        mirna_pattern: Optional compiled regex pattern for matching miRNA in attributes.
                      If None, only checks feature type and biotype fields.
    """
    removed_count = 0
    kept_count = 0
    
    try:
        with open(input_gtf, 'r') as fh_in, open(output_gtf, 'w') as fh_out:
            for line in fh_in:
                # Keep comment lines if requested
                if line.startswith('#'):
                    if keep_comments:
                        fh_out.write(line)
                    continue
                
                # Check if line is miRNA-related
                if is_mirna_line(line, mirna_pattern):
                    removed_count += 1
                else:
                    fh_out.write(line)
                    kept_count += 1
        
        print(f"Removed {removed_count} microRNA-related lines")
        print(f"Kept {kept_count} non-microRNA lines")
        print(f"Output written to: {output_gtf}")
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_gtf}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing GTF file: {e}", file=sys.stderr)
        sys.exit(1)


def compile_mirna_pattern(pattern_str):
    """
    Compile miRNA pattern regex from string.
    
    Args:
        pattern_str: Regex pattern string or None for default pattern
    
    Returns:
        Compiled regex pattern
    """
    if pattern_str:
        try:
            return re.compile(pattern_str)
        except re.error as e:
            print(f"Error: Invalid regular expression pattern: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # Default pattern for miRNA detection
        return re.compile(r'gene_name\s+".*mir-.*"|transcript_name\s+".*mir-.*"|gene_biotype\s+"miRNA"|transcript_biotype\s+"miRNA"')


def parse_arguments():
    """Parse command-line arguments. Order: required I/O, optional pattern/flags, version."""
    parser = argparse.ArgumentParser(
        description='Remove all microRNA entries from an Ensembl GTF file',
        usage='%(prog)s [-h] [-v,--version]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required I/O
    parser.add_argument('-i', '--input', action='store', dest='input_gtf', required=True,
                        metavar='', help='Input GTF file')
    parser.add_argument('-o', '--output', action='store', dest='output_gtf', required=True,
                        metavar='', help='Output GTF file (without microRNA entries)')
    parser.add_argument('-p', '--pattern', action='store', dest='mirna_pattern', default=None,
                        metavar='',
                        help="""Regex for matching miRNA in GTF attributes. If not set, only feature type and biotype are used.
Example: 'gene_name\\s+"[^"]*(?:[Mm][Ii][Rr][_-]|[^"]*[-_][Mm][Ii][Rr][_-])'""")
    parser.add_argument('--remove-comments', action='store_true', dest='remove_comments',
                        help='Remove comment lines from output (default: keep comments)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    
    return parser.parse_args()


def main():
    """Main function to orchestrate the miRNA removal workflow."""
    args = parse_arguments()
    
    # Compile regex pattern
    mirna_pattern = compile_mirna_pattern(args.mirna_pattern)
    
    # Remove miRNA entries
    remove_mirna_from_gtf(
        args.input_gtf,
        args.output_gtf,
        keep_comments=not args.remove_comments,
        mirna_pattern=mirna_pattern
    )


if __name__ == "__main__":
    main()

