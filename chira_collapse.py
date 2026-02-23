#!/usr/bin/env python
from collections import defaultdict
import argparse
import gzip
import chira_utilities

# Note: Bio.SeqIO is no longer used - raw file parsing is faster for large files


def print_configuration(args):
    """Print configuration summary."""
    print('Input FASTQ          : ' + args.fastq)
    print('Output FASTA         : ' + args.fasta)
    print('Length of the UMI    : ' + str(args.umi_len))


def collapse_fastq_to_fasta(fastq_file, fasta_file, umi_len):
    """
    Collapse FASTQ reads to FASTA format, deduplicating sequences.
    
    Args:
        fastq_file: Path to input FASTQ file (gzipped or not)
        fasta_file: Path to output FASTA file
        umi_len: Length of UMI to trim from 5' end (0 if no UMI)
    """
    d_uniq_reads = defaultdict(int)
    # Cache UMI length check to avoid repeated conditionals
    has_umi = umi_len > 0
    
    # OPTIMIZATION: Use larger buffer size (2MB) for I/O efficiency with large FASTQ files
    # This significantly reduces system calls and improves performance
    # 2MB is a good balance between memory usage and I/O efficiency
    BUFFER_SIZE = 2 * 1024 * 1024  # 2MB buffer
    
    # Use raw file parsing for better performance with large files
    # This is significantly faster than SeqIO.parse() for very large FASTQ files
    # FASTQ format: @header (line 1), sequence (line 2), +header (line 3), quality (line 4)
    # Note: gzip.open() already has internal buffering, so we only apply explicit buffering to regular files
    if fastq_file.endswith('.gz'):
        # gzip.open() in text mode already uses internal buffering for decompression
        # Additional buffering layer is not needed and could complicate the code
        fh_fastq = gzip.open(fastq_file, 'rt')
    else:
        # For regular files, use buffering parameter to reduce system calls
        fh_fastq = open(fastq_file, 'r', buffering=BUFFER_SIZE)
    
    with fh_fastq:
        line_num = 0
        for line in fh_fastq:
            line_num += 1
            # Sequence is on every 2nd line of each 4-line block (lines 2, 6, 10, ...)
            if line_num % 4 == 2:  # Sequence line
                seq_str = line.rstrip('\n\r')
                if has_umi:
                    umi = seq_str[:umi_len]
                    sequence = seq_str[umi_len:]
                else:
                    umi = ""
                    sequence = seq_str
                d_uniq_reads[(sequence, umi)] += 1
    
    # Write output with optimized string formatting
    # OPTIMIZATION: Use larger buffer size (2MB) for write operations
    # This reduces write system calls and improves performance for large output files
    c = 1
    with open(fasta_file, "w", buffering=BUFFER_SIZE) as fh_out:
        # Sort items once and iterate directly
        # Use f-strings for better performance than concatenation
        for (sequence, umi), readcount in sorted(d_uniq_reads.items()):
            if umi:
                seqid = f"{c}|{umi}|{readcount}"
            else:
                seqid = f"{c}|{readcount}"
            fh_out.write(f">{seqid}\n{sequence}\n")
            c += 1


def parse_arguments():
    """Parse command-line arguments. Order: required I/O, optional UMI, version."""
    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: collapse FASTQ reads to FASTA format',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required I/O
    parser.add_argument('-i', '--fastq', action='store', dest='fastq', required=True, metavar='',
                        help='Input FASTQ file (gzipped or not)')
    parser.add_argument('-o', '--fasta', action='store', dest='fasta', required=True, metavar='',
                        help='Output FASTA file')
    parser.add_argument("-u", '--umi_len', action='store', type=int, default=0, metavar='',
                        help="""Length of the UMI if present. Trimmed from the 5' end of each read and appended to the tag id.""")
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {chira_utilities.__version__}')

    return parser.parse_args()


def main():
    """Main function to orchestrate the FASTQ to FASTA collapse workflow."""
    args = parse_arguments()
    print_configuration(args)
    collapse_fastq_to_fasta(args.fastq, args.fasta, args.umi_len)


if __name__ == "__main__":
    main()
