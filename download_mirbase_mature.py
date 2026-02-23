#!/usr/bin/env python
"""
Download species-specific mature miRNA sequences from miRBase.

This script downloads the mature miRNA FASTA file from miRBase and extracts
sequences for a specific species based on the species code (e.g., hsa for 
human, mmu for mouse, bta for bovine, rno for rat, etc.).
"""

import argparse
import sys
import os
import requests

# Real-time output when stdout is not a TTY (e.g. batch jobs)
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)


def download_mirbase_mature(output_file, version=None, timeout=30):
    """
    Download mature miRNA FASTA file from miRBase.
    
    Args:
        output_file: Path to save the downloaded file (uncompressed FASTA)
        version: miRBase version (e.g., "22.1"). If None, downloads CURRENT version
        timeout: Timeout in seconds for download (default: 30)
    
    Returns:
        True if successful, False otherwise
    """
    if version:
        # Use specific version URL
        url = f"https://www.mirbase.org/download_version_files/{version}/mature.fa"
    else:
        # Use CURRENT version
        url = "https://www.mirbase.org/download/mature.fa"

    try:
        print(f"Downloading from: {url}")
        with requests.get(url, stream=True, timeout=timeout) as r:
            r.raise_for_status()
            
            total_size = int(r.headers.get('content-length', 0))
            downloaded = 0
            first_chunk = True
            
            with open(output_file, 'w', encoding='utf-8') as f:
                for chunk in r.iter_content(chunk_size=8192, decode_unicode=True):
                    if chunk:
                        if first_chunk:
                            # Validate it's FASTA (ensure chunk is string, not bytes)
                            chunk_str = chunk if isinstance(chunk, str) else chunk.decode('utf-8')
                            if not chunk_str.strip().startswith('>'):
                                print(f"Error: Not a valid FASTA file (starts with: {chunk_str[:100]})", file=sys.stderr)
                                return False
                            first_chunk = False
                            chunk = chunk_str
                        
                        # Ensure chunk is string for writing
                        if isinstance(chunk, bytes):
                            chunk = chunk.decode('utf-8')
                        
                        f.write(chunk)
                        downloaded += len(chunk.encode('utf-8'))
                        
                        if total_size > 0:
                            percent = (downloaded / total_size) * 100
                            print(f"\rProgress: {percent:.1f}% ({downloaded}/{total_size} bytes)", end='', flush=True)
            
            if total_size > 0:
                print()
            print(f"Successfully downloaded to: {output_file}")
            return True
                    
    except requests.exceptions.RequestException as e:
        print(f"Failed to download: {e}", file=sys.stderr)
        return False


def extract_species_mirnas(input_file, output_file, species_code):
    """
    Extract species-specific mature miRNA sequences from miRBase FASTA file.
    
    Converts all "U" (uracil) nucleotides to "T" (thymine) to be compatible with
    ChiRA analysis, which expects DNA sequences rather than RNA sequences.
    
    Args:
        input_file: Path to input FASTA file (uncompressed)
        output_file: Path to output FASTA file (species-specific sequences)
        species_code: Species code (e.g., 'hsa' for human, 'mmu' for mouse, 'bta' for bovine)
    
    Returns:
        Number of sequences extracted
    """
    sequence_count = 0
    current_header = None
    current_sequence = []
    in_species_sequence = False
    
    try:
        with open(input_file, 'r', encoding='utf-8') as fh_in, open(output_file, 'w') as fh_out:
            for line in fh_in:
                line = line.rstrip('\n\r')
                
                if line.startswith('>'):
                    # Write previous sequence if it was species-specific
                    if in_species_sequence and current_header and current_sequence:
                        # Convert U to T in sequence (both uppercase and lowercase)
                        sequence_str = ''.join(current_sequence).replace('U', 'T').replace('u', 't')
                        fh_out.write(f"{current_header}\n{sequence_str}\n")
                        sequence_count += 1
                    
                    # Check if this header is for the target species
                    # Pattern: >species_code- (e.g., >hsa-, >mmu-, >bta-)
                    current_header = line
                    in_species_sequence = f'>{species_code}-' in line
                    current_sequence = []
                elif in_species_sequence:
                    current_sequence.append(line)
            
            # Write last sequence if it was species-specific
            if in_species_sequence and current_header and current_sequence:
                # Convert U to T in sequence (both uppercase and lowercase)
                sequence_str = ''.join(current_sequence).replace('U', 'T').replace('u', 't')
                fh_out.write(f"{current_header}\n{sequence_str}\n")
                sequence_count += 1
        
        return sequence_count
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing FASTA file: {e}", file=sys.stderr)
        sys.exit(1)


def cleanup_temp_file(temp_file, keep_file):
    """
    Remove temporary file if not keeping it.
    
    Args:
        temp_file: Path to temporary file
        keep_file: If True, keep the file; if False, remove it
    """
    if not keep_file:
        try:
            os.remove(temp_file)
            print(f"Removed temporary file: {temp_file}")
        except Exception as e:
            print(f"Warning: Could not remove temporary file: {e}", file=sys.stderr)


def parse_arguments():
    """Parse command-line arguments. Order: required species/output, optional version/flags, version."""
    parser = argparse.ArgumentParser(
        description='Download species-specific mature miRNA sequences from miRBase',
        usage='%(prog)s [-h] [-v,--version]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-s', '--species', action='store', dest='species_code', required=True,
                        metavar='', help='Species code (e.g., hsa=human, mmu=mouse, bta=bovine, rno=rat)')
    parser.add_argument('-o', '--output', action='store', dest='output_file', required=True,
                        metavar='', help='Output FASTA for species-specific mature miRNAs')
    parser.add_argument('--mirbase-version', action='store', dest='mirbase_version', default=None,
                        metavar='', help='miRBase version (e.g., "22.1"). Default: CURRENT')
    parser.add_argument('--keep-full', action='store_true', dest='keep_full',
                        help='Keep full mature.fa after extraction (default: remove)')
    parser.add_argument('--timeout', action='store', type=int, dest='timeout', default=30,
                        metavar='', help='Download timeout in seconds (default: 30)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    
    return parser.parse_args()


def main():
    """Main function to orchestrate the miRBase mature miRNA download workflow."""
    args = parse_arguments()
    
    temp_file = 'mature.fa'
    
    print(f"Downloading mature miRNA sequences from miRBase...")
    print(f"Version: {args.mirbase_version or 'CURRENT'}")
    
    if not download_mirbase_mature(temp_file, args.mirbase_version, args.timeout):
        print("Error: Failed to download mature miRNA sequences from miRBase", file=sys.stderr)
        sys.exit(1)
    
    print(f"\nExtracting {args.species_code} mature miRNAs...")
    sequence_count = extract_species_mirnas(temp_file, args.output_file, args.species_code)
    
    print(f"Extracted {sequence_count} mature miRNA sequences for {args.species_code}")
    print(f"Output written to: {args.output_file}")
    
    cleanup_temp_file(temp_file, args.keep_full)


if __name__ == "__main__":
    main()

