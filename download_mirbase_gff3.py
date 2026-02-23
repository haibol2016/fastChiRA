#!/usr/bin/env python
"""
Download species-specific GFF3 files from miRBase.

This script downloads the GFF3 file from miRBase for a specific species
based on the species code (e.g., hsa for human, mmu for mouse, bta for 
bovine, rno for rat, etc.).

The GFF3 file contains chromosomal coordinates of microRNAs including:
- miRNA_primary_transcript (hairpin precursor sequences)
- miRNA (mature sequences): ChiRA can make use the mature sequences for the analysis
- You don't need to convert the GFF3 file to GTF file for the analysis.

Optional features:
- Coordinate liftover: Convert coordinates between genome assemblies using pyliftover
- Chromosome renaming: Rename chromosomes based on a mapping file
- Processing order: Download -> Liftover -> Rename chromosomes
"""

import argparse
import sys
import os
import requests

# Real-time output when stdout is not a TTY (e.g. batch jobs)
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)

# Optional dependency for coordinate liftover
try:
    from pyliftover import LiftOver
    PYLIFTOVER_AVAILABLE = True
except ImportError:
    PYLIFTOVER_AVAILABLE = False


def read_chromosome_mapping(mapping_file):
    """
    Read chromosome mapping from a tab-separated file.
    
    Args:
        mapping_file: Path to tab-separated file with two columns:
                      gff3_chromosome_name<tab>target_chromosome_name
    
    Returns:
        Dictionary mapping gff3 chromosome names to target chromosome names
    """
    chr_mapping = {}
    try:
        with open(mapping_file, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 2:
                    print(f"Warning: Skipping line {line_num} in mapping file (expected 2 columns, got {len(fields)})", 
                          file=sys.stderr)
                    continue
                
                gff3_chr = fields[0].strip()
                target_chr = fields[1].strip()
                
                if gff3_chr and target_chr:
                    chr_mapping[gff3_chr] = target_chr
                else:
                    print(f"Warning: Skipping line {line_num} in mapping file (empty chromosome name)", 
                          file=sys.stderr)
        
        print(f"Loaded {len(chr_mapping)} chromosome mappings from {mapping_file}")
        return chr_mapping
        
    except FileNotFoundError:
        print(f"Error: Mapping file '{mapping_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading mapping file: {e}", file=sys.stderr)
        sys.exit(1)


def perform_coordinate_liftover(input_file, output_file, chain_file, source_genome, target_genome):
    """
    Perform coordinate liftover on a GFF3 file using pyliftover.
    
    Args:
        input_file: Path to input GFF3 file
        output_file: Path to output GFF3 file with lifted coordinates
        chain_file: Path to chain file for coordinate conversion
        source_genome: Source genome assembly name (e.g., 'hg19', 'hg38/GRCh38')
        target_genome: Target genome assembly name (e.g., 'hg38/GRCh38', 'hg19')
    
    Returns:
        Tuple of (success: bool, converted_count: int, failed_count: int)
    """
    if not PYLIFTOVER_AVAILABLE:
        print("Error: pyliftover package is required for coordinate liftover. "
              "Install it with: pip install pyliftover", file=sys.stderr)
        return False, 0, 0
    
    try:
        # Initialize LiftOver converter
        print(f"Loading chain file: {chain_file}")
        lo = LiftOver(chain_file)
        
        converted_count = 0
        failed_count = 0
        total_lines = 0
        
        with open(input_file, 'r', encoding='utf-8') as f_in, \
             open(output_file, 'w', encoding='utf-8') as f_out:
            
            for line in f_in:
                total_lines += 1
                
                # Skip comment lines and empty lines
                if line.startswith('#') or not line.strip():
                    f_out.write(line)
                    continue
                
                # GFF3 format: tab-separated
                # Columns: seqid, source, type, start, end, score, strand, phase, attributes
                fields = line.rstrip('\n\r').split('\t')
                
                if len(fields) >= 9:
                    seqid = fields[0]
                    start_pos = fields[3]  # 1-based start (GFF3 uses 1-based coordinates)
                    end_pos = fields[4]    # 1-based end
                    strand = fields[6]
                    
                    # Convert coordinates (pyliftover uses 0-based coordinates)
                    try:
                        start_int = int(start_pos)
                        end_int = int(end_pos)
                        
                        # Convert to 0-based for pyliftover
                        start_0based = start_int - 1
                        end_0based = end_int - 1
                        
                        # Perform liftover for start position
                        start_converted = lo.convert_coordinate(seqid, start_0based)
                        # Perform liftover for end position
                        end_converted = lo.convert_coordinate(seqid, end_0based)
                        
                        if start_converted and end_converted and len(start_converted) > 0 and len(end_converted) > 0:
                            # Use first conversion result (pyliftover may return multiple)
                            new_start_0based = start_converted[0][1]
                            new_end_0based = end_converted[0][1]
                            new_seqid = start_converted[0][0]
                            
                            # Convert back to 1-based for GFF3
                            new_start = new_start_0based + 1
                            new_end = new_end_0based + 1
                            
                            # Validate: ensure start <= end
                            if new_start > new_end:
                                # Swap if needed (shouldn't happen, but handle edge cases)
                                new_start, new_end = new_end, new_start
                            
                            # Update fields
                            fields[0] = new_seqid  # Update chromosome name if changed
                            fields[3] = str(new_start)
                            fields[4] = str(new_end)
                            
                            converted_count += 1
                        else:
                            # Liftover failed for this feature
                            failed_count += 1
                            # Keep original coordinates if liftover fails
                            pass
                            
                    except (ValueError, IndexError) as e:
                        # Invalid coordinates, skip conversion
                        failed_count += 1
                        pass
                    
                    # Write the line (with or without conversion)
                    f_out.write('\t'.join(fields) + '\n')
                else:
                    # Write line as-is if it doesn't have expected format
                    f_out.write(line)
        
        print(f"Liftover completed: {converted_count} features converted, {failed_count} failed out of {total_lines} total lines")
        return True, converted_count, failed_count
        
    except FileNotFoundError:
        print(f"Error: Chain file '{chain_file}' not found", file=sys.stderr)
        return False, 0, 0
    except Exception as e:
        print(f"Error performing coordinate liftover: {e}", file=sys.stderr)
        return False, 0, 0


def convert_chromosome_names(input_file, output_file, chr_mapping):
    """
    Convert chromosome names in a GFF3 file using a mapping dictionary.
    
    Args:
        input_file: Path to input GFF3 file
        output_file: Path to output GFF3 file with converted chromosome names
        chr_mapping: Dictionary mapping source chromosome names to target names
    
    Returns:
        Number of lines converted
    """
    converted_count = 0
    total_lines = 0
    
    try:
        with open(input_file, 'r', encoding='utf-8') as f_in, \
             open(output_file, 'w', encoding='utf-8') as f_out:
            
            for line in f_in:
                total_lines += 1
                
                # Skip comment lines and empty lines
                if line.startswith('#') or not line.strip():
                    f_out.write(line)
                    continue
                
                # GFF3 format: tab-separated, first column is chromosome/sequence name
                fields = line.rstrip('\n\r').split('\t')
                
                if len(fields) >= 1:
                    original_chr = fields[0]
                    
                    # Apply mapping if available
                    if original_chr in chr_mapping:
                        fields[0] = chr_mapping[original_chr]
                        converted_count += 1
                    
                    # Write the line (with or without conversion)
                    f_out.write('\t'.join(fields) + '\n')
                else:
                    # Write line as-is if it doesn't have expected format
                    f_out.write(line)
        
        print(f"Converted {converted_count} chromosome names out of {total_lines} total lines")
        return converted_count
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error converting chromosome names: {e}", file=sys.stderr)
        sys.exit(1)


def download_mirbase_gff3(output_file, species_code, version=None, timeout=600, 
                          chr_mapping=None, liftover_params=None):
    """
    Download species-specific GFF3 file from miRBase.
    
    Args:
        output_file: Path to save the downloaded file
        species_code: Species code (e.g., 'hsa' for human, 'mmu' for mouse, 'bta' for bovine)
        version: miRBase version number (e.g., "21"). If None, downloads CURRENT version
        timeout: Timeout in seconds for download (default: 60)
        chr_mapping: Optional dictionary mapping source chromosome names to target names
        liftover_params: Optional dict with keys: chain_file, source_genome, target_genome
    
    Returns:
        True if successful, False otherwise
    """
    if version:
        # Use specific version URL
        # Pattern: https://www.mirbase.org/download_version_genome_files/{version}/{species_code}.gff3
        url = f"https://www.mirbase.org/download_version_genome_files/{version}/{species_code}.gff3"
    else:
        # Use CURRENT version
        # Pattern: https://www.mirbase.org/download/{species_code}.gff3
        url = f"https://www.mirbase.org/download/{species_code}.gff3"

    try:
        print(f"Downloading from: {url}")
        with requests.get(url, stream=True, timeout=timeout) as r:
            r.raise_for_status()
            
            # Check if we got a valid response (not 404, etc.)
            if r.status_code != 200:
                print(f"Error: HTTP {r.status_code} - {r.reason}", file=sys.stderr)
                return False
            
            total_size = int(r.headers.get('content-length', 0))
            downloaded = 0
            first_chunk = True
            
            with open(output_file, 'w', encoding='utf-8') as f:
                for chunk in r.iter_content(chunk_size=8192, decode_unicode=True):
                    if chunk:
                        if first_chunk:
                            # Validate it's GFF3 (check for GFF version header)
                            chunk_str = chunk if isinstance(chunk, str) else chunk.decode('utf-8')
                            if not (chunk_str.startswith('##gff-version') or 
                                    chunk_str.startswith('#') or
                                    'gff' in chunk_str.lower()[:100]):
                                print(f"Warning: File may not be a valid GFF3 file (starts with: {chunk_str[:100]})", 
                                      file=sys.stderr)
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
            
            # Process the file: liftover first, then chromosome renaming
            current_file = output_file
            temp_files = []
            
            # Step 1: Perform coordinate liftover if requested
            if liftover_params:
                temp_file = output_file + '.liftover.tmp'
                temp_files.append(temp_file)
                print(f"\nPerforming coordinate liftover...")
                print(f"  Source: {liftover_params['source_genome']}")
                print(f"  Target: {liftover_params['target_genome']}")
                print(f"  Chain file: {liftover_params['chain_file']}")
                
                success, converted, failed = perform_coordinate_liftover(
                    current_file, temp_file,
                    liftover_params['chain_file'],
                    liftover_params['source_genome'],
                    liftover_params['target_genome']
                )
                
                if not success:
                    print("Error: Coordinate liftover failed", file=sys.stderr)
                    # Clean up temp files
                    for tf in temp_files:
                        if os.path.exists(tf):
                            os.remove(tf)
                    return False
                
                # Update current file to the lifted version
                if current_file != output_file:
                    os.remove(current_file)
                current_file = temp_file
            
            # Step 2: Apply chromosome name mapping if provided
            if chr_mapping:
                final_file = output_file
                if current_file == output_file:
                    # Need to create temp file for renaming
                    temp_file = output_file + '.rename.tmp'
                    temp_files.append(temp_file)
                    # Rename current file to temp file
                    os.rename(current_file, temp_file)
                    current_file = temp_file
                else:
                    temp_file = current_file
                
                print(f"\nApplying chromosome name mapping...")
                convert_chromosome_names(temp_file, final_file, chr_mapping)
                
                # Clean up intermediate files
                if temp_file != output_file and os.path.exists(temp_file):
                    os.remove(temp_file)
            elif current_file != output_file:
                # Liftover was performed but no renaming, move temp file to final location
                os.rename(current_file, output_file)
            
            # Clean up any remaining temp files
            for tf in temp_files:
                if os.path.exists(tf) and tf != output_file:
                    try:
                        os.remove(tf)
                    except OSError:
                        pass
            
            if liftover_params or chr_mapping:
                print(f"Processing complete. Output written to: {output_file}")
            
            return True
                    
    except requests.exceptions.RequestException as e:
        print(f"Failed to download: {e}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return False


def validate_liftover_params(args):
    """
    Validate and prepare liftover parameters.
    
    Args:
        args: Parsed arguments
    
    Returns:
        dict: Liftover parameters dict or None if not provided
    """
    if args.chain_file:
        if not args.source_genome or not args.target_genome:
            print("Error: --source-genome and --target-genome are required when --chain-file is provided", 
                  file=sys.stderr)
            sys.exit(1)
        if not os.path.exists(args.chain_file):
            print(f"Error: Chain file '{args.chain_file}' not found", file=sys.stderr)
            sys.exit(1)
        return {
            'chain_file': args.chain_file,
            'source_genome': args.source_genome,
            'target_genome': args.target_genome
        }
    return None


def print_download_info(args, liftover_params, chr_mapping):
    """Print download information."""
    print(f"Downloading GFF3 file from miRBase...")
    print(f"Species: {args.species_code}")
    print(f"Version: {args.mirbase_version or 'CURRENT'}")
    if liftover_params:
        print(f"Coordinate liftover: {liftover_params['source_genome']} -> {liftover_params['target_genome']}")
    if chr_mapping:
        print(f"Chromosome mapping: {args.chr_mapping_file}")


def parse_arguments():
    """Parse command-line arguments. Order: required I/O, optional version/timeout, chromosome mapping, liftover, version."""
    parser = argparse.ArgumentParser(
        description='Download species-specific GFF3 file from miRBase',
        usage='%(prog)s [-h] [-v,--version]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required I/O
    parser.add_argument('-s', '--species', action='store', dest='species_code', required=True,
                        metavar='', help='Species code (e.g., hsa=human, mmu=mouse, bta=bovine, rno=rat)')
    parser.add_argument('-o', '--output', action='store', dest='output_file', required=True,
                        metavar='', help='Output GFF3 file for species-specific microRNA annotations')
    parser.add_argument('--mirbase-version', action='store', dest='mirbase_version', default=None,
                        metavar='', help='miRBase version (e.g., "21"). Default: CURRENT')
    parser.add_argument('--timeout', action='store', type=int, dest='timeout', default=60,
                        metavar='', help='Download timeout in seconds (default: 60)')
    parser.add_argument('-m', '--chromosome-mapping', action='store', dest='chr_mapping_file', default=None,
                        metavar='',
                        help="""Tab-separated file: gff3_chromosome_name<tab>target_chromosome_name.
Applied after liftover (if performed).""")
    parser.add_argument('--source-genome', action='store', dest='source_genome', default=None,
                        metavar='',
                        help="""Source genome assembly for liftover (e.g., hg19, mm9). Required if --chain-file given.""")
    parser.add_argument('--target-genome', action='store', dest='target_genome', default=None,
                        metavar='',
                        help="""Target genome assembly for liftover (e.g., hg38, mm10). Required if --chain-file given.""")
    parser.add_argument('--chain-file', action='store', dest='chain_file', default=None,
                        metavar='',
                        help="""Path to chain file for coordinate liftover (e.g. from UCSC).
Requires --source-genome and --target-genome. Needs pyliftover: pip install pyliftover""")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

    return parser.parse_args()


def main():
    """Main function to orchestrate the miRBase GFF3 download workflow."""
    args = parse_arguments()
    
    # Validate and prepare liftover parameters
    liftover_params = validate_liftover_params(args)
    
    # Read chromosome mapping if provided
    chr_mapping = None
    if args.chr_mapping_file:
        chr_mapping = read_chromosome_mapping(args.chr_mapping_file)
    
    # Print download information
    print_download_info(args, liftover_params, chr_mapping)
    
    # Download and process GFF3 file
    if not download_mirbase_gff3(args.output_file, args.species_code, args.mirbase_version, 
                                 args.timeout, chr_mapping, liftover_params):
        print("Error: Failed to download or process GFF3 file from miRBase", file=sys.stderr)
        sys.exit(1)
    
    print(f"Output written to: {args.output_file}")


if __name__ == "__main__":
    main()

