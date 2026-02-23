#!/usr/bin/env python
"""
Download GTF, and genome FASTA files from Ensembl.

This script downloads species-specific files from Ensembl FTP server
for a given genome version and GTF annotation version.
"""

import argparse
import sys
import os
import requests
import gzip
import shutil
from ftplib import FTP
from urllib.parse import urljoin
import time

# Real-time output when stdout is not a TTY (e.g. batch jobs)
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)


def download_file_http(url, local_path, timeout=60):
    """Download a file from HTTP/HTTPS URL with progress indication."""
    print(f"Downloading from: {url}")
    try:
        with requests.get(url, stream=True, timeout=timeout) as r:
            r.raise_for_status()
            total_size = int(r.headers.get('content-length', 0))
            downloaded = 0
            
            with open(local_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0:
                            percent = (downloaded / total_size) * 100
                            print(f"\rProgress: {percent:.1f}% ({downloaded}/{total_size} bytes)", end='', flush=True)
            
            if total_size > 0:
                print()  # New line after progress
            print(f"Successfully downloaded {local_path}")
            return True
    except requests.exceptions.RequestException as e:
        print(f"HTTP download failed: {e}", file=sys.stderr)
        return False


def download_file_ftp(host, remote_path, local_path, timeout=60):
    """Download a file from FTP server with progress indication."""
    print(f"Downloading from: ftp://{host}{remote_path}")
    try:
        with FTP(host, timeout=timeout) as ftp:
            ftp.login()
            
            # Get file size for progress
            try:
                file_size = ftp.size(remote_path)
            except:
                file_size = 0
            
            downloaded = 0
            with open(local_path, 'wb') as f:
                def callback(data):
                    nonlocal downloaded
                    f.write(data)
                    downloaded += len(data)
                    if file_size > 0:
                        percent = (downloaded / file_size) * 100
                        print(f"\rProgress: {percent:.1f}% ({downloaded}/{file_size} bytes)", end='', flush=True)
                
                ftp.retrbinary(f"RETR {remote_path}", callback)
            
            if file_size > 0:
                print()  # New line after progress
            print(f"Successfully downloaded {local_path}")
            return True
    except Exception as e:
        print(f"FTP download failed: {e}", file=sys.stderr)
        return False


def decompress_file(gz_file, output_file):
    """Decompress a gzipped file."""
    print(f"Decompressing {gz_file}...")
    try:
        # Use text mode for FASTA and GTF files (they are text files)
        with gzip.open(gz_file, 'rt', encoding='utf-8') as f_in:
            with open(output_file, 'w', encoding='utf-8') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"Decompressed to {output_file}")
        return True
    except Exception as e:
        print(f"Error decompressing file: {e}", file=sys.stderr)
        return False


def list_ftp_directory(host, remote_dir):
    """List files in an FTP directory."""
    try:
        with FTP(host, timeout=30) as ftp:
            ftp.login()
            # Use NLST which returns just filenames (simpler than LIST)
            try:
                filenames = ftp.nlst(remote_dir)
                # Remove directory path prefix if present
                filenames = [os.path.basename(f) for f in filenames if f not in ['.', '..']]
                return filenames
            except:
                # Fallback to LIST if NLST doesn't work
                files = []
                ftp.retrlines(f'LIST {remote_dir}', files.append)
                # Parse file list to get filenames
                filenames = []
                for line in files:
                    parts = line.split()
                    if len(parts) >= 9:
                        filename = ' '.join(parts[8:])
                        if filename not in ['.', '..']:
                            filenames.append(filename)
                return filenames
    except Exception as e:
        print(f"Warning: Could not list FTP directory: {e}", file=sys.stderr)
        return []


def find_ensembl_file(host, remote_dir, pattern):
    """
    Find a file in an FTP directory matching a pattern.
    
    Args:
        host: FTP hostname
        remote_dir: Remote directory path
        pattern: Pattern to match (e.g., '.cdna.all.fa.gz')
    
    Returns:
        Filename if found, None otherwise
    """
    files = list_ftp_directory(host, remote_dir)
    for filename in files:
        if pattern in filename:
            return filename
    return None


def get_ensembl_file_info(species, release_version, file_type, assembly=None):
    """
    Get Ensembl file information (filename and path).
    
    Args:
        species: Species name (e.g., 'homo_sapiens', 'mus_musculus')
        release_version: Ensembl release version (e.g., '110')
        file_type: Type of file ('cdna', 'ncrna', 'gtf', 'genome')
        assembly: Optional assembly name (e.g., 'GRCh38'). If not provided, will try to find it.
    
    Returns:
        Tuple of (filename, remote_path) or None if not found
    """
    base_path = f"/pub/release-{release_version}"
    species_cap = species.capitalize()
    
    if file_type == 'gtf':
        remote_dir = f"{base_path}/gtf/{species}/"
        pattern = f'.{release_version}.gtf.gz'
    elif file_type == 'genome':
        remote_dir = f"{base_path}/fasta/{species}/dna/"
        pattern = '.dna.primary_assembly.fa.gz'
    else:
        return None
    
    # Try to find the file by listing the directory
    ftp_host = "ftp.ensembl.org"
    filename = find_ensembl_file(ftp_host, remote_dir, pattern)
    
    if filename:
        remote_path = remote_dir + filename
        return (filename, remote_path)
    
    # Fallback: construct filename if assembly is provided
    if assembly:
        if file_type == 'gtf':
            filename = f"{species_cap}.{assembly}.{release_version}.gtf.gz"
        elif file_type == 'genome':
            filename = f"{species_cap}.{assembly}.dna.primary_assembly.fa.gz"
        remote_path = remote_dir + filename
        return (filename, remote_path)
    
    return None


def download_ensembl_files(species, genome_version, gtf_version, output_dir, 
                           assembly=None, keep_compressed=False, decompress=True, timeout=60):
    """
    Download Ensembl files for a given species and versions.
    
    Args:
        species: Species name (e.g., 'homo_sapiens', 'mus_musculus')
        genome_version: Ensembl release version for genome/cDNA/ncRNA (e.g., '110')
        gtf_version: Ensembl release version for GTF annotation (e.g., '110')
        output_dir: Directory to save downloaded files
        assembly: Optional assembly name (e.g., 'GRCh38'). If not provided, will try to auto-detect.
        keep_compressed: If True, keep compressed files after decompression
        decompress: If True, decompress gzipped files
        timeout: Download timeout in seconds
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Ensembl FTP server
    ftp_host = "ftp.ensembl.org"
    
    files_to_download = [
        ('gtf', gtf_version),
        ('genome', genome_version)
    ]
    
    downloaded_files = []
    
    for file_type, version in files_to_download:
        print(f"\n{'='*60}")
        print(f"Downloading {file_type.upper()} file (release {version})...")
        print(f"{'='*60}")
        
        file_info = get_ensembl_file_info(species, version, file_type, assembly)
        if not file_info:
            print(f"Error: Could not find {file_type} file for species {species}, version {version}", file=sys.stderr)
            if not assembly:
                print(f"  Hint: Try specifying --assembly if auto-detection fails", file=sys.stderr)
            continue
        
        filename, remote_path = file_info
        local_gz_file = os.path.join(output_dir, filename)
        local_file = local_gz_file[:-3] if local_gz_file.endswith('.gz') else local_gz_file
        
        # Download file
        download_success = False
        
        # Try HTTPS first
        https_url = f"https://{ftp_host}{remote_path}"
        if download_file_http(https_url, local_gz_file, timeout):
            download_success = True
        else:
            # Fallback to FTP
            print("HTTPS download failed, trying FTP...", file=sys.stderr)
            if download_file_ftp(ftp_host, remote_path, local_gz_file, timeout):
                download_success = True
        
        if not download_success:
            print(f"Failed to download {file_type} file", file=sys.stderr)
            continue
        
        # Decompress if requested
        if decompress and filename.endswith('.gz'):
            if decompress_file(local_gz_file, local_file):
                downloaded_files.append(local_file)
                # Remove compressed file if not keeping it
                if not keep_compressed:
                    os.remove(local_gz_file)
                    print(f"Removed compressed file: {local_gz_file}")
            else:
                print(f"Warning: Decompression failed, keeping compressed file: {local_gz_file}", file=sys.stderr)
                downloaded_files.append(local_gz_file)
        else:
            downloaded_files.append(local_gz_file)
    
    print(f"\n{'='*60}")
    print("Download Summary:")
    print(f"{'='*60}")
    for f in downloaded_files:
        if os.path.exists(f):
            size = os.path.getsize(f)
            size_mb = size / (1024 * 1024)
            print(f"  {os.path.basename(f)}: {size_mb:.2f} MB")
        else:
            print(f"  {os.path.basename(f)}: NOT FOUND")
    
    return downloaded_files


def parse_arguments():
    """Parse command-line arguments. Order: required species/version/output, optional assembly/flags, version."""
    parser = argparse.ArgumentParser(
        description='Download GTF and genome FASTA files from Ensembl',
        usage='%(prog)s [-h] [-v,--version]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-s', '--species', action='store', dest='species', required=True,
                        metavar='', help='Species (e.g., homo_sapiens, mus_musculus, bos_taurus)')
    parser.add_argument('-g', '--genome-version', action='store', dest='genome_version', required=True,
                        metavar='', help='Ensembl release versionfor genome (e.g., 110)')
    parser.add_argument('-t', '--gtf-version', action='store', dest='gtf_version', required=True,
                        metavar='', help='Ensembl release version for GTF (e.g., 110)')
    parser.add_argument('-o', '--output-dir', action='store', dest='output_dir', required=True,
                        metavar='', help='Output directory for downloaded files')
    parser.add_argument('-a', '--assembly', action='store', dest='assembly', default=None,
                        metavar='', help='Genome assembly (e.g., GRCh38). Auto-detect if not set.')
    parser.add_argument('--keep-compressed', action='store_true', dest='keep_compressed',
                        help='Keep compressed files after decompression (default: remove)')
    parser.add_argument('--no-decompress', action='store_true', dest='no_decompress',
                        help='Do not decompress gzipped files (default: decompress)')
    parser.add_argument('--timeout', action='store', type=int, default=60,
                        metavar='', help='Download timeout in seconds')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    
    return parser.parse_args()


def main():
    """Main function to orchestrate the Ensembl download workflow."""
    args = parse_arguments()
    
    # Normalize species name (lowercase, replace spaces with underscores)
    species = args.species.lower().replace(' ', '_')
    
    download_ensembl_files(
        species,
        args.genome_version,
        args.gtf_version,
        args.output_dir,
        assembly=args.assembly,
        keep_compressed=args.keep_compressed,
        decompress=not args.no_decompress,
        timeout=args.timeout
    )


if __name__ == "__main__":
    main()

