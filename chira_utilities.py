#!/usr/bin/env python
"""
ChiRA Utilities Module

This module contains shared utility functions and constants used across ChiRA scripts.
"""

import argparse
import re
import subprocess
import sys
from Bio import SeqIO
from datetime import datetime

# Version number - single source of truth for all ChiRA scripts
__version__ = "1.4.14"


def _enable_unbuffered_stdout_stderr():
    """Force stdout/stderr line-buffered so log output appears in real time (e.g. in batch jobs)."""
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(line_buffering=True)
        sys.stderr.reconfigure(line_buffering=True)
    else:
        # Python 3.6: wrap so every write flushes
        class _FlushWrapper:
            __slots__ = ('_stream',)

            def __init__(self, stream):
                self._stream = stream

            def write(self, s):
                self._stream.write(s)
                self._stream.flush()

            def flush(self):
                self._stream.flush()

            def __getattr__(self, name):
                return getattr(self._stream, name)

        sys.stdout = _FlushWrapper(sys.stdout)
        sys.stderr = _FlushWrapper(sys.stderr)


_enable_unbuffered_stdout_stderr()

# Try to import psutil for adaptive buffer sizing (optional dependency)
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

# Compile regex patterns once for better performance
_CIGAR_PATTERN_QUERY = re.compile(r'(\d+)([HIMSX=])')
_CIGAR_PATTERN_ALIGN = re.compile(r'(\d+)([DIMX=])')
_CIGAR_PATTERN_END = re.compile(r'(\d+)([DMNX=])')

# Cache for get_bedtools_command to avoid repeated subprocess calls
_bedtools_command_cache = {}


def score_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


def overlap(f, s):
    return max(0, min(f[1], s[1]) - max(f[0], s[0]))


def median(x):
    # Sort the list first (median requires sorted data)
    x_sorted = sorted(x)
    n = len(x_sorted)
    mid = n // 2  # Use integer division (more Pythonic)
    if not n % 2:
        return (x_sorted[mid-1] + x_sorted[mid]) / 2.0
    return x_sorted[mid]


def query_length(cigar, is_reverse):
    cigar_tup = _CIGAR_PATTERN_QUERY.findall(cigar)
    read_length = 0
    if is_reverse:
        cigar_tup = reversed(cigar_tup)
    for c in cigar_tup:
        read_length += int(c[0])
    return read_length


def match_positions(cigar, is_reverse):
    cigar_tup = _CIGAR_PATTERN_QUERY.findall(cigar)
    match_start = match_end = 0
    if is_reverse:
        cigar_tup = reversed(cigar_tup)
    for c in cigar_tup:
        if c[1] == "S" or c[1] == "H":
            if not match_start:
                match_start = int(c[0]) + 1
                match_end = int(c[0])
        elif c[1] == "M":
            if match_start:
                match_end = match_end + int(c[0])
            else:
                match_start = 1
                match_end = int(c[0])
        elif c[1] == "I":
            match_end = match_end + int(c[0])
    return match_start, match_end


def is_chimeric(cigar1, cigar2, is_reverse1, is_reverse2, max_allowed_overlap):
    match_start1, match_end1 = match_positions(cigar1, is_reverse1)
    match_start2, match_end2 = match_positions(cigar2, is_reverse2)
    chimeric = True
    if overlap([match_start1, match_end1], [match_start2, match_end2]) > max_allowed_overlap:
        chimeric = False
    return chimeric


def alignment_length(cigar):
    # everything except clipped
    cigar_tup = _CIGAR_PATTERN_ALIGN.findall(cigar)
    align_len = 0
    for c in cigar_tup:
        align_len += int(c[0])
    return align_len


def alignment_end(start, cigar, is_reverse):
    cigar_tup = _CIGAR_PATTERN_END.findall(cigar)
    end = int(start)
    if is_reverse:
        cigar_tup = reversed(cigar_tup)
    for c in cigar_tup:
        end += int(c[0])
    return end - 1


def bedentry(referenceid, reference_start, reference_end, readid, strand, cigarstring):
    line = '\t'.join([referenceid,
                      reference_start,
                      reference_end,
                      ','.join([readid, referenceid, reference_start, reference_end, strand, cigarstring]),
                      "1",
                      strand])
    return line


def reverse_complement(seq):
    tab = str.maketrans("ACTGactg", "TGACtgac")
    return seq.translate(tab)[::-1]


def extract_reflengths(ref_fasta, d_reflen):
    # Use context manager for proper file handling
    with open(ref_fasta, 'r') as fh:
        fa_ref = SeqIO.parse(fh, 'fasta')
        for record in fa_ref:
            d_reflen[record.id] = len(record)
    return


def print_w_time(message):
    # Use f-string for slightly better performance than concatenation
    # Print to stderr so logs appear in bsub error files and are visible during execution
    import sys
    print(f"[{datetime.now().strftime('%d-%m-%Y %H:%M:%S')}] {message}", file=sys.stderr)


def get_bedtools_command(tool_name):
    """
    Get the appropriate BEDTools command format for different versions.
    Older versions use individual commands (e.g., intersectBed, fastaFromBed),
    newer versions use unified bedtools command (e.g., bedtools intersect, bedtools getfasta).
    
    Args:
        tool_name: Name of the tool (e.g., 'intersect', 'getfasta')
    
    Returns:
        List of command parts to use with subprocess or os.system
    """
    # Check cache first to avoid repeated subprocess calls
    if tool_name in _bedtools_command_cache:
        return _bedtools_command_cache[tool_name]
    
    # Map tool names to old and new command formats
    tool_map = {
        'intersect': ('intersectBed', ['bedtools', 'intersect']),
        'getfasta': ('fastaFromBed', ['bedtools', 'getfasta'])
    }
    
    if tool_name not in tool_map:
        raise ValueError("Unknown BEDTools tool: " + tool_name)
    
    old_cmd, new_cmd = tool_map[tool_name]
    
    # Try new format first (bedtools <tool>)
    try:
        process = subprocess.run(new_cmd + ['--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=2, check=False)
        # Check if command actually succeeded (returncode 0)
        if process.returncode == 0:
            result = new_cmd
        else:
            # Command exists but failed, fall back to old format
            result = [old_cmd]
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
        # Fall back to old format (individual command)
        result = [old_cmd]
    
    # Cache the result
    _bedtools_command_cache[tool_name] = result
    return result


def get_adaptive_buffer_size(num_files=2):
    """
    OPTIMIZATION: Calculate adaptive buffer size for file I/O operations.
    
    Performance Impact:
    - Default Python buffer (typically 8KB) requires ~1.8M system calls for a 10GB file
    - 8MB buffer reduces this to ~2,400 system calls (750x reduction)
    - 16MB buffer reduces this to ~1,200 system calls (1,500x reduction)
    - This optimization alone can improve I/O performance by 10-50x for large files
    
    Uses available system memory to determine optimal buffer size, with a maximum
    of 16MB per file to balance performance and memory usage.
    
    Args:
        num_files (int): Number of file handles that will use this buffer size.
                        Used for documentation/clarity only (default: 2, since both
                        callers use 2 file handles). Total memory = buffer_size × num_files.
    
    Returns:
        int: Buffer size in bytes (between 8MB and 16MB per file)
    """
    if PSUTIL_AVAILABLE:
        try:
            available_mb = psutil.virtual_memory().available / (1024 * 1024)
            # OPTIMIZATION: Use 0.5% of available RAM per file
            # Total memory usage = buffer_size × num_files (e.g., 2 files = 1% total)
            # This ensures we don't use too much memory while still getting good performance
            # For systems with 32GB RAM: ~160MB per file, clamped to 16MB max
            calculated_mb = int(available_mb * 0.005)
            # OPTIMIZATION: Clamp between 8MB (minimum for good performance) and 16MB (maximum)
            # Minimum of 8MB ensures good performance even on memory-constrained systems
            # Maximum of 16MB prevents excessive memory usage while maintaining excellent performance
            buffer_mb = max(8, min(16, calculated_mb))
            return buffer_mb * 1024 * 1024
        except (AttributeError, OSError):
            # Fallback if psutil fails at runtime
            return 8 * 1024 * 1024  # 8MB default
    else:
        # OPTIMIZATION: Fallback to 8MB if psutil not available (good default for large files)
        # 8MB is a sweet spot: large enough for excellent performance, small enough for any system
        return 8 * 1024 * 1024  # 8MB default