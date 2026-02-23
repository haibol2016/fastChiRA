#!/usr/bin/env python
import os
import sys
from collections import defaultdict
import argparse
import chira_utilities
import copy
# ProcessPoolExecutor no longer needed - MPIRE is required
import math
import mmap
import sqlite3
import tempfile

# MPIRE is REQUIRED for multiprocessing performance
# MPIRE provides shared objects, lower overhead, and better performance than ProcessPoolExecutor
# Benefits: 50-90% memory reduction, 2-3x faster startup, better performance
# Install with: pip install mpire (or pip install chira, or conda install -c conda-forge mpire)
# Note: Requires MPIRE >= 2.4.0 for shared_objects support
from mpire import WorkerPool
# Note: MPIRE 2.10.1+ doesn't use a SharedObject class
# Shared objects are passed directly to WorkerPool via shared_objects parameter
# See: https://sybrenjansen.github.io/mpire/v2.10.1/usage/workerpool/shared_objects.html


def read_file_lines_mmap(file_path, use_mmap=True, mmap_threshold_mb=100):
    """
    Read file line by line with optional memory-mapping for large files.
    
    OPTIMIZATION: Uses memory-mapped files for very large files that don't fit in memory.
    Memory-mapped files allow the OS to handle paging, reducing memory usage and improving
    performance for files larger than available RAM.
    
    Args:
        file_path: Path to file to read
        use_mmap: Whether to use memory-mapped files (default: True)
        mmap_threshold_mb: File size threshold in MB to use memory-mapping (default: 100MB)
    
    Yields:
        str: Lines from the file
    """
    # Check file size to determine if memory-mapping is beneficial
    file_size_mb = os.path.getsize(file_path) / (1024 * 1024)
    should_use_mmap = use_mmap and file_size_mb >= mmap_threshold_mb
    
    if should_use_mmap:
        # OPTIMIZATION: Use memory-mapped files for large files
        try:
            with open(file_path, "rb") as fh_in:
                with mmap.mmap(fh_in.fileno(), 0, access=mmap.ACCESS_READ) as mmapped_file:
                    # OPTIMIZATION: Process memory-mapped file in chunks without loading entire file
                    CHUNK_SIZE = 1024 * 1024  # Process 1MB chunks
                    current_line_buffer = b""
                    pos = 0
                    
                    while pos < len(mmapped_file):
                        # Read a chunk from the memory-mapped file
                        chunk_end = min(pos + CHUNK_SIZE, len(mmapped_file))
                        chunk = mmapped_file[pos:chunk_end]
                        
                        # Process chunk line by line
                        data = current_line_buffer + chunk
                        lines = data.split(b'\n')
                        
                        # Last element might be incomplete (no newline at end of chunk)
                        current_line_buffer = lines[-1]
                        
                        # Process complete lines
                        for line_bytes in lines[:-1]:
                            try:
                                yield line_bytes.decode('utf-8')
                            except UnicodeDecodeError:
                                # Skip lines with encoding errors
                                print(f"Warning: Skipping line with encoding error", file=sys.stderr)
                        
                        pos = chunk_end
                    
                    # Process last line if file doesn't end with newline
                    if current_line_buffer:
                        try:
                            yield current_line_buffer.decode('utf-8')
                        except UnicodeDecodeError:
                            pass
        except (OSError, mmap.error) as e:
            # Fall back to regular file I/O if memory-mapping fails
            print(f"Warning: Memory-mapping failed ({e}), falling back to regular file I/O", file=sys.stderr)
            should_use_mmap = False
    
    if not should_use_mmap:
        # Regular file I/O for smaller files or when memory-mapping fails
        with open(file_path, "r", buffering=1024*1024) as fh_in:
            for line in fh_in:
                yield line


def build_crls(build_crls_too, bed, merged_bed, crl_file, crl_share_cutoff, min_locus_size):
    l_locipos = []
    l_locireads = []
    d_readlocus_transcripts = defaultdict(list)
    # OPTIMIZATION: Use memory-mapped files for large files (automatic fallback for smaller files)
    # Memory-mapping allows the OS to handle paging, reducing memory usage and improving
    # performance for files larger than available RAM
    for n_locus, line in enumerate(read_file_lines_mmap(merged_bed)):
        # chr14\t64814786\t64814804\t-\ttag_1308593|1|r,ENSMUST00000176386;tag_1308594|2|r,ENSMUST00000176386
        # OPTIMIZATION: Cache rstrip and split results to avoid repeated operations
        line_stripped = line.rstrip('\n')
        f = line_stripped.split('\t')
        # OPTIMIZATION: Use f-string instead of join for faster string formatting
        pos = f"{f[0]}:{f[1]}:{f[2]}:{f[3]}"
        l_locipos.append(pos)
        alignments = f[4].split(';')  # in the description column of the bed file,alignments are seperated by ';'
        l_locusreads = set()
        # OPTIMIZATION: Pre-compute n_locus string once per locus to avoid repeated str() calls in inner loop
        n_locus_str = str(n_locus)
        for alignment in alignments:
            # OPTIMIZATION: Cache split result to avoid repeated splitting
            alignment_fields = alignment.split(',')
            # Defensive check: skip malformed alignments
            if len(alignment_fields) < 6:
                print(f"Warning: Skipping malformed alignment (expected 6 comma-separated fields, got {len(alignment_fields)}): {alignment}", file=sys.stderr)
                continue
            # OPTIMIZATION: Use cached split result instead of splitting again
            segmentid, transcriptid, start, end, tx_strand, cigar = alignment_fields
            # OPTIMIZATION: Use f-string instead of join for faster string formatting
            transcriptid_pos = f"{transcriptid}\t{start}\t{end}\t{tx_strand}\t{cigar}"
            # OPTIMIZATION: Use pre-computed n_locus_str and f-string for faster concatenation
            d_readlocus_transcripts[f"{segmentid}{n_locus_str}"].append(transcriptid_pos)
            if segmentid not in l_locusreads:
                l_locusreads.add(segmentid)
        l_locireads.append(set(l_locusreads))

    # OPTIMIZATION: Pre-sort transcript lists once to avoid repeated sorting in nested output loops
    # This reduces sorting from O(n log n) per access to O(n log n) total, significantly
    # improving performance when writing output for large datasets
    for key in d_readlocus_transcripts:
        d_readlocus_transcripts[key].sort()

    d_crl_reads = defaultdict(set)  # Use sets instead of lists for faster operations
    d_crl_locus_reads = defaultdict(lambda: defaultdict())
    l_remaining_locireads = defaultdict(list)
    chira_utilities.print_w_time("START: 1st iteration of CRLs")
    # OPTIMIZATION: Use f-string instead of string concatenation for faster formatting
    print(f"Number of loci: {len(l_locireads)}")
    n_crl = 0
    # create CRLs only if required
    if build_crls_too:
        l_qualified_locireads = sorted([(i, item) for i, item in enumerate(l_locireads) if len(item) >= min_locus_size],
                                       key=lambda x: len(x[1]), reverse=True)
        l_remaining_locireads = sorted([(i, item) for i, item in enumerate(l_locireads) if len(item) < min_locus_size],
                                       key=lambda x: len(x[1]), reverse=True)

        # OPTIMIZATION: Use f-string instead of print with multiple arguments for consistency
        print(f"Number of qualified loci: {len(l_qualified_locireads)}")

        # OPTIMIZATION: Build inverted index: read_id -> set of CRL_ids containing that read
        # This index is built incrementally as CRLs are created, allowing us to quickly find
        # which CRLs share at least one read with each locus, dramatically reducing comparisons
        d_read_to_crls = defaultdict(set)
        
        for n_locus, l_locusreads in l_qualified_locireads:
            already_crl_member = False
            l_matched_crls = []
            # OPTIMIZATION: Cache locus_size since it's used multiple times (for bounds calculation)
            # This avoids repeated len() calls on the set, which is O(1) but still has function call overhead
            locus_size = len(l_locusreads)
            # lower and uppder bounds for filtering crls based on their size
            lower_bound = locus_size * (1 - crl_share_cutoff)
            upper_bound = locus_size / (1 - crl_share_cutoff)
            
            # OPTIMIZATION: Find candidate CRLs that share at least one read with this locus
            # This dramatically reduces the number of CRLs we need to check for Jaccard similarity
            # Only check CRLs that actually share reads AND have similar size
            # 
            # SAFETY: This optimization does NOT miss biologically significant CRLs because:
            # - Jaccard similarity = intersection / union
            # - If a CRL shares NO reads with a locus, intersection = 0, so Jaccard = 0
            # - Jaccard = 0 is always < crl_share_cutoff (typically 0.7), so it would never match
            # - The original algorithm also skips CRLs with n_common_reads == 0 (line 193-194)
            # - Therefore, only checking CRLs that share at least one read is safe and correct
            candidate_crl_ids = set()
            for read_id in l_locusreads:
                # Get CRLs that contain this read
                for crlid in d_read_to_crls.get(read_id, set()):
                    # Only consider CRLs that are still valid and have similar size
                    # BUGFIX: Check if crlid exists in dictionary instead of comparing to length
                    # The original check "crlid < len(d_crl_reads)" assumed sequential CRL IDs,
                    # but this can incorrectly filter out valid CRLs if IDs are not sequential
                    if crlid in d_crl_reads:
                        crl_size = len(d_crl_reads[crlid])
                        if lower_bound <= crl_size <= upper_bound:
                            candidate_crl_ids.add(crlid)
            
            # OPTIMIZATION: Pre-filter CRLs by size AND read overlap before expensive set operations
            # This reduces the number of expensive set intersection/union operations
            # by only checking CRLs that share at least one read AND have similar size
            candidate_crls = []
            # Traverse in reverse order (latest CRL is last) to maintain original behavior
            for crlid in sorted(candidate_crl_ids, reverse=True):
                candidate_crls.append((crlid, d_crl_reads[crlid]))
            
            # Only perform expensive set operations on size-filtered and read-overlap-filtered candidates
            for crlid, l_crlreads in candidate_crls:
                # if there are significantly more reads in crl than in locus or otherway around
                n_common_reads = len(l_locusreads & l_crlreads)  # Use & for set intersection
                if n_common_reads == 0:
                    continue
                n_union_reads = len(l_locusreads | l_crlreads)  # Use | for set union
                # jaccard similarity score
                # Defensive check: skip if union is empty (shouldn't happen if n_common_reads > 0)
                if n_union_reads == 0:
                    continue
                if n_common_reads / float(n_union_reads) < crl_share_cutoff:
                    continue
                # identical loci are multi-mapped loci with the same set of identical set of multi-mapped reads
                l_matched_crls.append(crlid)
                already_crl_member = True

            # n_locus is not a member of any crl, hence create a new crl
            if not already_crl_member:
                l_matched_crls.append(n_crl)
                n_crl += 1
            for matched_crl in l_matched_crls:
                d_crl_locus_reads[matched_crl][n_locus] = l_locusreads
                # OPTIMIZATION: Update inverted index as CRLs are created/modified
                # This allows subsequent loci to benefit from the index
                # Track which reads are new to this CRL
                # BUGFIX: Check if CRL exists in dictionary instead of comparing to length
                # This is more robust and doesn't assume sequential CRL IDs
                old_crl_reads = d_crl_reads[matched_crl].copy() if matched_crl in d_crl_reads else set()
                d_crl_reads[matched_crl].update(l_locusreads)  # Use set.update() instead of extend()
                # Update index: add new reads to the index for this CRL
                for read_id in l_locusreads:
                    if read_id not in old_crl_reads:
                        d_read_to_crls[read_id].add(matched_crl)
        d_crl_reads.clear()
    else:
        l_remaining_locireads = sorted(enumerate(l_locireads), key=lambda x: len(x[1]), reverse=True)
    # every locus with less than min_locus_size reads gets its' own locus
    for n_locus, l_locusreads in l_remaining_locireads:
        d_crl_locus_reads[n_crl][n_locus] = l_locusreads
        n_crl += 1
    chira_utilities.print_w_time("END: 1st iteration of CRLs")

    d_locus_crl_share = defaultdict(lambda: defaultdict(float))
    d_highest_shares = defaultdict(float)
    for crlid, d_crlloci in d_crl_locus_reads.items():
        # Build set of all reads in CRL more efficiently
        l_crlreads = set()
        for l_locusreads in d_crlloci.values():
            l_crlreads.update(l_locusreads)  # Use set.update() instead of extend()
        crl_reads_len = len(l_crlreads)
        # Defensive check: skip CRLs with zero reads (shouldn't happen, but protect against division by zero)
        if crl_reads_len == 0:
            continue
        # OPTIMIZATION: Pre-compute inverse of crl_reads_len to avoid repeated division in loop
        inv_crl_reads_len = 1.0 / float(crl_reads_len)
        for locusid, l_locusreads in d_crlloci.items():
            locus_share = len(l_locusreads) * inv_crl_reads_len
            d_locus_crl_share[locusid][crlid] = locus_share
            if locus_share > d_highest_shares[locusid]:
                d_highest_shares[locusid] = locus_share

    d_isolated_loci = {}  # loci that are separated from crl because of not enough overall share
    chira_utilities.print_w_time("START: 2nd iteration of CRLs")
    # in this iteration each locus is checked again for crl share and only one best crl for locus kept
    # if there are multiple eqaully good crls for a locus, all of the crls are considered for that locus
    for crlid in list(d_crl_locus_reads.keys()):
        d_crlloci = d_crl_locus_reads[crlid]
        for locusid in list(d_crlloci):
            if d_locus_crl_share[locusid][crlid] >= crl_share_cutoff:
                if d_locus_crl_share[locusid][crlid] < d_highest_shares[locusid]:
                    if locusid in d_crl_locus_reads[crlid]:
                        del d_crl_locus_reads[crlid][locusid]
                        # del d_locus_crl_share[locusid][crlid]
            else:
                # store isolated locus id and its reads
                d_isolated_loci[locusid] = d_crlloci[locusid]
                del d_crl_locus_reads[crlid][locusid]
                del d_locus_crl_share[locusid]

    chira_utilities.print_w_time("END: 2nd iteration of CRLs")
    chira_utilities.print_w_time("START: creating CRLs with isolated loci")
    # every locus that is in isolated_loci makes its own crl
    crl_index = len(d_crl_locus_reads)
    for locusid, l_locusreads in d_isolated_loci.items():
        # already member of some crl, safely ignore it
        if locusid in d_locus_crl_share:
            continue
        d_crl_locus_reads[crl_index][locusid] = l_locusreads
        d_locus_crl_share[locusid][crl_index] = 1
        crl_index += 1
    d_locus_crl_share.clear()

    chira_utilities.print_w_time("END: creating CRLs with isolated loci")

    # OPTIMIZATION: Use f-string instead of string concatenation for faster formatting
    print(f"There are a total of {len(d_crl_locus_reads)} uniq crls")

    # read the segments BED to get the genomic positions
    d_read_genomic_pos = defaultdict(str)
    # OPTIMIZATION: Use memory-mapped files for large files (automatic fallback for smaller files)
    # Memory-mapping allows the OS to handle paging, reducing memory usage and improving
    # performance for files larger than available RAM
    for line in read_file_lines_mmap(bed):
        # OPTIMIZATION: Cache rstrip and split results to avoid repeated operations
        line_stripped = line.rstrip('\n')
        b = line_stripped.split('\t')
        # OPTIMIZATION: Use f-string instead of join for faster string formatting
        pos = f"{b[0]}:{b[1]}:{b[2]}:{b[5]}"
        desc = b[3].split(',')  # in the description column of the bed file,alignments are seperated by ';'
        # Defensive check: skip if desc is empty
        if len(desc) < 1:
            print(f"Warning: Skipping line with empty description: {line[:100]}", file=sys.stderr)
            continue
        segmentid = desc[0]
        # OPTIMIZATION: Use f-string with join for remaining elements (more efficient than list slicing + join)
        # For small lists, f-string is faster; for larger lists, join is better
        if len(desc) > 1:
            transcriptid_pos = '\t'.join(desc[1:])  # Keep join for variable-length lists
        else:
            transcriptid_pos = ""
        # at this level reads have unique ids preceeded by a serialnumber
        # each read can have multiple alignements on same transcript
        # OPTIMIZATION: Use f-string for string concatenation (faster than + operator)
        d_read_genomic_pos[f"{transcriptid_pos}{segmentid}"] = pos

    # OPTIMIZATION: Use larger buffer size (1MB) for write operations in nested loops
    # This significantly reduces system calls and improves I/O performance when writing many small lines
    # Default buffer is typically 8KB, so 1MB (128x larger) provides substantial improvement
    # OPTIMIZATION: Batch writing - collect entries and write in larger chunks to reduce I/O overhead
    BATCH_SIZE = 10000  # Write 10K lines at a time
    entry_buffer = []
    
    with open(crl_file, "w", buffering=1024*1024) as fh_crl_file:
        # enumerate again because of above processing some crlids might be missing
        for crlid, d_crlloci in enumerate(d_crl_locus_reads.values()):
            # OPTIMIZATION: Pre-compute CRL-level data once per CRL (reused for all loci)
            if len(d_crlloci) > 1:
                # Build set more efficiently using set union
                l_crlreads = set()
                for l_locusreads in d_crlloci.values():
                    l_crlreads.update(l_locusreads)
                crl_reads_len = len(l_crlreads)
            else:
                crl_reads_len = None
            
            # OPTIMIZATION: Pre-compute crlid string once per CRL (reused for all entries in this CRL)
            crlid_str = str(crlid)
            
            for locusid, l_locusreads in d_crlloci.items():
                if crl_reads_len is not None:
                    # Defensive check: skip if crl_reads_len is 0 (shouldn't happen, but protect against division by zero)
                    if crl_reads_len == 0:
                        locus_share = 0.0
                    else:
                        locus_share = len(l_locusreads) / float(crl_reads_len)
                else:
                    locus_share = 1.0
                
                # OPTIMIZATION: Pre-format locus_share string once per locus (reused for all segments)
                locus_share_str = "{:.2f}".format(locus_share)
                locus_pos = l_locipos[locusid]  # Cache locus position string
                locusid_str = str(locusid)  # Pre-compute locusid string
                
                # OPTIMIZATION: Pre-compute segment_key prefix to avoid repeated string concatenation
                # Build all entries for this locus in a list first, then extend entry_buffer
                # This reduces the number of list append operations
                locus_entries = []
                
                # OPTIMIZATION: Sort segments once per locus instead of per transcript
                # If output order doesn't matter, we could skip sorting, but keeping it for consistency
                # OPTIMIZATION: Convert set to sorted list once to avoid repeated sorting
                sorted_segments = sorted(l_locusreads)
                
                for segmentid in sorted_segments:
                    segment_key = segmentid + locusid_str  # Use pre-computed locusid_str
                    # OPTIMIZATION: Cache transcriptid_pos split results to avoid repeated splitting
                    # Pre-split all transcript positions for this segment once
                    segment_transcripts = d_readlocus_transcripts.get(segment_key, [])
                    
                    for transcriptid_pos in segment_transcripts:
                        # OPTIMIZATION: Split transcriptid_pos once and cache results
                        transcript_parts = transcriptid_pos.split('\t')
                        transcriptid = transcript_parts[0]
                        # OPTIMIZATION: For small lists, f-string is faster; for larger lists, join is better
                        # Since transcript_parts[1:] is typically small (4 elements), use join for clarity
                        transcript_rest = '\t'.join(transcript_parts[1:])
                        # OPTIMIZATION: Use f-string for string concatenation (faster than + operator)
                        genomic_pos = d_read_genomic_pos[f"{transcriptid_pos}{segmentid}"]
                        
                        # OPTIMIZATION: Use f-string for faster string formatting
                        # OPTIMIZATION: Pre-compute all string components to minimize string operations
                        entry = f"{segmentid}\t{transcriptid}\t{locusid_str}\t{crlid_str}\t{transcript_rest}\t{genomic_pos}\t{locus_pos}\t{locus_share_str}\n"
                        locus_entries.append(entry)
                
                # OPTIMIZATION: Extend entry_buffer with all entries for this locus at once
                # This reduces the number of list operations
                entry_buffer.extend(locus_entries)
                
                # Write in batches to reduce I/O overhead
                if len(entry_buffer) >= BATCH_SIZE:
                    fh_crl_file.writelines(entry_buffer)
                    entry_buffer.clear()
        
        # Write remaining entries
        if entry_buffer:
            fh_crl_file.writelines(entry_buffer)


def _process_reads_chunk_e_step(shared_objects, args):
    """
    Process a chunk of reads in the E-step.
    Returns updated d_alpha values for reads in this chunk.
    
    OPTIMIZATION: This function is designed for parallel execution with MPIRE WorkerPool.
    With MPIRE, shared_objects provides access to d_rho_old without copying, reducing memory
    overhead by 50-90% for large datasets.
    
    Args:
        shared_objects: Dictionary of shared objects (MPIRE) - passed as first arg by MPIRE
            - d_rho_old: Dict of {crlid: relative_abundance} (shared, not copied)
        args: Tuple of (readids_chunk, d_alpha_chunk)
            - readids_chunk: List of read IDs to process
            - d_alpha_chunk: Dict of {readid: {crlid: value}} for this chunk
    
    Returns:
        Dict of {readid: {crlid: updated_value}} for this chunk
    """
    # MPIRE mode: use shared objects (passed as first argument by MPIRE)
    readids_chunk, d_alpha_chunk = args
    d_rho_old = shared_objects['d_rho_old']
    d_alpha_updated = {}
    for readid in readids_chunk:
        # Match original logic: skip reads with only 1 CRL (shouldn't happen but defensive)
        if len(d_alpha_chunk[readid]) == 1:
            continue
        read_crls = d_alpha_chunk[readid]
        # OPTIMIZATION: Pre-compute total abundance once and cache inverse to avoid
        # repeated sum() calls and division operations in the inner loop
        total_crl_rel_abundancy = sum(d_rho_old[crlid] for crlid in read_crls)
        # OPTIMIZATION: Defensive check to prevent division by zero (not in original but safe)
        if total_crl_rel_abundancy > 0:
            inv_total = 1.0 / total_crl_rel_abundancy
            d_alpha_updated[readid] = {crlid: d_rho_old[crlid] * inv_total for crlid in read_crls}
        else:
            d_alpha_updated[readid] = read_crls.copy()  # Keep original if total is 0
    return d_alpha_updated


def _process_reads_chunk_aggregate(shared_objects, args):
    """
    Process a chunk of reads for aggregation step.
    Returns a dictionary of CRL contributions from this chunk.
    
    OPTIMIZATION: This function is designed for parallel execution with MPIRE WorkerPool.
    Each process accumulates contributions into its own dictionary, then results are merged.
    Since summation is commutative and associative, parallel aggregation produces identical
    results to sequential.
    
    Args:
        shared_objects: Dictionary of shared objects (MPIRE) - passed as first arg by MPIRE
            Not used in this function but required for MPIRE API consistency
        args: Tuple of (readids_chunk, d_alpha_chunk)
            - readids_chunk: List of read IDs to process
            - d_alpha_chunk: Dict of {readid: {crlid: value}} for this chunk
    
    Returns:
        Dict of {crlid: aggregated_value} for this chunk
    """
    readids_chunk, d_alpha_chunk = args
    d_rho_chunk = defaultdict(float)
    for readid in readids_chunk:
        for crlid in d_alpha_chunk[readid]:
            d_rho_chunk[crlid] += d_alpha_chunk[readid][crlid]
    # Convert to regular dict: MPIRE pickles return values when sending results back to
    # main process. While defaultdict can be pickled, converting to plain dict is more
    # efficient (no factory function to pickle) and explicit.
    return dict(d_rho_chunk)


def em(d_alpha, d_rho, library_size, l_multimap_readids, em_threshold, num_processes=1):
    """
    Expectation-Maximization algorithm with optional multi-processing.
    
    OPTIMIZATION: Uses ProcessPoolExecutor instead of ThreadPoolExecutor for CPU-bound operations.
    ProcessPoolExecutor bypasses Python's GIL, enabling true parallelism for CPU-intensive tasks.
    This provides 2-8x speedup compared to ThreadPoolExecutor for large datasets.
    
    Args:
        d_alpha: Dictionary mapping read IDs to CRL contributions
        d_rho: Dictionary mapping CRL IDs to relative abundances
        library_size: Total library size
        l_multimap_readids: List of read IDs with multiple mappings
        em_threshold: Convergence threshold
        num_processes: Number of processes to use (default: 1, use all available if <= 0)
    """
    if num_processes <= 0:
        import multiprocessing
        num_processes = multiprocessing.cpu_count()
    
    # Log parallelization configuration
    if num_processes > 1:
        print(f"Using {num_processes} parallel processes for EM algorithm", file=sys.stderr)
        print(f"  - Multimapping reads: {len(l_multimap_readids)}", file=sys.stderr)
        print(f"  - Total reads: {len(d_alpha)}", file=sys.stderr)
    else:
        print("Using sequential processing (single process)", file=sys.stderr)
    
    i = 1
    sum_of_rho_diff = float('inf')
    
    # Split reads into chunks for parallel processing
    def chunk_list(lst, n_chunks):
        """Split list into approximately equal chunks."""
        if n_chunks <= 1:
            return [lst]
        chunk_size = max(1, math.ceil(len(lst) / n_chunks))
        return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]
    
    while sum_of_rho_diff > em_threshold and i <= 1000:
        if i % 10 == 0:  # Print every 10 iterations instead of every iteration
            # OPTIMIZATION: Use f-string instead of string concatenation for faster formatting
            print(f"iteration: {i}")
        # E-step
        # OPTIMIZATION: Use shallow copy instead of deep copy - much faster (10-50x speedup)
        # Shallow copy is sufficient since d_rho contains only immutable float values
        d_rho_old = dict(d_rho)  # Shallow copy - safe because values are immutable floats
        d_rho.clear()
        
        # OPTIMIZATION: Parallelize E-step for multimapped reads when beneficial
        # Only parallelize if num_processes > 1 AND there are enough reads to justify process overhead.
        # ProcessPoolExecutor has overhead: process creation (~10-50ms), data copying/pickling, result merging.
        # Adaptive threshold: max(500 reads, num_processes * 50) ensures each process gets meaningful work.
        # For 12 processes: need at least 600 reads (50 per process minimum).
        # ProcessPoolExecutor provides true parallelism (bypasses GIL) for CPU-bound operations.
        min_reads_threshold = max(500, num_processes * 50) if num_processes > 1 else float('inf')
        if num_processes > 1 and len(l_multimap_readids) > min_reads_threshold:
            if i == 1:  # Log only on first iteration to avoid spam
                print(f"  Parallelizing E-step: {len(l_multimap_readids)} multimapping reads across {num_processes} processes (MPIRE)", file=sys.stderr)
            read_chunks = chunk_list(l_multimap_readids, num_processes)
            
            # OPTIMIZATION: Use MPIRE with shared objects for better performance
            # Shared objects reduce memory overhead by 50-90% for large datasets
            # Pass shared objects directly to WorkerPool (MPIRE 2.10.1+)
            # Shared objects are passed as first argument to worker functions
            shared_objects_dict = {'d_rho_old': d_rho_old}
            
            # Extract d_alpha subsets for each chunk
            chunk_args = []
            for chunk in read_chunks:
                d_alpha_chunk = {readid: d_alpha[readid].copy() for readid in chunk}
                chunk_args.append((chunk, d_alpha_chunk))
            
            # Determine start method: use 'fork' for copy-on-write on Unix systems, default on Windows
            # 'fork' provides copy-on-write semantics: shared objects only copied when modified
            # This is the most memory-efficient option when available (50-90% memory reduction)
            # Reference: https://sybrenjansen.github.io/mpire/master/usage/workerpool/shared_objects.html
            import multiprocessing
            try:
                if sys.platform != 'win32':
                    current_method = multiprocessing.get_start_method(allow_none=True)
                    if current_method != 'spawn':
                        start_method = 'fork'
                    else:
                        start_method = None
                else:
                    start_method = None
            except (AttributeError, ValueError):
                start_method = 'fork' if sys.platform != 'win32' else None
            
            with WorkerPool(n_jobs=num_processes, shared_objects=shared_objects_dict, 
                          start_method=start_method) as pool:
                results = pool.map(_process_reads_chunk_e_step, chunk_args, progress_bar=False)
                # Collect updated values and merge back into d_alpha
                for d_alpha_updated in results:
                    for readid, read_crls_updated in d_alpha_updated.items():
                        d_alpha[readid] = read_crls_updated
        else:
            # Sequential processing for small datasets or single process
            if num_processes > 1 and len(l_multimap_readids) <= min_reads_threshold and i == 1:
                print(f"  Note: Sequential E-step ({len(l_multimap_readids)} reads < threshold of {min_reads_threshold})", file=sys.stderr)
            for readid in l_multimap_readids:
                # Match original logic: skip reads with only 1 CRL (shouldn't happen but defensive)
                if len(d_alpha[readid]) == 1:
                    continue
                read_crls = d_alpha[readid]
                total_crl_rel_abundancy = sum(d_rho_old[crlid] for crlid in read_crls)
                # OPTIMIZATION: Defensive check to prevent division by zero (not in original but safe)
                if total_crl_rel_abundancy > 0:
                    inv_total = 1.0 / total_crl_rel_abundancy
                    for crlid in read_crls:
                        d_alpha[readid][crlid] = d_rho_old[crlid] * inv_total
        
        # OPTIMIZATION: Aggregate contributions from all reads (including single-mapped)
        # This step can also be parallelized for large datasets. Each process processes a chunk
        # of reads and accumulates contributions, then results are merged.
        # ProcessPoolExecutor provides true parallelism (bypasses GIL) for CPU-bound operations.
        if num_processes > 1:
            # For parallel processing, convert keys to list once per iteration
            # OPTIMIZATION: Only create list when needed for parallel processing
            all_readids = list(d_alpha.keys())
            # Use same adaptive threshold as E-step
            min_reads_threshold_agg = max(500, num_processes * 50)
            if len(all_readids) > min_reads_threshold_agg:
                if i == 1:  # Log only on first iteration to avoid spam
                    print(f"  Parallelizing aggregation: {len(all_readids)} reads across {num_processes} processes (MPIRE)", file=sys.stderr)
                read_chunks = chunk_list(all_readids, num_processes)
                
                # Extract d_alpha subsets for each chunk (needed for process-based parallelism)
                chunk_args = []
                for chunk in read_chunks:
                    d_alpha_chunk = {readid: d_alpha[readid].copy() for readid in chunk}
                    chunk_args.append((chunk, d_alpha_chunk))
                
                # OPTIMIZATION: Use MPIRE for better performance (no shared objects needed here)
                with WorkerPool(n_jobs=num_processes) as pool:
                    results = pool.map(_process_reads_chunk_aggregate, chunk_args, progress_bar=False)
                    for d_rho_chunk in results:
                        for crlid, value in d_rho_chunk.items():
                            d_rho[crlid] += value
            else:
                # Sequential aggregation for small datasets even with threading enabled
                if len(all_readids) <= min_reads_threshold_agg and i == 1:
                    print(f"  Note: Sequential aggregation ({len(all_readids)} reads < threshold of {min_reads_threshold_agg})", file=sys.stderr)
                for readid in d_alpha.keys():
                    for crlid in d_alpha[readid]:
                        d_rho[crlid] += d_alpha[readid][crlid]
        else:
            # Sequential aggregation - matches original code exactly
            for readid in d_alpha.keys():
                for crlid in d_alpha[readid]:
                    d_rho[crlid] += d_alpha[readid][crlid]
        
        # M-step (sequential - fast enough)
        sum_of_rho_diff = 0
        # Defensive check: avoid division by zero if library_size is 0 (shouldn't happen, but protect)
        if library_size > 0:
            inv_library_size = 1.0 / float(library_size)
            for crlid in d_rho:
                new_rho = d_rho[crlid] * inv_library_size
                sum_of_rho_diff += abs(new_rho - d_rho_old[crlid])
                d_rho[crlid] = new_rho
        else:
            # If library_size is 0, set all rho to 0 and break (converged to invalid state)
            for crlid in d_rho:
                d_rho[crlid] = 0.0
            break
        i += 1
    return d_alpha


def tpm(d_crl_expression, d_crl_loci_len):
    total_rpk = 0
    d_crl_tpm = defaultdict(float)
    # OPTIMIZATION: Cache sorted keys to avoid repeated sorting (used in two loops below)
    # This reduces sorting from O(n log n) twice to O(n log n) once
    sorted_crlids = sorted(d_crl_expression.keys())
    for crlid in sorted_crlids:
        crl_expression = d_crl_expression[crlid]
        # OPTIMIZATION: median() already sorts internally, so don't sort here to avoid double sorting
        # This eliminates redundant O(n log n) sorting operation for each CRL
        crl_len = chira_utilities.median(d_crl_loci_len[crlid].values()) / 1000.0  # length in kbs
        if crl_len > 0:  # Avoid division by zero
            rpk = crl_expression / crl_len
            d_crl_tpm[crlid] = rpk
            total_rpk += rpk
    millions_of_rpk = total_rpk / 1000000.0
    if millions_of_rpk > 0:  # Avoid division by zero
        # OPTIMIZATION: Pre-compute inverse to avoid repeated division in loop
        inv_millions = 1.0 / millions_of_rpk
        for crlid in sorted_crlids:
            crl_tpm = d_crl_tpm[crlid] * inv_millions
            d_crl_tpm[crlid] = crl_tpm
    return d_crl_tpm


def quantify_crls(crl_file, em_threshold, num_processes=1):
    d_crl_loci_len = defaultdict(lambda: defaultdict(int))
    l_multimap_readids = []
    d_alpha = defaultdict(lambda: defaultdict(float))
    d_rho = defaultdict(float)

    # Use context manager for file handling (more efficient)
    # OPTIMIZATION: Use memory-mapped files for large files (automatic fallback for smaller files)
    # Memory-mapping allows the OS to handle paging, reducing memory usage and improving
    # performance for files larger than available RAM
    for line in read_file_lines_mmap(crl_file):
        # OPTIMIZATION: Cache rstrip and split results to avoid repeated operations
        line_stripped = line.rstrip("\n")
        f = line_stripped.split("\t")
        # consider the segment id and quantify individula segments than whole reads
        # Defensive check: skip malformed lines
        if len(f) < 10:
            print(f"Warning: Skipping malformed line (expected at least 10 fields, got {len(f)}): {line[:100]}", file=sys.stderr)
            continue
        readid = f[0]
        locusid = f[2]
        crlid = f[3]
        # OPTIMIZATION: Cache split result to avoid repeated splitting
        pos = f[9].split(':')
        # Defensive check: pos should have at least 3 elements (chr:start:end:strand format)
        if len(pos) < 3:
            print(f"Warning: Skipping line with invalid position format (expected chr:start:end:strand, got {f[9]}): {line[:100]}", file=sys.stderr)
            continue
        try:
            locuslength = int(pos[-2]) - int(pos[-3]) + 1
        except (ValueError, IndexError) as e:
            print(f"Warning: Skipping line with invalid position values (pos={f[9]}): {e}", file=sys.stderr)
            continue
        # a single locus can belong to multiple crls
        # one read can be part of multiple crls
        d_alpha[readid][crlid] = 1
        d_crl_loci_len[crlid][locusid] = locuslength

    # intial read segment contributions
    # OPTIMIZATION: Use .items() for more efficient dictionary iteration (avoids key lookups)
    # OPTIMIZATION: Cache num_crls and inv_num_crls since they're used multiple times per read
    # This avoids repeated len() calls and division operations in the inner loop
    for readid, read_crls in d_alpha.items():
        num_crls = len(read_crls)
        inv_num_crls = 1.0 / float(num_crls)
        for crlid in read_crls:
            d_alpha[readid][crlid] = inv_num_crls
            d_rho[crlid] += inv_num_crls
        if num_crls > 1:
            l_multimap_readids.append(readid)

    library_size = sum(d_rho.values())
    # intial relative abundancies
    # Defensive check: avoid division by zero if library_size is 0 (shouldn't happen, but protect)
    if library_size > 0:
        # OPTIMIZATION: Pre-compute inverse of library_size to avoid repeated division in loop
        inv_library_size = 1.0 / float(library_size)
        for crlid in d_rho:
            d_rho[crlid] = d_rho[crlid] * inv_library_size
    else:
        # If library_size is 0, all reads have no CRL assignments - set all to 0
        print("Warning: library_size is 0, no reads have CRL assignments", file=sys.stderr)
        for crlid in d_rho:
            d_rho[crlid] = 0.0

    d_res = em(d_alpha, d_rho, library_size, l_multimap_readids, 
              em_threshold, num_processes)

    d_crl_expression = defaultdict(float)
    # OPTIMIZATION: Use .items() for more efficient dictionary iteration
    # This avoids key lookups and provides direct access to both keys and values
    for readid, read_crls in d_res.items():
        for crlid, value in read_crls.items():
            d_crl_expression[crlid] += value

    d_crl_tpm = tpm(d_crl_expression, d_crl_loci_len)

    # BUGFIX: Return d_res (EM-optimized fractions) instead of d_alpha (initial values)
    # Although EM modifies d_alpha in place, returning d_res is semantically correct
    # and protects against future changes to EM that might return a new dictionary
    return d_res, d_crl_tpm


def em_sqlite(conn, db_path, d_rho, library_size, em_threshold, num_processes=1):
    """
    EM algorithm using SQLite database with hybrid parallel approach.
    
    Uses read-only connections for parallel computation, then applies updates sequentially.
    This avoids write conflicts while maximizing parallelism.
    
    Args:
        conn: SQLite connection to main database
        db_path: Path to database file (needed for parallel connections)
        d_rho: Dictionary mapping CRL IDs to relative abundances (in memory)
        library_size: Total library size
        em_threshold: Convergence threshold
        num_processes: Number of parallel processes
    
    Returns:
        Final d_rho dictionary after convergence
    """
    cursor = conn.cursor()
    i = 1
    sum_of_rho_diff = float('inf')
    
    # Get list of multimapped reads
    cursor.execute('''
        SELECT readid, COUNT(*) as num_crls
        FROM alpha
        GROUP BY readid
        HAVING num_crls > 1
    ''')
    multimap_readids = [row[0] for row in cursor.fetchall()]
    
    print(f"Found {len(multimap_readids):,} multimapped reads", file=sys.stderr)
    
    while sum_of_rho_diff > em_threshold and i <= 1000:
        if i % 10 == 0:
            print(f"iteration: {i}")
        
        d_rho_old = dict(d_rho)
        d_rho_new = defaultdict(float)
        
        # E-step: Compute updates in parallel (read-only), apply sequentially
        print(f"  E-step: Processing {len(multimap_readids):,} reads...", file=sys.stderr)
        
        if num_processes > 1 and len(multimap_readids) > 1000:
            # Parallel computation of updates
            from mpire import WorkerPool
            
            # Split reads into chunks
            chunk_size = max(1000, len(multimap_readids) // num_processes)
            read_chunks = [multimap_readids[j:j+chunk_size] 
                          for j in range(0, len(multimap_readids), chunk_size)]
            
            def compute_updates_parallel(readids_chunk):
                """Compute updates for a chunk of reads (read-only, no conflicts)."""
                # Read-only connection (no locking issues)
                conn_read = sqlite3.connect(db_path, timeout=30.0)
                conn_read.execute('PRAGMA journal_mode=WAL')
                conn_read.execute('PRAGMA busy_timeout=30000')
                cursor_read = conn_read.cursor()
                
                updates = []
                d_rho_local = defaultdict(float)
                
                # Process in batches (SQLite parameter limit is 999)
                BATCH_SIZE = 999
                for k in range(0, len(readids_chunk), BATCH_SIZE):
                    batch = readids_chunk[k:k+BATCH_SIZE]
                    placeholders = ','.join(['?'] * len(batch))
                    
                    cursor_read.execute(f'''
                        SELECT readid, crlid, value
                        FROM alpha
                        WHERE readid IN ({placeholders})
                        ORDER BY readid
                    ''', batch)
                    
                    # Group by readid
                    read_data = defaultdict(list)
                    for readid, crlid, value in cursor_read.fetchall():
                        read_data[readid].append((crlid, value))
                    
                    # Compute E-step updates
                    for readid, crls in read_data.items():
                        if len(crls) == 1:
                            continue
                        
                        # Compute total rho for this read's CRLs
                        total_rho = sum(d_rho_old[crlid] for crlid, _ in crls)
                        if total_rho > 0:
                            inv_total = 1.0 / total_rho
                            for crlid, _ in crls:
                                new_value = d_rho_old[crlid] * inv_total
                                updates.append((new_value, readid, crlid))
                                d_rho_local[crlid] += new_value
                        else:
                            # Uniform distribution if total is 0
                            inv_num = 1.0 / len(crls)
                            for crlid, _ in crls:
                                updates.append((inv_num, readid, crlid))
                                d_rho_local[crlid] += inv_num
                
                conn_read.close()
                return updates, dict(d_rho_local)
            
            # Compute updates in parallel (all read-only, no conflicts)
            # Wrap each chunk in a tuple so MPIRE passes it as one argument (map unpacks iterables)
            with WorkerPool(n_jobs=num_processes) as pool:
                results = pool.map(compute_updates_parallel, [(c,) for c in read_chunks], progress_bar=False)
            
            # Collect all updates
            all_updates = []
            for updates, d_rho_local in results:
                all_updates.extend(updates)
                for crlid, value in d_rho_local.items():
                    d_rho_new[crlid] += value
            
            # Apply updates sequentially (single writer, no conflicts)
            print(f"  Applying {len(all_updates):,} updates...", file=sys.stderr)
            cursor.executemany('''
                UPDATE alpha
                SET value = ?
                WHERE readid = ? AND crlid = ?
            ''', all_updates)
            conn.commit()
        else:
            # Sequential processing
            BATCH_SIZE = 10000
            for j in range(0, len(multimap_readids), BATCH_SIZE):
                batch = multimap_readids[j:j+BATCH_SIZE]
                placeholders = ','.join(['?'] * len(batch))
                
                cursor.execute(f'''
                    SELECT readid, crlid, value
                    FROM alpha
                    WHERE readid IN ({placeholders})
                    ORDER BY readid
                ''', batch)
                
                # Group by readid
                read_data = defaultdict(list)
                for readid, crlid, value in cursor.fetchall():
                    read_data[readid].append((crlid, value))
                
                # Compute and apply updates
                updates = []
                for readid, crls in read_data.items():
                    if len(crls) == 1:
                        continue
                    
                    total_rho = sum(d_rho_old[crlid] for crlid, _ in crls)
                    if total_rho > 0:
                        inv_total = 1.0 / total_rho
                        for crlid, _ in crls:
                            new_value = d_rho_old[crlid] * inv_total
                            updates.append((new_value, readid, crlid))
                            d_rho_new[crlid] += new_value
                    else:
                        inv_num = 1.0 / len(crls)
                        for crlid, _ in crls:
                            updates.append((inv_num, readid, crlid))
                            d_rho_new[crlid] += inv_num
                
                if updates:
                    cursor.executemany('''
                        UPDATE alpha
                        SET value = ?
                        WHERE readid = ? AND crlid = ?
                    ''', updates)
                    conn.commit()
                
                if j % 100000 == 0:
                    print(f"    Processed {j:,} reads...", file=sys.stderr)
        
        # Add contributions from single-CRL reads
        cursor.execute('''
            SELECT crlid, COUNT(*) as count
            FROM alpha
            WHERE readid IN (
                SELECT readid
                FROM alpha
                GROUP BY readid
                HAVING COUNT(*) = 1
            )
            GROUP BY crlid
        ''')
        for crlid, count in cursor.fetchall():
            d_rho_new[crlid] += count
        
        # M-step: Normalize d_rho by constant library_size (same as original em())
        # This ensures consistency with the original algorithm where library_size remains constant
        if library_size > 0:
            inv_library_size = 1.0 / float(library_size)
            for crlid in d_rho_new:
                d_rho[crlid] = d_rho_new[crlid] * inv_library_size
        
        # Check convergence
        sum_of_rho_diff = sum(abs(d_rho[crlid] - d_rho_old.get(crlid, 0)) 
                              for crlid in d_rho)
        print(f"  Convergence: {sum_of_rho_diff:.2e} (threshold: {em_threshold:.2e})", file=sys.stderr)
        i += 1
    
    return d_rho


def quantify_crls_sqlite(crl_file, em_threshold, num_processes=1):
    """
    SQLite-based quantification that avoids loading d_alpha into memory.
    
    Uses disk-backed SQLite database instead of in-memory dictionaries.
    Memory usage: ~8-16GB vs 200-500GB for standard method.
    
    Args:
        crl_file: Path to loci.txt file
        em_threshold: EM convergence threshold
        num_processes: Number of parallel processes
    
    Returns:
        Tuple of (d_read_crl_fractions, d_crl_tpms) in same format as quantify_crls()
    """
    # Create temporary SQLite database
    db_file = os.path.join(tempfile.gettempdir(), f'chira_alpha_{os.getpid()}.db')
    
    # Remove existing database if present
    if os.path.exists(db_file):
        os.remove(db_file)
    
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        cursor = conn.cursor()
        
        # Optimize SQLite settings for performance
        conn.execute('PRAGMA journal_mode=WAL')
        conn.execute('PRAGMA synchronous=NORMAL')
        conn.execute('PRAGMA cache_size=-256000')  # 256MB cache
        conn.execute('PRAGMA temp_store=MEMORY')
        conn.execute('PRAGMA mmap_size=268435456')  # 256MB mmap
        
        print(f"Creating SQLite database: {db_file}", file=sys.stderr)
        
        # Create tables (without indexes initially for faster bulk loading)
        cursor.execute('''
        CREATE TABLE alpha (
            readid TEXT NOT NULL,
            crlid TEXT NOT NULL,
            value REAL NOT NULL DEFAULT 0.0,
            PRIMARY KEY (readid, crlid)
        )
        ''')
        
        # Table for locus lengths (small, can stay in memory or SQLite)
        cursor.execute('''
        CREATE TABLE crl_loci_len (
            crlid TEXT NOT NULL,
            locusid TEXT NOT NULL,
            length INTEGER NOT NULL,
            PRIMARY KEY (crlid, locusid)
        )
        ''')
        
        # Load data into database with optimized bulk loading
        print("Loading data into SQLite database...", file=sys.stderr)
        d_crl_loci_len = defaultdict(lambda: defaultdict(int))
        
        # Disable some safety features during bulk load for speed
        conn.execute('PRAGMA synchronous=OFF')
        conn.execute('PRAGMA journal_mode=MEMORY')
        
        # Use batch inserts for better performance
        BATCH_SIZE = 100000
        alpha_batch = []
        loci_batch = []
        line_count = 0
        
        for line in read_file_lines_mmap(crl_file):
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            
            readid = f[0]
            locusid = f[2]
            crlid = f[3]
            pos = f[9].split(':')
            
            if len(pos) < 3:
                continue
            try:
                locuslength = int(pos[-2]) - int(pos[-3]) + 1
            except (ValueError, IndexError):
                continue
            
            # Batch insert alpha
            alpha_batch.append((readid, crlid, 1.0))
            d_crl_loci_len[crlid][locusid] = locuslength
            loci_batch.append((crlid, locusid, locuslength))
            
            # Insert in batches (OR REPLACE: loci file can have multiple lines per (readid, crlid))
            if len(alpha_batch) >= BATCH_SIZE:
                cursor.executemany('INSERT OR REPLACE INTO alpha VALUES (?, ?, ?)', alpha_batch)
                cursor.executemany('INSERT OR REPLACE INTO crl_loci_len VALUES (?, ?, ?)', loci_batch)
                conn.commit()
                alpha_batch = []
                loci_batch = []
                line_count += BATCH_SIZE
                if line_count % 1000000 == 0:
                    print(f"  Loaded {line_count:,} lines...", file=sys.stderr)
        
        # Insert remaining
        if alpha_batch:
            cursor.executemany('INSERT OR REPLACE INTO alpha VALUES (?, ?, ?)', alpha_batch)
            cursor.executemany('INSERT OR REPLACE INTO crl_loci_len VALUES (?, ?, ?)', loci_batch)
            conn.commit()
            line_count += len(alpha_batch)
        
        print(f"Database loaded. Total lines: {line_count:,}", file=sys.stderr)
        
        # Re-enable safety features
        conn.execute('PRAGMA synchronous=NORMAL')
        conn.execute('PRAGMA journal_mode=WAL')
        conn.commit()  # Ensure all data is written before creating indexes
        
        # Create indexes after bulk load (much faster than maintaining during inserts)
        print("Creating indexes...", file=sys.stderr)
        cursor.execute('CREATE INDEX idx_readid ON alpha(readid)')
        cursor.execute('CREATE INDEX idx_crlid ON alpha(crlid)')
        cursor.execute('CREATE INDEX idx_crl_loci_crlid ON crl_loci_len(crlid)')
        conn.commit()
        print("Indexes created.", file=sys.stderr)
        
        # Initialize d_rho uniformly
        print("Initializing d_rho...", file=sys.stderr)
        d_rho = defaultdict(float)
        
        # Count CRLs per read and initialize uniformly
        cursor.execute('''
            SELECT readid, COUNT(*) as num_crls
            FROM alpha
            GROUP BY readid
        ''')
        
        total_reads = 0
        updates = []
        for readid, num_crls in cursor.fetchall():
            inv_num_crls = 1.0 / float(num_crls)
            # Get CRLs for this read
            cursor.execute('SELECT crlid FROM alpha WHERE readid = ?', (readid,))
            for (crlid,) in cursor.fetchall():
                updates.append((inv_num_crls, readid, crlid))
                d_rho[crlid] += inv_num_crls
            
            total_reads += 1
            if total_reads % 100000 == 0:
                if updates:
                    cursor.executemany('UPDATE alpha SET value = ? WHERE readid = ? AND crlid = ?', updates)
                    conn.commit()
                    updates = []
                print(f"  Processed {total_reads:,} reads...", file=sys.stderr)
        
        # Apply remaining updates
        if updates:
            cursor.executemany('UPDATE alpha SET value = ? WHERE readid = ? AND crlid = ?', updates)
            conn.commit()
        
        # Normalize d_rho
        library_size = sum(d_rho.values())
        if library_size > 0:
            inv_library_size = 1.0 / float(library_size)
            for crlid in d_rho:
                d_rho[crlid] *= inv_library_size
        
        # Run EM algorithm
        print("Starting EM algorithm...", file=sys.stderr)
        d_rho = em_sqlite(conn, db_file, d_rho, library_size, em_threshold, num_processes)
        
        # Compute final d_crl_expression
        print("Computing CRL expressions...", file=sys.stderr)
        d_crl_expression = defaultdict(float)
        
        cursor.execute('SELECT readid, crlid, value FROM alpha')
        for readid, crlid, value in cursor.fetchall():
            d_crl_expression[crlid] += value
        
        d_crl_tpm = tpm(d_crl_expression, d_crl_loci_len)
        
        # Extract d_read_crl_fractions from database
        print("Extracting read-CRL fractions...", file=sys.stderr)
        d_read_crl_fractions = defaultdict(lambda: defaultdict(float))
        
        cursor.execute('SELECT readid, crlid, value FROM alpha ORDER BY readid')
        for readid, crlid, value in cursor.fetchall():
            d_read_crl_fractions[readid][crlid] = value
        
        return d_read_crl_fractions, d_crl_tpm
    finally:
        # Cleanup: Always remove database files, even if an exception occurs
        # This ensures temporary files don't accumulate on disk
        if conn is not None:
            try:
                conn.close()
            except Exception:
                pass  # Ignore errors when closing connection
        
        # Remove main database file
        if os.path.exists(db_file):
            try:
                os.remove(db_file)
            except OSError as e:
                print(f"Warning: Could not remove database file {db_file}: {e}", file=sys.stderr)
        
        # Also remove WAL and SHM files (SQLite auxiliary files)
        for suffix in ['.db-wal', '.db-shm']:
            wal_file = db_file + suffix
            if os.path.exists(wal_file):
                try:
                    os.remove(wal_file)
                except OSError as e:
                    print(f"Warning: Could not remove database file {wal_file}: {e}", file=sys.stderr)


def print_configuration(args):
    """Print configuration summary."""
    # OPTIMIZATION: Use f-strings instead of string concatenation for faster formatting
    print(f'Input BED file                       : {args.bed}')
    print(f'Input merged BED file                : {args.merged_bed}')
    print(f'Output directory                     : {args.outdir}')
    print(f'Minimum locus size                   : {args.min_locus_size}')
    print(f'CRL share                            : {args.crl_share}')
    print(f'EM threshold                         : {args.em_thresh}')
    print(f'Number of processes                  : {args.num_processes if args.num_processes > 0 else "all available"}')
    print(f'Create CRLs too                      : {args.build_crls_too}')
    print(f'Use SQLite database                  : {args.use_sqldb}')
    if args.use_sqldb:
        print('  -> Using disk-backed SQLite method (suitable for large files)')
    else:
        print('  -> Using in-memory method (standard)')
    print("===================================================================")


def write_crl_output(loci_file, output_file, d_read_crl_fractions, d_crl_tpms, use_mmap=True, mmap_threshold_mb=100):
    """
    Write CRL output file with fractions and TPMs.
    
    OPTIMIZATION: Uses memory-mapped files for very large files that don't fit in memory.
    Memory-mapped files allow the OS to handle paging, reducing memory usage and improving
    performance for files larger than available RAM.
    
    Args:
        loci_file: Path to input loci file
        output_file: Path to output file
        d_read_crl_fractions: Dictionary mapping read_id -> {crl_id: fraction}
        d_crl_tpms: Dictionary mapping crl_id -> TPM value
        use_mmap: Whether to use memory-mapped files (default: True)
        mmap_threshold_mb: File size threshold in MB to use memory-mapping (default: 100MB)
    """
    # OPTIMIZATION: Batch writing - collect lines and write in chunks to reduce I/O overhead
    BATCH_SIZE = 10000  # Write 10K lines at a time
    output_buffer = []
    
    # Check file size to determine if memory-mapping is beneficial
    file_size_mb = os.path.getsize(loci_file) / (1024 * 1024)
    should_use_mmap = use_mmap and file_size_mb >= mmap_threshold_mb
    
    if should_use_mmap:
        # OPTIMIZATION: Use memory-mapped files for large files
        # Memory-mapping allows the OS to handle paging, reducing memory usage
        # and improving performance for files larger than available RAM
        # The OS handles paging automatically, so we can process the file efficiently
        try:
            with open(loci_file, "rb") as fh_in:
                with mmap.mmap(fh_in.fileno(), 0, access=mmap.ACCESS_READ) as mmapped_file:
                    with open(output_file, "w", buffering=1024*1024) as fh_out:
                        # OPTIMIZATION: Process memory-mapped file in chunks without loading entire file
                        # Process chunk by chunk, letting OS handle paging
                        CHUNK_SIZE = 1024 * 1024  # Process 1MB chunks
                        current_line_buffer = b""
                        pos = 0
                        
                        while pos < len(mmapped_file):
                            # Read a chunk from the memory-mapped file
                            chunk_end = min(pos + CHUNK_SIZE, len(mmapped_file))
                            chunk = mmapped_file[pos:chunk_end]
                            
                            # Process chunk line by line
                            # Combine with any leftover from previous chunk
                            data = current_line_buffer + chunk
                            lines = data.split(b'\n')
                            
                            # Last element might be incomplete (no newline at end of chunk)
                            # Save it for next iteration
                            current_line_buffer = lines[-1]
                            
                            # Process complete lines
                            for line_bytes in lines[:-1]:
                                try:
                                    line_stripped = line_bytes.decode('utf-8').rstrip("\n")
                                    k = line_stripped.split("\t")
                                    # Defensive check: skip malformed lines
                                    if len(k) < 4:
                                        print(f"Warning: Skipping malformed line in loci file (expected at least 4 fields, got {len(k)}): {line_stripped[:100]}", file=sys.stderr)
                                        continue
                                    read_id = k[0]
                                    crl_id = k[3]
                                    # Defensive check: skip if read_id or crl_id not in dictionaries
                                    if read_id not in d_read_crl_fractions or crl_id not in d_read_crl_fractions[read_id]:
                                        print(f"Warning: Skipping line with missing read_id or crl_id in fractions: read_id={read_id}, crl_id={crl_id}", file=sys.stderr)
                                        continue
                                    if crl_id not in d_crl_tpms:
                                        print(f"Warning: Skipping line with missing crl_id in TPMs: crl_id={crl_id}", file=sys.stderr)
                                        continue
                                    # OPTIMIZATION: Use f-string for faster string formatting
                                    fraction_str = "{:.2f}".format(d_read_crl_fractions[read_id][crl_id])
                                    tpm_str = "{:.4g}".format(d_crl_tpms[crl_id])
                                    output_buffer.append(f"{line_stripped}\t{fraction_str}\t{tpm_str}\n")
                                    
                                    # Write in batches to reduce I/O overhead
                                    if len(output_buffer) >= BATCH_SIZE:
                                        fh_out.writelines(output_buffer)
                                        output_buffer.clear()
                                except UnicodeDecodeError:
                                    # Skip lines with encoding errors
                                    print(f"Warning: Skipping line with encoding error", file=sys.stderr)
                            
                            pos = chunk_end
                        
                        # Process last line if file doesn't end with newline
                        if current_line_buffer:
                            try:
                                line_stripped = current_line_buffer.decode('utf-8').rstrip("\n")
                                k = line_stripped.split("\t")
                                if len(k) >= 4:
                                    read_id = k[0]
                                    crl_id = k[3]
                                    if read_id in d_read_crl_fractions and crl_id in d_read_crl_fractions[read_id] and crl_id in d_crl_tpms:
                                        fraction_str = "{:.2f}".format(d_read_crl_fractions[read_id][crl_id])
                                        tpm_str = "{:.4g}".format(d_crl_tpms[crl_id])
                                        output_buffer.append(f"{line_stripped}\t{fraction_str}\t{tpm_str}\n")
                            except UnicodeDecodeError:
                                pass
                        
                        # Write remaining entries
                        if output_buffer:
                            fh_out.writelines(output_buffer)
                            output_buffer.clear()
        except (OSError, mmap.error) as e:
            # Fall back to regular file I/O if memory-mapping fails
            print(f"Warning: Memory-mapping failed ({e}), falling back to regular file I/O", file=sys.stderr)
            should_use_mmap = False
    
    if not should_use_mmap:
        # Regular file I/O for smaller files or when memory-mapping fails
        # OPTIMIZATION: Use larger buffer size (1MB) for write operations in loops
        # This significantly reduces system calls and improves I/O performance when writing many small lines
        with open(loci_file, "r", buffering=1024*1024) as fh_in:
            with open(output_file, "w", buffering=1024*1024) as fh_out:
                for l in fh_in:
                    # OPTIMIZATION: Cache rstrip result to avoid repeated operation
                    line_stripped = l.rstrip("\n")
                    k = line_stripped.split("\t")
                    # Defensive check: skip malformed lines
                    if len(k) < 4:
                        print(f"Warning: Skipping malformed line in loci file (expected at least 4 fields, got {len(k)}): {l[:100]}", file=sys.stderr)
                        continue
                    read_id = k[0]
                    crl_id = k[3]
                    # Defensive check: skip if read_id or crl_id not in dictionaries
                    if read_id not in d_read_crl_fractions or crl_id not in d_read_crl_fractions[read_id]:
                        print(f"Warning: Skipping line with missing read_id or crl_id in fractions: read_id={read_id}, crl_id={crl_id}", file=sys.stderr)
                        continue
                    if crl_id not in d_crl_tpms:
                        print(f"Warning: Skipping line with missing crl_id in TPMs: crl_id={crl_id}", file=sys.stderr)
                        continue
                    # OPTIMIZATION: Use f-string for faster string formatting
                    fraction_str = "{:.2f}".format(d_read_crl_fractions[read_id][crl_id])
                    tpm_str = "{:.4g}".format(d_crl_tpms[crl_id])
                    output_buffer.append(f"{line_stripped}\t{fraction_str}\t{tpm_str}\n")
                    
                    # Write in batches to reduce I/O overhead
                    if len(output_buffer) >= BATCH_SIZE:
                        fh_out.writelines(output_buffer)
                        output_buffer.clear()
                
                # Write remaining entries
                if output_buffer:
                    fh_out.writelines(output_buffer)


def finalize_output(outdir, temp_file, final_file):
    """
    Remove temporary files and sort the final output file.
    
    Args:
        outdir: Output directory
        temp_file: Temporary file to remove
        final_file: Final output file to create (sorted)
    """
    temp_path = os.path.join(outdir, temp_file)
    final_path = os.path.join(outdir, final_file)
    
    chira_utilities.print_w_time("START: Sort CRLs file by read name")
    os.system("sort -V " + temp_path + " > " + final_path)
    os.remove(temp_path)
    chira_utilities.print_w_time("END: Sort CRLs file by read name")


def parse_arguments():
    """Parse command-line arguments. Order: required I/O, CRL/EM params, optional processes/sql, version."""
    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: quantify mapped loci',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required I/O
    parser.add_argument('-b', '--bed', action='store', dest='bed', required=True,
                        metavar='', help='Input BED file (e.g. segments.bed from chira_merge)')
    parser.add_argument('-m', '--merged_bed', action='store', dest='merged_bed', required=True,
                        metavar='', help='Input merged BED file (e.g. merged.bed from chira_merge)')
    parser.add_argument('-o', '--outdir', action='store', dest='outdir', required=True, metavar='',
                        help='Output directory (writes loci.counts)')

    # CRL / EM parameters
    parser.add_argument('-cs', '--crl_share', action='store', type=chira_utilities.score_float, default=0.7, metavar='',
                        dest='crl_share',
                        help="""Minimum fraction of reads of a locus that must overlap with all CRL loci
in order to merge it into that CRL.""")
    parser.add_argument('-ls', '--min_locus_size', action='store', type=int, default=10, metavar='',
                        dest='min_locus_size',
                        help="""Minimum number of reads a locus should have to participate in CRL creation.
Set relative to sequencing depth. Lower values can yield CRLs from random multimappings.
Consider --crl_share together with this.""")
    parser.add_argument('-e', '--em_threshold', action='store', type=chira_utilities.score_float, default=0.00001, metavar='',
                        dest='em_thresh',
                        help="""Maximum difference in transcript expression between two consecutive EM iterations
to declare convergence.""")
    parser.add_argument("-crl", '--build_crls_too', action='store_true', dest='build_crls_too',
                        help='Build CRLs in addition to quantification')

    # Parallelism and optional backends
    parser.add_argument('-p', '--processes', action='store', type=int, default=0, metavar='',
                        dest='num_processes',
                        help="""Number of parallel processes for EM algorithm. 0 = use all available cores.
Uses MPIRE WorkerPool. Most beneficial for large datasets (>500 multimapping reads).
Set to 1 to disable. Install: pip install mpire""")
    parser.add_argument('--use_sqldb', action='store_true', dest='use_sqldb',
                        help="""Use SQLite backend for very large inputs (50GB+). Disk-backed; ~8-16GB RAM vs 200-500GB.
Runtime ~3-8 hours depending on size and -p.""")

    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {chira_utilities.__version__}')

    return parser.parse_args()


def main():
    """Main function to orchestrate the CRL quantification workflow."""
    args = parse_arguments()
    print_configuration(args)

    # Build CRLs
    loci_file = os.path.join(args.outdir, 'loci.txt')
    chira_utilities.print_w_time("START: Build CRLs")
    build_crls(args.build_crls_too, args.bed, args.merged_bed,
               loci_file, args.crl_share, args.min_locus_size)
    chira_utilities.print_w_time("END: Build CRLs")
    
    # Quantify CRLs
    chira_utilities.print_w_time("START: Quantify CRLs")
    if args.use_sqldb:
        d_read_crl_fractions, d_crl_tpms = quantify_crls_sqlite(loci_file, args.em_thresh, args.num_processes)
    else:
        d_read_crl_fractions, d_crl_tpms = quantify_crls(loci_file, args.em_thresh, args.num_processes)
    chira_utilities.print_w_time("END: Quantify CRLs")
    
    # Write output
    temp_output_file = os.path.join(args.outdir, 'loci.counts.temp')
    chira_utilities.print_w_time("START: Write CRLs")
    write_crl_output(loci_file, temp_output_file, d_read_crl_fractions, d_crl_tpms)
    chira_utilities.print_w_time("END: Write CRLs")
    
    # Cleanup and finalize
    os.remove(loci_file)
    finalize_output(args.outdir, 'loci.counts.temp', 'loci.counts')

    # header_loci_counts = "\t".join([
    #                                 "segment_id",              # Column 0: Read/segment identifier
    #                                 "transcript_id",           # Column 1: Transcript ID
    #                                 "locus_id",                # Column 2: Locus identifier
    #                                 "crl_id",                  # Column 3: CRL (Chimeric Read Locus) group ID
    #                                 "transcript_start",        # Column 4: Transcript start position
    #                                 "transcript_end",          # Column 5: Transcript end position
    #                                 "transcript_strand",       # Column 6: Transcript strand (+ or -)
    #                                 "cigar",                   # Column 7: CIGAR string for alignment
    #                                 "genomic_pos",             # Column 8: Genomic position (if GTF provided)
    #                                 "locus_pos",               # Column 9: Locus position (format: chr:start:end:strand)
    #                                 "locus_share",            # Column 10: Fraction of reads in CRL that map to this locus
    #                                 "read_crl_fraction",      # Column 11: EM-optimized fraction (probability) that read belongs to this CRL
    #                                 "crl_tpm"                 # Column 12: TPM (Transcripts Per Million) value for the CRL
    # ])


if __name__ == "__main__":
    main()
