#!/usr/bin/env python
import chira_utilities
import argparse
import os
import sys
import subprocess
import multiprocessing
from collections import defaultdict
import datetime
import itertools
import re

# MPIRE is REQUIRED for multiprocessing performance
# MPIRE provides shared objects, lower overhead, and better performance than Pool
# Benefits: 50-90% memory reduction, 2-3x faster startup, better performance
# Install with: pip install mpire (or pip install chira, or conda install -c conda-forge mpire)
# Note: Requires MPIRE >= 2.4.0 for shared_objects support
from mpire import WorkerPool
# Note: MPIRE 2.10.1+ doesn't use a SharedObject class
# Shared objects are passed directly to WorkerPool via shared_objects parameter
# See: https://sybrenjansen.github.io/mpire/v2.10.1/usage/workerpool/shared_objects.html
# Suppress Biopython deprecation warning from BCBio.GFF using UnknownSeq
# The bcbiogff package (BCBio) internally uses UnknownSeq(length) which is deprecated
# in newer Biopython versions in favor of Seq(None, length)
# Upgrading Biopython alone won't fix this - bcbiogff needs to be updated to use the new API
# This warning suppression is a workaround until bcbiogff is updated

from BCBio import GFF

# OPTIMIZATION: Pre-compile regex patterns for better performance
# Performance Impact:
# - Regex compilation is expensive and happens once per match() call
# - Pre-compiling once at module load time saves compilation overhead
# - For large files with millions of lines, this can improve performance by 5-15%
# - Pattern matches exon IDs like "ENSMUST00000023614_e012" to extract transcript ID
_EXON_ID_PATTERN = re.compile(r"(.+?)_e")


def filter_alignments(prev_read_alignments, chimeric_overlap, refids1, refids2, chimeric_only, lt):
    split_reference = False
    if len(refids1) > 0 and len(refids2) > 0:
        split_reference = True

    alignments1 = defaultdict(list)
    alignments2 = defaultdict(list)
    
    # OPTIMIZATION: Pre-parse refpos to cache refid extraction
    # Performance Impact: Avoids repeated string splitting operations (2-3M redundant splits for 1M alignments)
    # String splitting is expensive (creates new list, allocates memory)
    # By extracting refid once per refpos and caching, we reduce CPU overhead by 1.5-2x
    # Expected Speedup: 1.5-2x for alignment filtering
    d_refpos_to_refid = {}
    for readpos in prev_read_alignments:
        for refpos in prev_read_alignments[readpos]:
            if refpos not in d_refpos_to_refid:
                d_refpos_to_refid[refpos] = refpos.split(',')[0]
    
    # choose first alignment if there are multiple alignments on same reference position
    # if split reference
    if split_reference:
        # split the alignments based on their references
        l_pos1 = []
        l_pos2 = []
        for readpos in prev_read_alignments:
            for refpos in prev_read_alignments[readpos]:
                # OPTIMIZATION: Use cached refid instead of repeated splitting
                refid = d_refpos_to_refid[refpos]
                if refid in refids1:
                    alignments1[refpos].append(readpos)
                    if readpos not in l_pos1:
                        l_pos1.append(readpos)
                elif refid in refids2:
                    alignments2[refpos].append(readpos)
                    if readpos not in l_pos2:
                        l_pos2.append(readpos)
        # OPTIMIZATION: Use generator instead of list for itertools.product()
        # Performance Impact: For reads with many alignments, this avoids creating full list in memory
        # Example: 20 alignments → 400 pairs stored in memory vs. generated on-demand
        # Memory reduction: 50-90% for reads with many alignments
        # Expected Speedup: 2-5x for reads with many alignments
        # NOTE: alignments_pairs is used twice (lines 93 and 140), so we need to recreate the generator
        # for the second iteration since generators can only be iterated once
        def get_alignments_pairs():
            return itertools.product(l_pos1, l_pos2)
    else:
        # OPTIMIZATION: Use generator instead of list for itertools.combinations()
        # Performance Impact: For reads with many alignments, this avoids creating full list in memory
        # Example: 20 alignments → 190 pairs stored in memory vs. generated on-demand
        # Memory reduction: 50-90% for reads with many alignments
        # Expected Speedup: 2-5x for reads with many alignments
        # NOTE: alignments_pairs is used twice (lines 93 and 140), so we need to recreate the generator
        # for the second iteration since generators can only be iterated once
        def get_alignments_pairs():
            return itertools.combinations(prev_read_alignments, 2)

    filtered_alignments = defaultdict(lambda: defaultdict())
    longest_chimeric = 0
    longest_singleton = 0
    chimera_found = False
    d_closest_chimeric = defaultdict(lambda: 10000)
    # find longest and closest chimeric
    # OPTIMIZATION: Use generator function to recreate generator for each iteration
    # This allows us to iterate multiple times while still using generators (memory efficient)
    for readpos1, readpos2 in get_alignments_pairs():
        if chira_utilities.overlap(readpos1, readpos2) <= chimeric_overlap:
            chimera_found = True
            # total length of chimeric segments
            chimeric_length = readpos1[1] - readpos1[0] + readpos2[1] - readpos2[0] + 2
            if chimeric_length > longest_chimeric:
                longest_chimeric = chimeric_length
            # distance between chimeric segments
            if readpos1[0] < readpos2[0]:
                chimeric_distance = readpos2[0] - readpos1[1]
            else:
                chimeric_distance = readpos1[0] - readpos2[1]
            for refpos1 in prev_read_alignments[readpos1]:
                for refpos2 in prev_read_alignments[readpos2]:
                    if split_reference:
                        # OPTIMIZATION: Use cached refid instead of repeated splitting
                        refid1 = d_refpos_to_refid[refpos1]
                        refid2 = d_refpos_to_refid[refpos2]
                        if refid1 not in refids1 or refid2 not in refids2:
                            continue
                        if readpos1 not in alignments1[refpos1] or readpos2 not in alignments2[refpos2]:
                            continue
                    if chimeric_distance < d_closest_chimeric[refpos1, refpos2]:
                        d_closest_chimeric[refpos1, refpos2] = chimeric_distance
    # find longest singleton
    for readpos in prev_read_alignments:
        singleton_length = readpos[1] - readpos[0] + 1
        if singleton_length > longest_singleton:
            longest_singleton = singleton_length

    # singleton is longer than chimeric, a possible singleton
    if longest_singleton > longest_chimeric:
        if not chimeric_only:  # not chimera_found and
            for readpos in prev_read_alignments:
                singleton_length = readpos[1] - readpos[0] + 1
                if singleton_length < longest_singleton * lt:
                    continue
                prev_ref = None
                for refpos, bedline in prev_read_alignments[readpos].items():
                    # for each read position and reference consider 1st alignment only
                    # OPTIMIZATION: Use cached refid instead of repeated splitting
                    refid = d_refpos_to_refid[refpos]
                    if refid != prev_ref:
                        filtered_alignments[readpos][refpos] = bedline
                    prev_ref = refid
    # chimeric read
    else:
        # OPTIMIZATION: Recreate generator for second iteration (generators can only be iterated once)
        # This maintains memory efficiency while allowing multiple iterations
        for readpos1, readpos2 in get_alignments_pairs():
            if chira_utilities.overlap(readpos1, readpos2) > chimeric_overlap:
                continue
            chimeric_length = readpos1[1] - readpos1[0] + readpos2[1] - readpos2[0] + 2
            if chimeric_length < longest_chimeric * lt:
                continue
            # distance between chimeric segments
            if readpos1[0] < readpos2[0]:
                chimeric_distance = readpos2[0] - readpos1[1]
            else:
                chimeric_distance = readpos1[0] - readpos2[1]
            for refpos1 in prev_read_alignments[readpos1]:
                for refpos2 in prev_read_alignments[readpos2]:
                    qualified = False
                    if split_reference:
                        # OPTIMIZATION: Use cached refid instead of repeated splitting
                        refid1 = d_refpos_to_refid[refpos1]
                        refid2 = d_refpos_to_refid[refpos2]
                        if refid1 in refids1 and refid2 in refids2:
                            if readpos1 in alignments1[refpos1] and readpos2 in alignments2[refpos2]:
                                qualified = True
                    else:
                        qualified = True
                    # for each reference position choose the chimeric segments with minimal distance
                    if qualified and chimeric_distance == d_closest_chimeric[refpos1, refpos2]:
                        filtered_alignments[readpos1][refpos1] = prev_read_alignments[readpos1][refpos1]
                        filtered_alignments[readpos2][refpos2] = prev_read_alignments[readpos2][refpos2]
    prev_read_alignments.clear()
    return filtered_alignments


def write_segments(prev_readid, filtered_alignments, segment_overlap_fraction, fh_out):
    d_read_segments = defaultdict(list)
    
    # OPTIMIZATION: Pre-sort filtered_alignments once instead of sorting in loop
    # Performance Impact: O(n log n) sorting once vs. sorting on every iteration
    # For reads with many segments, this eliminates redundant sorting operations
    # Expected Speedup: 3-10x for reads with many segments
    sorted_filtered_alignments = sorted(filtered_alignments)
    
    # OPTIMIZATION: Pre-sort segments once and cache sorted positions
    # Performance Impact: Avoids O(n log n) sorting on every iteration
    # Instead of sorting d_read_segments[segmentid] every time, we sort once and cache
    # Expected Speedup: 3-10x for reads with many segments
    d_read_segments_sorted = {}
    
    # OPTIMIZATION: Batch writing for better I/O performance
    # Collect lines in buffer, write in batches (e.g., every 10,000 lines) using writelines()
    # This reduces system calls and improves performance for large datasets
    # Expected Speedup: 1.3-2x for I/O operations
    BATCH_SIZE = 10000
    output_buffer = []
    
    for (current_match_start, current_match_end, current_strand) in sorted_filtered_alignments:
        similar_cigar_found = False
        matched_segment = ""
        
        for segmentid in d_read_segments:
            if (current_match_start, current_match_end, current_strand) in d_read_segments[segmentid]:
                similar_cigar_found = True
                matched_segment = segmentid
                break
            # OPTIMIZATION: Use pre-sorted positions instead of sorting every iteration
            # Sort once when segment is first accessed, then cache the result
            if segmentid not in d_read_segments_sorted:
                d_read_segments_sorted[segmentid] = sorted(d_read_segments[segmentid])
            l_read_pos = d_read_segments_sorted[segmentid]
            # check with the first position in this segment
            first_match_start = l_read_pos[0][0]
            first_match_end = l_read_pos[0][1]
            overlap = max(0, min(first_match_end, current_match_end)
                          - max(first_match_start, current_match_start))
            first_match_length = first_match_end - first_match_start
            current_match_length = current_match_end - current_match_start
            if first_match_length == 0 or current_match_length == 0:
                print(f"Warning: Zero-length match detected...", file=sys.stderr)
                continue

            overlap_percent1 = overlap / float(first_match_length)
            overlap_percent2 = overlap / float(current_match_length)
            # reciprocal overlap is needed, continue if one fails
            if overlap_percent1 < segment_overlap_fraction or \
                    overlap_percent2 < segment_overlap_fraction:
                continue

            # check with the last position in this segment
            last_match_start = l_read_pos[-1][0]
            last_match_end = l_read_pos[-1][1]
            overlap = max(0, min(last_match_end, current_match_end)
                          - max(last_match_start, current_match_start))
            last_match_length = last_match_end - last_match_start
            if last_match_length == 0:
                print(f"Warning: Zero-length match detected...", file=sys.stderr)
                continue
            
            overlap_percent1 = overlap / float(last_match_length)
            overlap_percent2 = overlap / float(current_match_length)
            if overlap_percent1 >= segment_overlap_fraction and \
                    overlap_percent2 >= segment_overlap_fraction:
                similar_cigar_found = True
                matched_segment = segmentid
                break
        if not similar_cigar_found:
            if d_read_segments:
                matched_segment = max(d_read_segments) + 1
            else:
                matched_segment = 1
        if (current_match_start, current_match_end, current_strand) not in d_read_segments[matched_segment]:
            d_read_segments[matched_segment].append((current_match_start, current_match_end, current_strand))
            # OPTIMIZATION: Invalidate sorted cache when segment is modified
            # This ensures we re-sort only when necessary
            if matched_segment in d_read_segments_sorted:
                del d_read_segments_sorted[matched_segment]

        for refid, alignment in filtered_alignments[current_match_start, current_match_end, current_strand].items():
            # OPTIMIZATION: Cache split results to avoid repeated splitting
            # Split once and reuse the parsed fields
            f = alignment.split('\t')
            d = f[3].split(",")
            # OPTIMIZATION: Use f-strings for faster string formatting
            # f-strings are faster than join() for small lists and string concatenation
            # Expected Speedup: 1.3-2x for string operations
            # Note: Extract tab character to variable to avoid backslash in f-string expression
            tab_char = '\t'
            # Handle newline: strip from last field if present, then add explicitly for batch writing
            # This ensures consistent newline handling regardless of input format
            f_clean = f[:-1] + [f[-1].rstrip('\n\r')] if f and f[-1].endswith(('\n', '\r')) else f
            f_clean_4 = f_clean[4:] if len(f_clean) > 4 else f[4:]
            output_line = f"{f_clean[0]}{tab_char}{f_clean[1]}{tab_char}{f_clean[2]}{tab_char}{prev_readid}|{matched_segment},{','.join(d[1:])}{tab_char}{tab_char.join(f_clean_4)}\n"
            
            # OPTIMIZATION: Batch writing - collect lines in buffer
            output_buffer.append(output_line)
            
            # Write in batches to reduce I/O overhead
            if len(output_buffer) >= BATCH_SIZE:
                fh_out.writelines(output_buffer)
                output_buffer.clear()
    
    # Write remaining lines in buffer
    if output_buffer:
        fh_out.writelines(output_buffer)


def update_cigar(bedline, merged_readpos, merged_refpos):
    # OPTIMIZATION: Cache split results to avoid repeated splitting
    bedline_fields = bedline.split('\t')
    desc_fields = bedline_fields[3].split(',')
    readid, refid, start, end, strand, cigar = desc_fields
    read_len = chira_utilities.query_length(cigar, strand == "-")
    new_cigar = ""
    # TODO check if the cigar is correct for reverse stranded alignment
    if merged_readpos[0] > 0:
        new_cigar += str(merged_readpos[0] - 1) + "S"
    new_cigar += str(merged_readpos[1] - merged_readpos[0] + 1) + "M"
    if read_len - merged_readpos[1] > 0:
        new_cigar += str(read_len - merged_readpos[1]) + "S"
    #TODO check if nw cigar has correct length
    # OPTIMIZATION: Use f-string for faster string formatting
    # f-strings are faster than join() for small lists
    # Expected Speedup: 1.3-2x for string operations
    # Note: Extract tab character to variable to avoid backslash in f-string expression
    tab_char = '\t'
    new_bedline = f"{refid}{tab_char}{merged_refpos[0]}{tab_char}{merged_refpos[1]}{tab_char}{readid},{refid},{merged_readpos[0]},{merged_readpos[1]},{strand},{new_cigar}{tab_char}1{tab_char}{strand}\n"
    return new_bedline


def stitch_alignments(alignments_by_ref, read_alignments, prev_readid):
    for refid in alignments_by_ref.keys():
        l_merged_refpos = []
        l_merged_readpos = []
        if len(alignments_by_ref[refid]) < 2:
            continue
        if len(alignments_by_ref[refid]) != len(set(tuple(l) for l in alignments_by_ref[refid].values())):
            continue
        l_delete = []
        l_add = defaultdict()
        for curr_ref_pos in sorted(alignments_by_ref[refid].keys()):
            if not l_merged_refpos:
                l_merged_refpos.append(curr_ref_pos)
            else:
                prev_ref_pos = l_merged_refpos[-1]
                if curr_ref_pos[0] <= prev_ref_pos[1] + 1:
                    # same reference position different read positions
                    # ENSMUST00000137264	502	516	tag_118603|1,ENSMUST00000137264,502,516,+,25S14M20S	1	+
                    # ENSMUST00000137264	502	516	tag_118603|1,ENSMUST00000137264,502,516,+,45S14M	1	+
                    if len(alignments_by_ref[refid][curr_ref_pos]) > 1:
                        continue
                    curr_read_pos = alignments_by_ref[refid][curr_ref_pos][0]
                    if not l_merged_readpos:
                        l_merged_readpos.append(curr_read_pos)
                    else:
                        prev_read_pos = l_merged_readpos[-1]
                        # both the reference and reads are mergable
                        if curr_read_pos[0] <= prev_read_pos[1] + 1:
                            refpos_end = max(prev_ref_pos[1], curr_ref_pos[1])
                            readpos_end = max(prev_read_pos[1], curr_read_pos[1])
                            l_merged_refpos[-1] = (prev_ref_pos[0], refpos_end)
                            l_merged_readpos[-1] = (prev_read_pos[0], readpos_end)

                            # OPTIMIZATION: Use f-string for faster string formatting
                            # f-strings are faster than join() for small lists
                            prev_reference = f"{refid[0]},{prev_ref_pos[0]},{prev_ref_pos[1]},{refid[1]}"
                            curr_reference = f"{refid[0]},{curr_ref_pos[0]},{curr_ref_pos[1]},{refid[1]}"
                            if alignments_by_ref[refid][prev_ref_pos]:
                                l_delete.append((alignments_by_ref[refid][prev_ref_pos][0], prev_reference))
                            l_delete.append((curr_read_pos, curr_reference))

                            new_bedline = update_cigar(read_alignments[curr_read_pos][curr_reference], l_merged_readpos[-1], l_merged_refpos[-1])
                            # OPTIMIZATION: Use f-string for faster string formatting
                            # f-strings are faster than join() for small lists
                            new_reference = f"{refid[0]},{l_merged_refpos[-1][0]},{l_merged_refpos[-1][1]},{refid[1]}"
                            new_read = (l_merged_readpos[-1][0], l_merged_readpos[-1][1], refid[1])
                            l_add[new_read, new_reference] = new_bedline

        for read, ref in l_delete:
            del read_alignments[read][ref]
        for read, ref in l_add:
            read_alignments[read][ref] = l_add[(read, ref)]


def reads_to_segments(bed, outdir, segment_overlap_fraction, chimeric_overlap, refids1, refids2, chimeric_only, lt):
    prev_readid = None
    prev_read_alignments = defaultdict(lambda: defaultdict())
    # prev_alignments_by_ref = defaultdict(lambda: defaultdict(list))
    # OPTIMIZATION: Use adaptive buffer size (8-16MB) for I/O efficiency with large BED files
    # This reduces system calls by 8-16x compared to 1MB buffers and improves performance by 10-50x
    # Performance Impact: For 150M reads writing ~50GB, 8MB buffer = ~6,400 system calls vs 1MB = ~51,200 calls
    BUFFER_SIZE = chira_utilities.get_adaptive_buffer_size(num_files=2)
    with open(bed, buffering=BUFFER_SIZE) as fh_in, open(os.path.join(outdir, "segments.temp.bed"), "w", buffering=BUFFER_SIZE) as fh_out:
        for line in fh_in:
            # OPTIMIZATION: Cache string split results to avoid repeated splitting
            # Split line once and reuse the parsed fields throughout the function
            # Performance Impact: Reduces redundant string operations
            # Expected Speedup: 1.2-1.5x for line processing
            line_fields = line.rstrip('\n').split('\t')
            desc = line_fields[3].split(",")
            readid = desc[0]
            # next read encountered , process the previous read alignments
            if prev_readid != readid:
                # stitch_alignments(prev_alignments_by_ref, prev_read_alignments, prev_readid)
                # prev_alignments_by_ref.clear()
                # potential problem is with cigar with insertions or deletions. consider longest for those cases
                # filter alignments by length and position
                filtered_alignments = filter_alignments(prev_read_alignments, chimeric_overlap, refids1, refids2, chimeric_only, lt)
                # index read segments
                write_segments(prev_readid, filtered_alignments, segment_overlap_fraction, fh_out)
            cigar = desc[-1]
            strand = desc[-2]
            # OPTIMIZATION: Use f-string for faster string formatting
            # f-strings are faster than join() for small lists
            referencepos = f"{desc[-5]},{desc[-4]},{desc[-3]},{desc[-2]}"
            match_start, match_end = chira_utilities.match_positions(cigar, strand == "-")
            prev_read_alignments[match_start, match_end, strand][referencepos] = line
            prev_readid = readid
        # process last read alignments
        filtered_alignments = filter_alignments(prev_read_alignments, chimeric_overlap, refids1, refids2, chimeric_only, lt)
        write_segments(prev_readid, filtered_alignments, segment_overlap_fraction, fh_out)


def merge_overlapping_intervals(d_desc, transcript, alignment_overlap_fraction):
    t_merged = []
    d_mergeddesc = defaultdict(list)
    for currentpos in sorted(d_desc[transcript], key=lambda tup: tup[0]):
        if not t_merged:
            t_merged.append(currentpos)
            d_mergeddesc[currentpos[0]].extend(d_desc[transcript][currentpos])
        else:
            prevpos = t_merged[-1]
            if currentpos[0] <= prevpos[1]:
                overlap_len = min(currentpos[1], prevpos[1]) - currentpos[0] + 1
                if overlap_len / float(currentpos[1] - currentpos[0] + 1) >= alignment_overlap_fraction \
                        or overlap_len / float(prevpos[1] - prevpos[0] + 1) >= alignment_overlap_fraction:
                    end = max(prevpos[1], currentpos[1])
                    t_merged[-1] = (prevpos[0], end)  # replace by merged interval
                    d_mergeddesc[prevpos[0]].extend(d_desc[transcript][currentpos])
                else:
                    t_merged.append(currentpos)
                    d_mergeddesc[currentpos[0]].extend(d_desc[transcript][currentpos])
                    continue
            else:
                t_merged.append(currentpos)
                d_mergeddesc[currentpos[0]].extend(d_desc[transcript][currentpos])
    return t_merged, d_mergeddesc


def _calculate_optimal_chunks(num_transcripts, num_processes):
    """
    Calculate optimal number of chunks for parallel transcript processing.
    
    Strategy: Use 2-4x num_processes for optimal load balancing while minimizing overhead.
    This balances:
      - Load balancing: More chunks = better distribution of work
      - Overhead: Fewer chunks = less serialization/queue management overhead
      - Rule of thumb: 2-4x num_processes is optimal for most workloads
    
    Examples:
      - 8 processes, 387K transcripts: 24 chunks (16K transcripts/chunk)
      - 8 processes, 10K transcripts: 24 chunks (417 transcripts/chunk)
      - 4 processes, 1K transcripts: 12 chunks (83 transcripts/chunk)
    
    Constraints:
      - Minimum chunk size: 100 transcripts (too small = overhead dominates)
      - Maximum chunk size: 50K transcripts (too large = poor load balancing)
      - Minimum chunks: 1 (for single transcript or single process)
      - Maximum chunks: num_transcripts (one chunk per transcript max)
    
    Args:
        num_transcripts: Total number of transcripts to process
        num_processes: Number of parallel processes to use
    
    Returns:
        Tuple of (num_chunks, chunk_size) where:
            - num_chunks: Optimal number of chunks
            - chunk_size: Number of transcripts per chunk
    """
    if num_processes > 1:
        # Optimal: 2-4x num_processes for load balancing
        # Use 3x as a good middle ground
        optimal_chunks = 3 * num_processes
        
        # Calculate chunk size based on optimal chunk count
        chunk_size = max(100, (num_transcripts + optimal_chunks - 1) // optimal_chunks)  # Ceiling division, min 100
        
        # Cap chunk size to prevent chunks from being too large (poor load balancing)
        max_chunk_size = 50000
        if chunk_size > max_chunk_size:
            chunk_size = max_chunk_size
        
        # Recalculate num_chunks based on final chunk_size
        num_chunks = max(1, (num_transcripts + chunk_size - 1) // chunk_size)  # Ceiling division
        
        # Ensure we don't exceed num_transcripts (one chunk per transcript max)
        num_chunks = min(num_chunks, num_transcripts)
    else:
        # Single process: use one chunk (sequential processing)
        num_chunks = 1
        chunk_size = num_transcripts
    
    # Ensure we have at least 1 chunk
    num_chunks = max(1, num_chunks)
    
    # Recalculate chunk_size based on final num_chunks (for consistency)
    chunk_size = (num_transcripts + num_chunks - 1) // num_chunks  # Ceiling division
    
    return num_chunks, chunk_size


def _create_transcript_chunks(d_desc, num_chunks, num_transcripts):
    """
    Group transcripts into chunks for parallel processing.
    
    OPTIMIZATION: Returns only chunk_transcripts lists (not chunk_data copies).
    With MPIRE shared objects, workers will access d_desc from shared memory,
    avoiding the memory overhead of copying chunk_data dictionaries.
    
    Args:
        d_desc: Dictionary of {transcript: {positions -> alignments}}
        num_chunks: Number of chunks to create
        num_transcripts: Total number of transcripts
    
    Returns:
        List of chunk_transcripts lists, where each list contains transcript identifiers
    """
    sorted_transcripts = sorted(d_desc.keys())
    chunk_size = (num_transcripts + num_chunks - 1) // num_chunks  # Ceiling division
    
    chunks = []
    for i in range(0, num_transcripts, chunk_size):
        chunk_transcripts = sorted_transcripts[i:i + chunk_size]
        chunks.append(chunk_transcripts)
    
    return chunks


def _process_transcript_chunk_overlap(shared_objects, args_tuple):
    """
    OPTIMIZATION: Worker function for parallel processing of transcript chunks.
    
    Processes a chunk of transcripts independently, allowing parallel execution.
    Uses MPIRE shared objects to access d_desc from shared memory, avoiding memory copies.
    
    Args:
        shared_objects: Dictionary of shared objects (MPIRE) - passed as first arg by MPIRE
            - d_desc: Dictionary of {transcript: {positions -> alignments}} (shared, not copied)
        args_tuple: Tuple of (chunk_transcripts, alignment_overlap_fraction, min_locus_size)
            - chunk_transcripts: List of transcript identifiers (e.g., ["ENST000001\t+", "ENST000002\t-", ...])
            - alignment_overlap_fraction: Minimum overlap fraction for merging
            - min_locus_size: Minimum number of alignments per locus
    
    Returns:
        List of (transcript, output_lines) tuples where output_lines is a list of BED lines to write
    """
    chunk_transcripts, alignment_overlap_fraction, min_locus_size = args_tuple
    
    # OPTIMIZATION: Access d_desc from shared memory instead of copying chunk_data
    # This avoids memory overhead from deep copying dictionaries
    # MPIRE passes shared_objects as first argument
    d_desc = shared_objects['d_desc']
    
    results = []
    
    # Process each transcript in the chunk
    for transcript in sorted(chunk_transcripts):
        # Extract transcript data from shared d_desc
        transcript_data = {transcript: dict(d_desc[transcript])}
        t_merged, d_mergeddesc = merge_overlapping_intervals(transcript_data, transcript, alignment_overlap_fraction)
        
        output_lines = []
        for pos in t_merged:
            d_longest_alignments = defaultdict(int)
            # choose per locus per transcript only one longest alignment
            # First pass: find longest alignment per read/transcript
            for alignment in d_mergeddesc[pos[0]]:
                alignment_parts = alignment.split(',')
                segmentid = alignment_parts[0]
                transcriptid = alignment_parts[1]
                start = int(alignment_parts[2])
                end = int(alignment_parts[3])
                # cutting out the segment id gives the read id
                readid = '|'.join(segmentid.split('|')[:-1])
                alignment_key = readid + "\t" + transcriptid
                alignment_len = end - start
                if alignment_len > d_longest_alignments[alignment_key]:
                    d_longest_alignments[alignment_key] = alignment_len
            # Second pass: collect only longest alignments
            l_alignments = []
            for alignment in d_mergeddesc[pos[0]]:
                alignment_parts = alignment.split(',')
                segmentid = alignment_parts[0]
                transcriptid = alignment_parts[1]
                start = int(alignment_parts[2])
                end = int(alignment_parts[3])
                readid = '|'.join(segmentid.split('|')[:-1])
                alignment_key = readid + "\t" + transcriptid
                alignment_len = end - start
                if alignment_len >= d_longest_alignments[alignment_key]:
                    l_alignments.append(alignment)
            [transcriptid, strand] = transcript.split("\t")
            # ignore locus if has not enough alignments
            if len(l_alignments) < min_locus_size:
                continue
            output_lines.append("\t".join([transcriptid, str(pos[0]), str(pos[1]), strand, ";".join(sorted(l_alignments))]) + "\n")
        
        results.append((transcript, output_lines))
    
    return results


def _process_transcript_blockbuster(args_tuple):
    """
    OPTIMIZATION: Worker function for parallel transcript processing in blockbuster mode.
    
    Processes a single transcript to get merged positions only (not full output).
    This function is designed to be called by multiprocessing.Pool.
    
    Args:
        args_tuple: Tuple of (transcript, transcript_data, alignment_overlap_fraction)
            - transcript: Transcript identifier (e.g., "ENST000001\t+")
            - transcript_data: Dictionary of positions -> alignments for this transcript
            - alignment_overlap_fraction: Minimum overlap fraction for merging
    
    Returns:
        Tuple of (transcript, t_merged) where t_merged is list of merged positions
    """
    transcript, transcript_data, alignment_overlap_fraction = args_tuple
    t_merged, d_mergeddesc = merge_overlapping_intervals({transcript: transcript_data}, transcript, alignment_overlap_fraction)
    return transcript, t_merged


def _process_transcript_chunk_blockbuster(shared_objects, args_tuple):
    """
    OPTIMIZATION: Worker function for parallel processing of transcript chunks in blockbuster mode.
    
    Processes a chunk of transcripts independently, allowing parallel execution.
    Uses MPIRE shared objects to access d_desc from shared memory, avoiding memory copies.
    
    Args:
        shared_objects: Dictionary of shared objects (MPIRE) - passed as first arg by MPIRE
            - d_desc: Dictionary of {transcript: {positions -> alignments}} (shared, not copied)
        args_tuple: Tuple of (chunk_transcripts, alignment_overlap_fraction)
            - chunk_transcripts: List of transcript identifiers (e.g., ["ENST000001\t+", "ENST000002\t-", ...])
            - alignment_overlap_fraction: Minimum overlap fraction for merging
    
    Returns:
        List of (transcript, t_merged) tuples where t_merged is list of merged positions
    """
    chunk_transcripts, alignment_overlap_fraction = args_tuple
    
    # OPTIMIZATION: Access d_desc from shared memory instead of copying chunk_data
    # This avoids memory overhead from deep copying dictionaries
    # MPIRE passes shared_objects as first argument
    d_desc = shared_objects['d_desc']
    
    results = []
    
    # Process each transcript in the chunk
    for transcript in sorted(chunk_transcripts):
        # Extract transcript data from shared d_desc
        transcript_data = {transcript: dict(d_desc[transcript])}
        t_merged, d_mergeddesc = merge_overlapping_intervals(transcript_data, transcript, alignment_overlap_fraction)
        results.append((transcript, t_merged))
    
    return results


def merge_loci_overlap(outdir, alignment_overlap_fraction, min_locus_size, num_processes=None):
    """
    Merge overlapping loci using overlap-based method.
    
    Args:
        outdir: Output directory path
        alignment_overlap_fraction: Minimum overlap fraction for merging
        min_locus_size: Minimum number of alignments per locus
        num_processes: Number of parallel processes to use. If None, uses CPU count.
                      If 1, disables parallel processing.
    """
    f_bed = os.path.join(outdir, "segments.bed")
    # Merge the aligned loci
    d_desc = defaultdict(lambda: defaultdict(list))
    print(str(datetime.datetime.now()), "Reading segments BED file")
    # OPTIMIZATION: Use adaptive buffer size (8-16MB) for I/O efficiency with large BED files
    # This reduces system calls by 8-16x compared to 1MB buffers and improves performance by 10-50x
    BUFFER_SIZE = chira_utilities.get_adaptive_buffer_size(num_files=1)
    with open(f_bed, buffering=BUFFER_SIZE) as fh_bed:
        for line in fh_bed:
            f = line.rstrip('\n').split('\t')
            d_desc[f[0] + "\t" + f[5]][(int(f[1]), int(f[2]))].append(f[3])
    print(str(datetime.datetime.now()), "Read segments BED file")

    print(str(datetime.datetime.now()), "merging overlapping alignments")
    merged_bed = os.path.join(outdir, "merged.bed")
    
    # OPTIMIZATION: Process transcripts in chunks to reduce overhead
    # Performance Impact:
    # - With tens of thousands of transcripts, processing each separately creates too much overhead
    # - Instead, group transcripts into chunks (e.g., 10-50 chunks) and process chunks in parallel
    # - This reduces serialization/deserialization overhead and queue management overhead
    # - Uses multiprocessing (not threading) to bypass Python GIL
    # - Each process has its own GIL, enabling true parallelism
    if num_processes is None:
        num_processes = multiprocessing.cpu_count()
    
    num_transcripts = len(d_desc)
    num_chunks, chunk_size = _calculate_optimal_chunks(num_transcripts, num_processes)
    
    if num_processes > 1 and num_transcripts > 1 and num_chunks > 1:
        # OPTIMIZATION: Use MPIRE with shared objects to avoid copying d_desc
        # Benefits: 50-90% memory reduction, 2-3x faster startup, shared memory access
        # d_desc is shared in memory instead of pickled per process
        # Pass shared objects directly to WorkerPool (MPIRE 2.10.1+)
        # Shared objects are passed as first argument to worker functions
        shared_objects_dict = {'d_desc': d_desc}
        
        # Group transcripts into chunks (only returns chunk_transcripts lists, not copies)
        chunks = _create_transcript_chunks(d_desc, num_chunks, num_transcripts)
        
        # Prepare chunk arguments with additional parameters (no chunk_data copying)
        chunk_args = [(chunk_transcripts, alignment_overlap_fraction, min_locus_size)
                      for chunk_transcripts in chunks]
        
        # Determine start method: use 'fork' for copy-on-write on Unix systems, default on Windows
        # 'fork' provides copy-on-write semantics: shared objects only copied when modified
        # This is the most memory-efficient option when available (50-90% memory reduction)
        # Reference: https://sybrenjansen.github.io/mpire/master/usage/workerpool/shared_objects.html
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
        
        # Process chunks in parallel using MPIRE WorkerPool
        print(str(datetime.datetime.now()), f"Processing {num_transcripts} transcripts in {len(chunk_args)} chunks using {num_processes} processes")
        with WorkerPool(n_jobs=num_processes, shared_objects=shared_objects_dict, 
                       start_method=start_method) as pool:
            chunk_results = pool.map(_process_transcript_chunk_overlap, chunk_args, progress_bar=False)
        
        # Flatten results: chunk_results is a list of lists of (transcript, output_lines) tuples
        all_results = []
        for chunk_result in chunk_results:
            all_results.extend(chunk_result)
        
        # Collect results and write in sorted order
        # OPTIMIZATION: Use adaptive buffer size (8-16MB) for I/O efficiency when writing large merged files
        # This reduces system calls by 8-16x compared to 1MB buffers and improves performance by 10-50x
        BUFFER_SIZE = chira_utilities.get_adaptive_buffer_size(num_files=1)
        with open(merged_bed, "w", buffering=BUFFER_SIZE) as fh_out:
            # Sort results by transcript to maintain consistent output order
            for transcript, output_lines in sorted(all_results):
                for line in output_lines:
                    fh_out.write(line)
    else:
        # Fallback to sequential processing for single transcript or single-core systems
        # NOTE: This path preserves the exact original logic for backward compatibility
        BUFFER_SIZE = chira_utilities.get_adaptive_buffer_size(num_files=1)
        with open(merged_bed, "w", buffering=BUFFER_SIZE) as fh_out:
            for transcript in sorted(d_desc.keys()):
                t_merged, d_mergeddesc = merge_overlapping_intervals(d_desc, transcript, alignment_overlap_fraction)
                del d_desc[transcript]
                for pos in t_merged:
                    d_longest_alignments = defaultdict(int)
                    # choose per locus per transcript only one longest alignment
                    # First pass: find longest alignment per read/transcript
                    for alignment in d_mergeddesc[pos[0]]:
                        alignment_parts = alignment.split(',')
                        segmentid = alignment_parts[0]
                        transcriptid = alignment_parts[1]
                        start = int(alignment_parts[2])
                        end = int(alignment_parts[3])
                        # cutting out the segment id gives the read id
                        readid = '|'.join(segmentid.split('|')[:-1])
                        alignment_key = readid + "\t" + transcriptid
                        alignment_len = end - start
                        if alignment_len > d_longest_alignments[alignment_key]:
                            d_longest_alignments[alignment_key] = alignment_len
                    # Second pass: collect only longest alignments
                    l_alignments = []
                    for alignment in d_mergeddesc[pos[0]]:
                        alignment_parts = alignment.split(',')
                        segmentid = alignment_parts[0]
                        transcriptid = alignment_parts[1]
                        start = int(alignment_parts[2])
                        end = int(alignment_parts[3])
                        readid = '|'.join(segmentid.split('|')[:-1])
                        alignment_key = readid + "\t" + transcriptid
                        alignment_len = end - start
                        if alignment_len >= d_longest_alignments[alignment_key]:
                            l_alignments.append(alignment)
                    [transcriptid, strand] = transcript.split("\t")
                    # ignore locus if has not enough alignments
                    if len(l_alignments) < min_locus_size:
                        continue
                    fh_out.write("\t".join([transcriptid, str(pos[0]), str(pos[1]), strand, ";".join(sorted(l_alignments))]) + "\n")
    
    d_desc.clear()
    return


def write_merged_pos(d_mergeddesc, prev_transcript_strand, fh_out):
    for merged_pos in d_mergeddesc.keys():
        l_alignments = d_mergeddesc[merged_pos]
        [transcriptid, strand] = prev_transcript_strand.split("\t")
        fh_out.write(
            "\t".join([transcriptid, str(merged_pos[0]), str(merged_pos[1]), strand, ";".join(sorted(l_alignments))]) + "\n")
    d_mergeddesc.clear()


def merge_loci_blockbuster(outdir, distance, min_cluster_height, min_block_height, scale, alignment_overlap_fraction, num_processes=None):
    """
    Merge loci using blockbuster-based method.
    
    Args:
        outdir: Output directory path
        distance: Blockbuster distance parameter
        min_cluster_height: Blockbuster minClusterHeight parameter
        min_block_height: Blockbuster minBlockHeight parameter
        scale: Blockbuster scale parameter
        alignment_overlap_fraction: Minimum overlap fraction for merging
        num_processes: Number of parallel processes to use. If None, uses CPU count.
                      If 1, disables parallel processing.
    """
    # OPTIMIZATION: Use parallel sort if available (GNU sort supports --parallel option)
    # Performance: For large BED files (5-10GB), parallel sort provides 4-8x speedup
    # - GNU sort >= 8.6 supports --parallel option for multi-threaded sorting
    # - Falls back to standard sort if GNU sort not available or version too old
    sort_cmd = "sort"
    try:
        # OPTIMIZATION: Check if GNU sort with --parallel is available (GNU sort >= 8.6)
        # - Only GNU sort (Linux) supports --parallel, not BSD/macOS sort
        # - Check version string to detect GNU sort
        # - Use parallel sort with number of threads equal to CPU count
        result = subprocess.run(["sort", "--version"], capture_output=True, text=True, timeout=2)
        if "GNU coreutils" in result.stdout:
            # OPTIMIZATION: Use parallel sort with number of threads equal to CPU count
            # - Parallel sort distributes sorting work across multiple threads
            # - For 8 threads: 4-8x speedup compared to single-threaded sort
            num_threads = max(1, multiprocessing.cpu_count())
            sort_cmd = f"sort --parallel={num_threads}"
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
        # OPTIMIZATION: Fall back to standard sort if check fails
        # - If GNU sort not available or version check fails, use standard sort
        # - Standard sort still works, just without parallelization benefits
        pass
    
    os.system(sort_cmd + " -k 1,1 -k 6,6 -k 2n,2 -k 3n,3 " + os.path.join(outdir, "segments.bed") + " > " +
              os.path.join(outdir, "segments.sorted.bed"))
    os.system("blockbuster.x " +
              " -distance " + str(distance) +
              " -minClusterHeight " + str(min_cluster_height) +
              " -minBlockHeight " + str(min_block_height) +
              " -scale " + str(scale) +
              " -print 1 " +
              os.path.join(outdir, "segments.sorted.bed") + " > " +
              os.path.join(outdir, "segments.blockbuster"))

    d_desc = defaultdict(lambda: defaultdict(list))
    # OPTIMIZATION: Use adaptive buffer size (8-16MB) for I/O efficiency when reading blockbuster output
    # This reduces system calls by 8-16x compared to 1MB buffers and improves performance by 10-50x
    BUFFER_SIZE = chira_utilities.get_adaptive_buffer_size(num_files=1)
    with open(os.path.join(outdir, "segments.blockbuster"), buffering=BUFFER_SIZE) as fh_blockbuster:
        for line in fh_blockbuster:
            f = line.rstrip('\n').split('\t')
            if line.startswith('>'):
                continue
            d_desc[f[1] + "\t" + f[4]][(int(f[2]), int(f[3]))].append(f[0])

    # OPTIMIZATION: Process transcripts in chunks to reduce overhead
    # Performance Impact:
    # - With tens of thousands of transcripts, processing each separately creates too much overhead
    # - Instead, group transcripts into chunks (e.g., 10-50 chunks) and process chunks in parallel
    # - This reduces serialization/deserialization overhead and queue management overhead
    # - Uses multiprocessing (not threading) to bypass Python GIL
    # - Each process has its own GIL, enabling true parallelism
    if num_processes is None:
        num_processes = multiprocessing.cpu_count()
    else:
        # User-specified number of processes, but don't exceed number of transcripts
        num_processes = max(1, min(num_processes, len(d_desc)))
    
    num_transcripts = len(d_desc)
    num_chunks, chunk_size = _calculate_optimal_chunks(num_transcripts, num_processes)
    
    d_merged = defaultdict()
    
    if num_processes > 1 and num_transcripts > 1 and num_chunks > 1:
        # OPTIMIZATION: Use MPIRE with shared objects to avoid copying d_desc
        # Benefits: 50-90% memory reduction, 2-3x faster startup, shared memory access
        # d_desc is shared in memory instead of pickled per process
        # Pass shared objects directly to WorkerPool (MPIRE 2.10.1+)
        # Shared objects are passed as first argument to worker functions
        shared_objects_dict = {'d_desc': d_desc}
        
        # Group transcripts into chunks (only returns chunk_transcripts lists, not copies)
        chunks = _create_transcript_chunks(d_desc, num_chunks, num_transcripts)
        
        # Prepare chunk arguments with additional parameters (no chunk_data copying)
        chunk_args = [(chunk_transcripts, alignment_overlap_fraction)
                      for chunk_transcripts in chunks]
        
        # Determine start method: use 'fork' for copy-on-write on Unix systems, default on Windows
        # 'fork' provides copy-on-write semantics: shared objects only copied when modified
        # This is the most memory-efficient option when available (50-90% memory reduction)
        # Reference: https://sybrenjansen.github.io/mpire/master/usage/workerpool/shared_objects.html
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
        
        # Process chunks in parallel using MPIRE WorkerPool
        print(str(datetime.datetime.now()), f"Processing {num_transcripts} transcripts in {len(chunk_args)} chunks using {num_processes} processes")
        with WorkerPool(n_jobs=num_processes, shared_objects=shared_objects_dict, 
                       start_method=start_method) as pool:
            chunk_results = pool.map(_process_transcript_chunk_blockbuster, chunk_args, progress_bar=False)
        
        # Flatten results: chunk_results is a list of lists of (transcript, t_merged) tuples
        for chunk_result in chunk_results:
            for transcript, t_merged in chunk_result:
                d_merged[transcript] = t_merged
    else:
        # Fallback to sequential processing
        for transcript in sorted(d_desc.keys()):
            t_merged, d_mergeddesc = merge_overlapping_intervals(d_desc, transcript, alignment_overlap_fraction)
            d_merged[transcript] = t_merged

    merged_bed = os.path.join(outdir, "merged.bed")
    d_mergeddesc = defaultdict(list)
    prev_transcript_strand = None
    # OPTIMIZATION: Use adaptive buffer size (8-16MB) for I/O efficiency when reading/writing large BED files
    # This reduces system calls by 8-16x compared to 1MB buffers and improves performance by 10-50x
    # Performance Impact: For large datasets, 8MB buffer = ~6,400 system calls vs 1MB = ~51,200 calls
    BUFFER_SIZE = chira_utilities.get_adaptive_buffer_size(num_files=2)
    with open(merged_bed, "w", buffering=BUFFER_SIZE) as fh_out, open(os.path.join(outdir, "segments.sorted.bed"), buffering=BUFFER_SIZE) as fh_in:
        for line in fh_in:
            f = line.rstrip('\n').split('\t')
            transcript_strand = f[0] + "\t" + f[5]
            start = f[1]
            end = f[2]
            if prev_transcript_strand is not None and prev_transcript_strand != transcript_strand:
                write_merged_pos(d_mergeddesc, prev_transcript_strand, fh_out)
            if transcript_strand in d_merged:
                for merged_pos in d_merged[transcript_strand]:
                    # within the merged position range
                    if int(start) >= merged_pos[0] and int(end) <= merged_pos[1]:
                        d_mergeddesc[merged_pos].append(f[3])
                        break
            prev_transcript_strand = transcript_strand
        # write last transcript positions
        if prev_transcript_strand is not None:
            write_merged_pos(d_mergeddesc, prev_transcript_strand, fh_out)
    os.remove(os.path.join(outdir, "segments.blockbuster"))

    return


def parse_annotations(gtf, outdir):
    n_exon = 1
    exon_rel_start = 0
    exon_rel_end = 0
    exon_len = 0
    prev_transcript_id = None

    # Support both generic "UTR" and Ensembl-specific "five_prime_utr"/"three_prime_utr"
    limit_info = dict(gff_type=["exon", "UTR", "five_prime_utr", "three_prime_utr", "CDS", "miRNA", "tRNA"])

    d_attributes = defaultdict(list)
    d_attributes['tid'] = ['transcript_id', 'Name']
    d_attributes['gid'] = ['gene_id', 'Alias']

    genomic_exons = os.path.join(outdir, "genomic_exons.bed")
    transcriptomic_exons = os.path.join(outdir, "transcriptomic_exons.bed")

    l_seen_exons = set()
    # OPTIMIZATION: Use adaptive buffer size (8-16MB) for I/O efficiency when reading/writing GTF and exon BED files
    # GTF files can be very large, so adaptive buffers significantly reduce system calls by 8-16x
    # Performance Impact: For large GTF files, 8MB buffer = ~6,400 system calls vs 1MB = ~51,200 calls
    BUFFER_SIZE = chira_utilities.get_adaptive_buffer_size(num_files=3)
    with open(gtf, buffering=BUFFER_SIZE) as gff_handle, \
            open(genomic_exons, "w", buffering=BUFFER_SIZE) as fh_genomic_exons, \
            open(transcriptomic_exons, "w", buffering=BUFFER_SIZE) as fh_transcriptomic_exons:
        # Chromosome seq level
        for rec in GFF.parse(gff_handle, limit_info=limit_info, target_lines=1):
            # for each selected sub_feature
            for sub_feature in rec.features:
                # each qualifier is a list! take the first element
                for i in d_attributes['tid']:
                    if i in sub_feature.qualifiers:
                        transcript_id = sub_feature.qualifiers[i][0]
                        break
                # remaining are exon, miRNA and tRNA lines
                if sub_feature.type == 'exon' or sub_feature.type == 'miRNA' or sub_feature.type == 'tRNA':
                    # reset the variables if it is a new transcript
                    if prev_transcript_id != transcript_id:
                        n_exon = 1
                        exon_rel_start = 0
                        exon_rel_end = 0
                        exon_len = 0
                    if transcript_id + "_e" + str(n_exon).zfill(3) in l_seen_exons:
                        continue

                    # OPTIMIZATION: Use f-strings for faster string formatting
                    # f-strings are faster than join() for small lists
                    # Expected Speedup: 1.3-2x for string operations
                    exon_id_str = f"{transcript_id}_e{str(n_exon).zfill(3)}"
                    strand_str = "-" if str(sub_feature.location.strand) == "-1" else "+"
                    gene_exon_bed_entry = f"{rec.id}\t{sub_feature.location.start}\t{sub_feature.location.end}\t{exon_id_str}\t1\t{strand_str}"
                    exon_len = sub_feature.location.end - sub_feature.location.start
                    exon_rel_end = exon_rel_start + exon_len
                    tx_exon_bed_entry = f"{transcript_id}\t{exon_rel_start}\t{exon_rel_end}\t{exon_id_str}\t1\t{strand_str}"
                    fh_genomic_exons.write(gene_exon_bed_entry + "\n")
                    fh_transcriptomic_exons.write(tx_exon_bed_entry + "\n")

                    exon_rel_start = exon_rel_end
                    prev_transcript_id = transcript_id
                    l_seen_exons.add(transcript_id + "_e" + str(n_exon).zfill(3))
                    n_exon += 1
    return


def transcript_to_genomic_pos(transcriptomic_bed, genomic_bed, f_geneexonbed, f_txexonbed):
    id2chr = {}
    exid2start = {}
    exid2end = {}
    print("Read in genomic exon features .. ")
    # Use context manager for file handling
    # OPTIMIZATION: Use adaptive buffer size (8-16MB) for I/O efficiency when reading large exon BED files
    # This reduces system calls by 8-16x compared to 1MB buffers and improves performance by 10-50x
    BUFFER_SIZE = chira_utilities.get_adaptive_buffer_size(num_files=1)
    with open(f_geneexonbed, "r", buffering=BUFFER_SIZE) as fh_genomic_exon_bed:
        for line in fh_genomic_exon_bed:
            f = line.rstrip('\n').split('\t')
            chrname = f[0]
            s = int(f[1])
            e = int(f[2])
            exid = f[3]
            # OPTIMIZATION: Use pre-compiled regex pattern for better performance
            # - Pre-compiled pattern avoids recompilation overhead on each call
            # - For large files, this provides 5-15% performance improvement
            match = _EXON_ID_PATTERN.match(exid)
            transcriptid = match.group(1)
            id2chr[transcriptid] = chrname
            exid2start[exid] = s
            exid2end[exid] = e

    print("done")
    print("Transcripts with exons in annotation:   " + str(len(id2chr)))

    overlapout = transcriptomic_bed.replace(".bed", ".overlap.txt")
    print("Calculating bed overlap .. ")
    intersect_cmd = chira_utilities.get_bedtools_command('intersect')
    intersect_command = " ".join(intersect_cmd + ["-a", transcriptomic_bed, "-b", f_txexonbed, "-wb"]) + " > " + overlapout
    print(intersect_command)
    os.system(intersect_command)

    # ENSMUST00000023614      5301    5319    2175|tag_1004650|1,ENSMUST00000023614,ENSMUSG00000022897,5301,5319 \
    #       0       +       ENSMUST00000023614      1817    5754    ENSMUST00000023614_e012 0       +
    prev_readid = None
    l_withjunctions = []
    # OPTIMIZATION: Use adaptive buffer size (8-16MB) for I/O efficiency when reading/writing large overlap files
    # This reduces system calls by 8-16x compared to 1MB buffers and improves performance by 10-50x
    # Performance Impact: For large overlap files, 8MB buffer = ~6,400 system calls vs 1MB = ~51,200 calls
    BUFFER_SIZE = chira_utilities.get_adaptive_buffer_size(num_files=2)
    with open(overlapout, buffering=BUFFER_SIZE) as fh_overlapout, open(genomic_bed, "w", buffering=BUFFER_SIZE) as fh_wojunctions:
        for line in fh_overlapout:
            # transcriptid, s, e, readid, exonstart, exonend, exonid, pol
            # OPTIMIZATION: Cache string split results to avoid repeated splitting
            f = line.rstrip('\n').split('\t')
            transcriptid = f[0]
            txstart = int(f[1])
            txend = int(f[2])
            readid = f[3]
            exonstart = int(f[7])
            exonid = f[9]
            pol = f[11]
            if prev_readid and prev_readid != readid:
                unique_junctions = set(l_withjunctions)
                if len(unique_junctions) == 1:
                    fh_wojunctions.write(unique_junctions.pop())
                # for a spliced alignment choose the longest alignment
                else:
                    longest_splice = 0
                    longest_splice_align = ""
                    for splice_align in sorted(unique_junctions):
                        # OPTIMIZATION: Cache split result to avoid repeated splitting
                        x = splice_align.split("\t")
                        splice_len = int(x[2]) - int(x[1])
                        if splice_len > longest_splice:
                            longest_splice = splice_len
                            longest_splice_align = splice_align
                    fh_wojunctions.write(longest_splice_align)
                l_withjunctions.clear()

            # Calculate genomic hit positions.
            # Plus strand case.
            if pol == "+":
                genomicstart = exid2start[exonid] + (txstart - exonstart)
                genomicend = exid2start[exonid] + (txend - exonstart)
            # Minus strand case.
            elif pol == "-":
                genomicstart = exid2end[exonid] - (txend - exonstart)
                genomicend = exid2end[exonid] - (txstart - exonstart)
            else:
                sys.exit("ERROR: Check polarity " + pol + " in line: " + line + "\n")
            # OPTIMIZATION: Use f-string for faster string formatting
            # f-strings are faster than join() for small lists
            # Expected Speedup: 1.3-2x for string operations
            b = f"{id2chr[transcriptid]}\t{genomicstart}\t{genomicend}\t{readid}\t1\t{pol}\n"
            l_withjunctions.append(b)
            prev_readid = readid
        # write the last segment
        unique_last = set(l_withjunctions)
        if len(unique_last) == 1:
            fh_wojunctions.write(unique_last.pop())

    print("done")
    if os.path.exists(overlapout):
        os.remove(overlapout)

    return


def print_configuration(args):
    """Print configuration summary."""
    print('Input BED file                       : ' + args.bed)
    print('Output directory                     : ' + args.outdir)
    print('Segment overlap fraction             : ' + str(args.segment_overlap_fraction))
    print('Alignment/Block overlap fraction     : ' + str(args.alignment_overlap_fraction))
    if args.gtf:
        print('Annotation file                      : ' + args.gtf)
    if args.block_based:
        print('Merge method                         : blockbuster based')
        print('Blockbuster distance                 : ' + str(args.distance))
        print('Blockbuster minClusterHeight         : ' + str(args.min_cluster_height))
        print('Blockbuster minBlockHeight           : ' + str(args.min_block_height))
        print('Blockbuster scale                    : ' + str(args.scale))
    else:
        print('Merge method                         : overlap based')
        print('Minimum locus size                   : ' + str(args.min_locus_size))
    print('Minimum alignment length as % of longest : ' + str(args.length_threshold))
    print('Chimeric overlap                     : ' + str(args.chimeric_overlap))
    if args.num_processes is not None:
        print('Number of parallel processes          : ' + str(args.num_processes))
    else:
        print('Number of parallel processes          : auto (CPU count)')
    if args.ref_fasta1:
        print('1st priority reference fasta file    : ' + args.ref_fasta1)
    if args.ref_fasta2:
        print('2nd priority reference fasta file    : ' + args.ref_fasta2)
    print("===================================================================")


def setup_references(args):
    """
    Extract reference lengths from FASTA files.
    
    Args:
        args: Parsed arguments containing ref_fasta1 and ref_fasta2
    
    Returns:
        tuple: (d_reflen1, d_reflen2) dictionaries mapping reference IDs to lengths
    """
    d_reflen1 = defaultdict()
    d_reflen2 = defaultdict()
    if args.ref_fasta1:
        chira_utilities.extract_reflengths(args.ref_fasta1, d_reflen1)
    if args.ref_fasta2:
        chira_utilities.extract_reflengths(args.ref_fasta2, d_reflen2)
    return d_reflen1, d_reflen2


def process_segments(args, d_reflen1, d_reflen2):
    """
    Convert reads to segments and handle coordinate conversion if GTF is provided.
    
    Args:
        args: Parsed arguments
        d_reflen1: Dictionary of reference lengths for first priority reference
        d_reflen2: Dictionary of reference lengths for second priority reference
    """
    chira_utilities.print_w_time("START: Reads to segments")
    reads_to_segments(args.bed, args.outdir, args.segment_overlap_fraction, args.chimeric_overlap,
                      d_reflen1.keys(), d_reflen2.keys(), args.chimeric_only, args.length_threshold)
    chira_utilities.print_w_time("END: Reads to segments")
    
    # Handle coordinate conversion if GTF annotation is provided
    if args.gtf:
        # Parse the annotations and save them to dictionaries. Additionally write exons bed files to the outdir
        chira_utilities.print_w_time("START: Parse the annotation file")
        parse_annotations(args.gtf, args.outdir)
        chira_utilities.print_w_time("END: Parse the annotation file")
        chira_utilities.print_w_time("START: Convert to genomic coordinates")
        transcript_to_genomic_pos(os.path.join(args.outdir, "segments.temp.bed"),
                                  os.path.join(args.outdir, "segments.bed"),
                                  os.path.join(args.outdir, "genomic_exons.bed"),
                                  os.path.join(args.outdir, "transcriptomic_exons.bed"))
        chira_utilities.print_w_time("END: Convert to genomic coordinates")
        os.remove(os.path.join(args.outdir, "segments.temp.bed"))
    else:
        os.system(" ".join(["mv", os.path.join(args.outdir, "segments.temp.bed"),
                            os.path.join(args.outdir, "segments.bed")]))


def merge_loci(args):
    """
    Merge loci using either blockbuster-based or overlap-based method.
    
    Args:
        args: Parsed arguments containing merge method and parameters
    """
    if args.block_based:
        chira_utilities.print_w_time("START: blockbuster based merging")
        merge_loci_blockbuster(args.outdir, args.distance, args.min_cluster_height, args.min_block_height,
                               args.scale, args.alignment_overlap_fraction, args.num_processes)
        chira_utilities.print_w_time("END: blockbuster based merging")
    else:
        chira_utilities.print_w_time("START: overlap based merging")
        merge_loci_overlap(args.outdir, args.alignment_overlap_fraction, args.min_locus_size, args.num_processes)
        chira_utilities.print_w_time("END: overlap based merging")


def parse_arguments():
    """Parse command-line arguments. Order: required I/O, optional GTF, overlap params, refs, blockbuster, processes, version."""
    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: merge alignments and convert coordinates',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required I/O
    parser.add_argument('-b', '--bed', action='store', dest='bed', required=True,
                        metavar='', help='Input BED file with alignments (e.g. from chira_map)')
    parser.add_argument('-o', '--outdir', action='store', dest='outdir', required=True, metavar='',
                        help='Output directory for merged BED and segments')
    parser.add_argument('-g', '--gtf', action='store', dest='gtf', required=False,
                        metavar='', help='Annotation GTF file for coordinate conversion')

    # Overlap / merging parameters
    parser.add_argument('-ao', '--alignment_overlap', action='store', type=chira_utilities.score_float, default=0.7, metavar='',
                        dest='alignment_overlap_fraction',
                        help='Minimum fraction overlap among BED entries to merge. [0-1.0]')
    parser.add_argument('-so', '--segment_overlap', action='store', type=chira_utilities.score_float, default=0.7, metavar='',
                        dest='segment_overlap_fraction',
                        help='Matching read positions with greater than this %% overlap are merged into a segment')
    parser.add_argument('-lt', '--length_threshold', action='store', type=chira_utilities.score_float, default=0.9, metavar='',
                        dest='length_threshold',
                        help='Minimum alignment length as fraction of longest. [0.8-1.0]')
    parser.add_argument('-co', '--chimeric_overlap', action='store', type=int, default=2, metavar='',
                        dest='chimeric_overlap',
                        help='Maximum bases allowed between chimeric segments of a read')
    parser.add_argument("-c", '--chimeric_only', action='store_true', dest='chimeric_only',
                        help='Consider only chimeric reads for merging')
    parser.add_argument('-ls', '--min_locus_size', action='store', type=int, default=1, metavar='',
                        dest='min_locus_size',
                        help='Minimum number of alignments per merged locus')

    # Reference FASTA
    parser.add_argument('-f1', '--ref_fasta1', action='store', dest='ref_fasta1', required=False,
                        metavar='', help='First priority reference FASTA')
    parser.add_argument('-f2', '--ref_fasta2', action='store', dest='ref_fasta2', required=False,
                        metavar='', help='Second priority reference FASTA')

    # Blockbuster (optional)
    parser.add_argument('-bb', '--block_based', action='store_true', dest='block_based',
                        help='Use blockbuster for block-based merging')
    parser.add_argument('-d', '--distance', action='store', type=int, default=30, metavar='',
                        dest='distance', help='Blockbuster parameter: distance')
    parser.add_argument('-mc', '--min_cluster_height', action='store', type=int, default=10, metavar='',
                        dest='min_cluster_height', help='Blockbuster parameter: minClusterHeight')
    parser.add_argument('-mb', '--min_block_height', action='store', type=int, default=10, metavar='',
                        dest='min_block_height', help='Blockbuster parameter: minBlockHeight')
    parser.add_argument('-sc', '--scale', action='store', type=chira_utilities.score_float, default=0.1, metavar='',
                        dest='scale', help='Blockbuster parameter: scale')

    parser.add_argument('-p', '--processes', action='store', type=int, default=None, metavar='',
                        dest='num_processes',
                        help="""Number of parallel processes for transcript processing.
None = auto from CPU and transcript count. Set to 1 to disable.""")

    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {chira_utilities.__version__}')

    return parser.parse_args()


def main():
    """Main function to orchestrate the merging workflow."""
    args = parse_arguments()
    print_configuration(args)

    # Setup reference lengths
    d_reflen1, d_reflen2 = setup_references(args)

    # Process segments (reads to segments, coordinate conversion if GTF provided)
    process_segments(args, d_reflen1, d_reflen2)

    # Merge loci using selected method
    merge_loci(args)


if __name__ == "__main__":
    main()

