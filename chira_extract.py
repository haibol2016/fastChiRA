#!/usr/bin/env python
import chira_utilities
import argparse
import os
import sys
from collections import defaultdict
import itertools
import datetime
from BCBio import GFF
import subprocess
import multiprocessing
import shutil
import math
import json
import re
import time
import gzip
import warnings
import pickle

# OPTIMIZATION: Try to import intervaltree for efficient interval queries
# intervaltree provides O(log n) lookup instead of O(n) linear search
try:
    from intervaltree import IntervalTree
    INTERVALTREE_AVAILABLE = True
except ImportError:
    INTERVALTREE_AVAILABLE = False
    IntervalTree = None

try:
    from Bio import BiopythonDeprecationWarning
    warnings.filterwarnings("ignore", message=".*UnknownSeq.*deprecated.*", category=BiopythonDeprecationWarning)
except ImportError:
    # Fallback for older Biopython versions that don't have BiopythonDeprecationWarning
    warnings.filterwarnings("ignore", message=".*UnknownSeq.*deprecated.*", category=DeprecationWarning)

# MPIRE is REQUIRED for multiprocessing performance
# MPIRE provides shared objects, lower overhead, and better performance than Process
# Benefits: 50-90% memory reduction, 2-3x faster startup, better performance
# Install with: pip install mpire (or pip install chira, or conda install -c conda-forge mpire)
# Note: Requires MPIRE >= 2.4.0 for shared_objects support
# Note: MPIRE 2.10.1+ doesn't use a SharedObject class
# Shared objects are passed directly to WorkerPool via shared_objects parameter
# See: https://sybrenjansen.github.io/mpire/v2.10.1/usage/workerpool/shared_objects.html

from mpire import WorkerPool



d_gene_annotations = defaultdict(lambda: defaultdict(str))
d_transcript_annotations = defaultdict(lambda: defaultdict(list))
# OPTIMIZATION: Store interval trees for UTR/CDS features per transcript
# This enables O(log n) lookup instead of O(n) linear search in guess_region()
d_transcript_interval_trees = defaultdict(dict)  # {transcript_id: {'UTR': IntervalTree, 'CDS': IntervalTree}}

# Constants for chimera list indices
# *.chimeras.n and *.chimeras.txt are tab-separated (no header on .n chunks).
# Columns 0-19: read_id, transcript/gene/region info, genomic_coordinates_1/2; 20-21: locus_id_1, locus_id_2 (e.g. chr11:73780216:73780507:+); 22-28: crl_group_id, tpm, scores; 29-34: hybridized_sequences, structure, positions, mfe, mirna_position, hybridized_subsequences.
CHIMERA_IDX_LOCUS1 = 20
CHIMERA_IDX_LOCUS2 = 21
CHIMERA_IDX_SEQUENCES = 29
CHIMERA_IDX_HYBRID = 30
CHIMERA_IDX_HYBRID_POS = 31
CHIMERA_IDX_MFE = 32
CHIMERA_IDX_MIRNA_POSITION = 33
CHIMERA_IDX_HYBRID_SUBSEQS = 34

# miRNA region types
MIRNA_REGION_TYPES = ["miRNA", "3p_mature_mir", "5p_mature_mir", "mature_mir"]

# I/O buffer size (2MB) used for GTF, CRL, FASTA, chimera and merge files
BUFFER_SIZE = 2 * 1024 * 1024


# --- Region/strand and GTF annotation (used by extraction) ---

def strandardize(strand):
    if strand == '-1':
        strand = '-'
    elif strand == '1':
        strand = '+'
    return strand



def guess_region(transcriptid, read_pos):
    [read_chr, read_start, read_end, read_strand] = read_pos.split(':')[-4:]
    region = 'NA'
    overlap_length = 0
    utr_start = utr_end = first_cds_start = last_cds_end = strand = None
    read_strand = strandardize(read_strand)
    read_start_int = int(read_start)
    read_end_int = int(read_end)
    
    # OPTIMIZATION: Use interval trees for O(log n) lookup instead of O(n) linear search
    # Falls back to linear search if intervaltree is not available
    if INTERVALTREE_AVAILABLE and transcriptid in d_transcript_interval_trees:
        # Use interval tree for UTR lookup
        if 'UTR' in d_transcript_interval_trees[transcriptid]:
            utr_tree = d_transcript_interval_trees[transcriptid]['UTR']
            # Query interval tree for overlapping UTR features
            # IntervalTree uses [start, end) (end is exclusive), so we use read_end_int + 1
            # Query returns intervals that overlap with [read_start_int, read_end_int + 1)
            overlapping_utrs = utr_tree[read_start_int:read_end_int + 1]
            for interval in overlapping_utrs:
                # interval.data contains: (chrom, strand, utr_type, original_pos_list)
                data = interval.data
                if not data or len(data) < 3:
                    continue
                chrom, feat_strand, utr_type = data[0], data[1], data[2]
                if read_chr != chrom or read_strand != feat_strand:
                    continue
                # Calculate overlap with actual read coordinates
                # interval.begin and interval.end are already in the correct format
                feat_start = interval.begin
                feat_end = interval.end - 1  # Convert back to inclusive end (interval tree uses exclusive end)
                overlap = chira_utilities.overlap([feat_start, feat_end], [read_start_int, read_end_int])
                if overlap > overlap_length:
                    overlap_length = overlap
                    if utr_type == 'five_prime_utr':
                        region = '5_prime_UTR'
                    elif utr_type == 'three_prime_utr':
                        region = '3_prime_UTR'
                    else:
                        region = 'UTR'
                    utr_start = feat_start
                    utr_end = feat_end
                    # Set strand from UTR features (needed for UTR 5'/3' determination if no CDS)
                    # All features of a transcript have the same strand
                    if strand is None:
                        strand = feat_strand
        
        # Use interval tree for CDS lookup
        if 'CDS' in d_transcript_interval_trees[transcriptid]:
            cds_tree = d_transcript_interval_trees[transcriptid]['CDS']
            overlapping_cdss = cds_tree[read_start_int:read_end_int + 1]
            for interval in overlapping_cdss:
                data = interval.data
                if not data or len(data) < 2:
                    continue
                chrom, feat_strand = data[0], data[1]
                if read_chr != chrom or read_strand != feat_strand:
                    continue
                feat_start = interval.begin
                feat_end = interval.end - 1  # Convert back to inclusive end
                overlap = chira_utilities.overlap([feat_start, feat_end], [read_start_int, read_end_int])
                if overlap > overlap_length:
                    overlap_length = overlap
                    region = 'CDS'
                # Track CDS boundaries for UTR determination
                if not first_cds_start or feat_start < first_cds_start:
                    first_cds_start = feat_start
                if not last_cds_end or feat_end > last_cds_end:
                    last_cds_end = feat_end
                # Set strand from CDS features (needed for UTR 5'/3' determination)
                # All features of a transcript have the same strand
                if strand is None:
                    strand = feat_strand
    else:
        # Fallback to linear search if intervaltree is not available
        if transcriptid in d_transcript_annotations['UTR']:
            for pos in d_transcript_annotations['UTR'][transcriptid]:
                chrom = pos[0]
                start = int(pos[1])
                end = int(pos[2])
                strand = strandardize(pos[3])
                # UTR type is stored as 5th element (if available), otherwise defaults to 'UTR'
                utr_type = pos[4] if len(pos) > 4 else 'UTR'
                if read_chr != chrom or read_strand != strand:
                    continue

                if chira_utilities.overlap([start, end], [read_start_int, read_end_int]) > overlap_length:
                    overlap_length = chira_utilities.overlap([start, end], [read_start_int, read_end_int])
                    # Use specific UTR type if available, otherwise use generic 'UTR'
                    if utr_type == 'five_prime_utr':
                        region = '5_prime_UTR'
                    elif utr_type == 'three_prime_utr':
                        region = '3_prime_UTR'
                    else:
                        region = 'UTR'
                    utr_start = start
                    utr_end = end

        if transcriptid in d_transcript_annotations['CDS']:
            for pos in d_transcript_annotations['CDS'][transcriptid]:
                chrom = pos[0]
                start = int(pos[1])
                end = int(pos[2])
                strand = strandardize(pos[3])
                if read_chr != chrom or read_strand != strand:
                    continue
                if chira_utilities.overlap([start, end], [read_start_int, read_end_int]) > overlap_length:
                    overlap_length = chira_utilities.overlap([start, end], [read_start_int, read_end_int])
                    region = 'CDS'
                if not first_cds_start or start < first_cds_start:
                    first_cds_start = start
                if not last_cds_end or end > last_cds_end:
                    last_cds_end = end

    # if region is still a generic UTR (not already determined from feature type)
    # determine 5' vs 3' based on position relative to CDS
    # Note: strand should be set from CDS features (all features of a transcript have the same strand)
    if region == 'UTR' and first_cds_start is not None and last_cds_end is not None and strand is not None:
        if utr_end <= first_cds_start:
            region = '5_prime_UTR'
            if strand == '-':
                region = '3_prime_UTR'
        elif utr_start >= last_cds_end:
            region = '3_prime_UTR'
            if strand == '-':
                region = '5_prime_UTR'

    # works for mirbase gff3
    if transcriptid in d_transcript_annotations['gid']:
        geneid = d_transcript_annotations['gid'][transcriptid]
        if geneid in d_transcript_annotations['mature']:
            for pos in d_transcript_annotations['mature'][geneid]:
                chrom = pos[0]
                start = int(pos[1])
                end = int(pos[2])
                strand = strandardize(pos[3])
                name = pos[4]
                if read_chr != chrom or read_strand != strand:
                    continue
                if chira_utilities.overlap([start, end], [int(read_start), int(read_end)]) > overlap_length:
                    overlap_length = chira_utilities.overlap([start, end], [int(read_start), int(read_end)])
                    if name.endswith("-3p"):
                        region = "3p_mature_mir"
                    elif name.endswith("-5p"):
                        region = "5p_mature_mir"
                    else:
                        region = "mature_mir"
    return region



def build_interval_trees():
    """
    Build interval trees for UTR and CDS features per transcript.
    
    OPTIMIZATION: This enables O(log n) lookup instead of O(n) linear search
    in guess_region(). Interval trees are built once after parsing annotations,
    then reused for all region lookups.
    """
    # Build UTR interval trees
    for transcript_id, utr_list in d_transcript_annotations['UTR'].items():
        utr_tree = IntervalTree()
        for pos in utr_list:
            chrom = pos[0]
            start = int(pos[1])
            end = int(pos[2]) + 1  # IntervalTree uses [start, end) (end is exclusive)
            strand = strandardize(pos[3])
            utr_type = pos[4] if len(pos) > 4 else 'UTR'
            # Store chrom, strand, utr_type, and original pos for later use
            # interval.data can store any data associated with the interval
            utr_tree[start:end] = (chrom, strand, utr_type, pos)
        if utr_tree:
            d_transcript_interval_trees[transcript_id]['UTR'] = utr_tree
    
    # Build CDS interval trees
    for transcript_id, cds_list in d_transcript_annotations['CDS'].items():
        cds_tree = IntervalTree()
        for pos in cds_list:
            chrom = pos[0]
            start = int(pos[1])
            end = int(pos[2]) + 1  # IntervalTree uses [start, end) (end is exclusive)
            strand = strandardize(pos[3])
            # Store chrom, strand, and original pos for later use
            cds_tree[start:end] = (chrom, strand, pos)
        if cds_tree:
            d_transcript_interval_trees[transcript_id]['CDS'] = cds_tree



def parse_annotations(f_gtf):
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
    d_attributes['name'] = ['gene_name', 'Name']
    # Use transcript_biotype first (transcript-level), then gene_biotype as fallback
    d_attributes['type'] = ['transcript_biotype', 'gene_biotype', 'Type']

    l_seen_exons = set()
    with open(f_gtf, "r", buffering=BUFFER_SIZE) as gtf_handle:
        # Chromosome seq level
        for rec in GFF.parse(gtf_handle, limit_info=limit_info, target_lines=1):
            # for each selected sub_feature
            for sub_feature in rec.features:

                # for some reason each qualifier is a list! take the first element
                for i in d_attributes['gid']:
                    if i in sub_feature.qualifiers:
                        gene_id = sub_feature.qualifiers[i][0]
                        break
                for i in d_attributes['tid']:
                    if i in sub_feature.qualifiers:
                        transcript_id = sub_feature.qualifiers[i][0]
                        break

                # NOTE: some mature miRs have multiple locations on genome
                # this is to select only the first location for mature miR

                # Handle UTR features: support both generic "UTR" and Ensembl-specific types
                if sub_feature.type in ['UTR', 'five_prime_utr', 'three_prime_utr']:
                    if transcript_id not in d_transcript_annotations['UTR']:
                        d_transcript_annotations['UTR'][transcript_id] = []
                    # Store UTR type information for later use in region assignment
                    utr_type = sub_feature.type
                    d_transcript_annotations['UTR'][transcript_id].append([rec.id,
                                                                           str(sub_feature.location.start),
                                                                           str(sub_feature.location.end),
                                                                           str(sub_feature.location.strand),
                                                                           utr_type])
                elif sub_feature.type == 'CDS':
                    if transcript_id not in d_transcript_annotations['CDS']:
                        d_transcript_annotations['CDS'][transcript_id] = []
                    d_transcript_annotations['CDS'][transcript_id].append([rec.id,
                                                                           str(sub_feature.location.start),
                                                                           str(sub_feature.location.end),
                                                                           str(sub_feature.location.strand)])
                # remaining are exon, miRNA and tRNA lines
                else:
                    if transcript_id + "_e" + str(n_exon).zfill(3) in l_seen_exons:
                        continue

                    # reset the variables if it is a new transcript
                    if prev_transcript_id != transcript_id:
                        d_transcript_annotations['len'][transcript_id] = 0
                        n_exon = 1
                        exon_rel_start = 0
                        exon_rel_end = 0
                        exon_len = 0

                    d_transcript_annotations['gid'][transcript_id] = gene_id
                    # biotype
                    biotype = 'NA'
                    for i in d_attributes['type']:
                        if i in sub_feature.qualifiers:
                            biotype = sub_feature.qualifiers[i][0]
                    if sub_feature.type == "miRNA":
                        biotype = "miRNA"
                    if sub_feature.type == "tRNA":
                        biotype = "tRNA"
                    d_gene_annotations['type'][gene_id] = biotype

                    # name
                    try:
                        for i in d_attributes['name']:
                            if i in sub_feature.qualifiers:
                                d_gene_annotations['name'][gene_id] = sub_feature.qualifiers[i][0]
                    except KeyError:
                        d_gene_annotations['name'][gene_id] = "NA"

                    exon_len = sub_feature.location.end - sub_feature.location.start
                    exon_rel_end = exon_rel_start + exon_len
                    # TODO: check +1 to length ?
                    d_transcript_annotations['len'][transcript_id] += exon_len + 1
                    exon_rel_start = exon_rel_end
                    prev_transcript_id = transcript_id
                    l_seen_exons.add(transcript_id + "_e" + str(n_exon).zfill(3))
                    n_exon += 1
    
    # OPTIMIZATION: Build interval trees for UTR/CDS features after parsing annotations
    # This enables O(log n) lookup instead of O(n) linear search in guess_region()
    if INTERVALTREE_AVAILABLE:
        build_interval_trees()




# --- Alignment parsing and locus helpers ---

def parse_alignment_line(alignment_line):
    """
    Parse an alignment line once into a structured dictionary.
    
    OPTIMIZATION: This avoids repeated string splitting operations and float conversions.
    Parse once, reuse the parsed structure throughout the function.
    Converts numeric fields to appropriate types (float, int) to avoid repeated conversions.
    """
    fields = alignment_line.rstrip('\n').split('\t')
    if len(fields) < 13:
        return None
    # BUGFIX: Field indices must match correct loci.counts header format
    # Column 10: locus_share (fraction of reads in CRL that map to this locus)
    # Column 11: read_crl_fraction (EM-optimized probability that read belongs to this CRL)
    # Column 12: crl_tpm (TPM value for the CRL)
    # OPTIMIZATION: Convert numeric fields once to avoid repeated float() calls
    # Cache float values to avoid repeated conversions in extract_and_write()
    try:
        locus_share_float = float(fields[10])      # Field 10: locus_share
        read_crl_fraction_float = float(fields[11]) # Field 11: read_crl_fraction (EM-optimized probability)
        tpm_float = float(fields[12])               # Field 12: crl_tpm
    except (ValueError, IndexError):
        # Fallback: set to 0.0 if conversion fails
        locus_share_float = 0.0
        read_crl_fraction_float = 0.0
        tpm_float = 0.0
    
    return {
        'segmentid': fields[0],
        'transcriptid': fields[1],
        'locusid': fields[2],
        'crlid': fields[3],
        'tx_pos_start': fields[4],
        'tx_pos_end': fields[5],
        'tx_pos_strand': fields[6],
        'cigar': fields[7],
        'genomic_pos': fields[8],
        'locuspos': fields[9],
        'locus_share': fields[10],  # Field 10: locus_share
        'locus_share_float': locus_share_float,  # Cached float value
        'read_crl_fraction': fields[11],  # Field 11: read_crl_fraction (EM-optimized probability)
        'read_crl_fraction_float': read_crl_fraction_float,  # Cached float value
        # Keep 'prob' and 'prob_float' as aliases for backward compatibility in code
        'prob': fields[11],  # Alias: read_crl_fraction is the probability
        'prob_float': read_crl_fraction_float,  # Alias: cached float value
        # Keep 'locusshare' and 'locusshare_float' as aliases for backward compatibility
        'locusshare': fields[10],  # Alias: locus_share
        'locusshare_float': locus_share_float,  # Alias: cached float value
        'tpm': fields[12],  # Field 12: crl_tpm
        'tpm_float': tpm_float,  # Cached float value
        'raw_line': alignment_line  # Keep original for compatibility if needed
    }



def filter_alignments(lines, tpm_threshold, score_cutoff):
    l_read_aligns = []
    for line in lines:
        f = line.rstrip('\n').split('\t')
        crl_tpm = f[12]
        # Correct field indices based on loci.counts header:
        # Column 10: locus_share (fraction of reads in CRL that map to this locus)
        # Column 11: read_crl_fraction (EM-optimized probability that read belongs to this CRL)
        # Column 12: crl_tpm (TPM value for the CRL)
        locus_share = f[10]
        read_crl_fraction = f[11]
        if float(crl_tpm) < tpm_threshold:
            continue
        locus_score = float(read_crl_fraction) * float(locus_share)
        if locus_score < score_cutoff:
            continue
        l_read_aligns.append(line)
    return l_read_aligns



def add_locus_to_set(locus, l_loci_bed):
    pos = locus.split(":")
    bed = "\t".join([":".join(pos[0:-3]), pos[-3], pos[-2], locus, "1", pos[-1]])
    if bed not in l_loci_bed:
        l_loci_bed.add(bed)



def update_best_hits(l_best_hits, hit_type):
    best_tpm = 0
    for i, hit in enumerate(l_best_hits):
        if hit_type == "chimera":
            tpm = float(hit[24]) + float(hit[25])
        elif hit_type == "singleton":
            tpm = float(hit[13])
        if tpm > best_tpm:
            best_tpm = tpm
    l_best_hits = [n for n in l_best_hits if n is not None]
    for i, hit in enumerate(l_best_hits):
        if hit_type == "chimera":
            tpm = float(hit[24]) + float(hit[25])
        elif hit_type == "singleton":
            tpm = float(hit[13])
        if tpm < best_tpm:
            l_best_hits[i] = None
    l_best_hits = [n for n in l_best_hits if n is not None]
    return l_best_hits



def extract_annotations(transcriptid, genomic_pos, d_regions, f_gtf):
    # OPTIMIZATION: Use .get() method for dictionary lookups
    # More Pythonic and slightly faster than conditional checks
    # .get() with default value avoids KeyError and is more efficient
    geneid = d_transcript_annotations['gid'].get(transcriptid, 'NA')
    name = d_gene_annotations['name'].get(geneid, 'NA')
    biotype = d_gene_annotations['type'].get(geneid, 'NA')
    if biotype == 'miRNA' or biotype == 'tRNA':
        # OPTIMIZATION: Use .get() with default to avoid KeyError check
        name = d_gene_annotations['family'].get(geneid, name)

    # OPTIMIZATION: Cache region lookups in d_regions dictionary
    # This avoids repeated calls to guess_region() for the same (transcriptid, genomic_pos) combination
    # Key format: transcriptid + '\t' + genomic_pos
    cache_key = transcriptid + '\t' + genomic_pos
    if cache_key in d_regions:
        region = d_regions[cache_key]
    else:
        region = "NA"
        if f_gtf:
            region = guess_region(transcriptid, genomic_pos)
        if region == "NA":
            region = biotype
        # Cache the result for future lookups
        d_regions[cache_key] = region

    # OPTIMIZATION: Use .get() method for dictionary lookup
    tx_length = d_transcript_annotations['len'].get(transcriptid, 'NA')
    return geneid, name, region, tx_length



def extract_and_write(readid, l_read_alignments, l_loci_bed, d_ref_lengths1, d_ref_lengths2, f_gtf,
                      d_regions, chimeric_overlap, fh_chimeras, fh_singletons, hybridize):
    chimera_found = False
    l_best_chimeras = []
    l_best_singletons = []
    
    # OPTIMIZATION: Parse ALL alignments once upfront into structured dictionaries
    # This avoids repeated string splitting operations (rstrip, split) throughout the function
    # Pre-strip and split all alignments at the start, then cache parsed results for reuse
    # Benefits: Each alignment is parsed exactly once, regardless of how many times it's used
    parsed_alignments = []  # List of all parsed alignments (for filtering)
    alignment_to_parsed = {}  # Map original alignment string to parsed dict (for lookups)
    
    # Parse all alignments once at the start
    for alignment in l_read_alignments:
        parsed = parse_alignment_line(alignment)
        if parsed is not None:
            parsed_alignments.append(parsed)
            alignment_to_parsed[alignment] = parsed
    
    # OPTIMIZATION: Pre-filter alignments before generating pairs to reduce O(n²) combinations
    # Filter out alignments with read_crl_fraction == 0 early, as they will never form valid chimeras
    # This dramatically reduces the number of pairs generated, especially for reads with many alignments
    # Example: 20 alignments with 5 having read_crl_fraction==0: 190 pairs -> 105 pairs (45% reduction)
    # Use parsed structures for filtering (faster than re-parsing)
    filtered_alignments = []
    for parsed in parsed_alignments:
        # Only keep alignments with non-zero read_crl_fraction (EM-optimized probability)
        # Use cached float value if available, otherwise parse once
        prob_float = parsed.get('prob_float', float(parsed['prob']))  # prob_float is alias for read_crl_fraction_float
        if prob_float != 0:
            # Get original alignment string from parsed dict (stored as 'raw_line')
            filtered_alignments.append(parsed['raw_line'])
    
    # OPTIMIZATION: Use generator instead of list to avoid creating full list in memory
    # This reduces memory usage, especially for reads with many alignments (O(n²) pairs)
    # Generator is lazy-evaluated, so pairs are created on-demand during iteration
    alignment_pairs = itertools.combinations(filtered_alignments, 2)

    for alignment1, alignment2 in alignment_pairs:
        # OPTIMIZATION: Use cached parsed structures instead of re-splitting
        # This avoids repeated rstrip() and split() operations on the same strings
        p1 = alignment_to_parsed[alignment1]
        p2 = alignment_to_parsed[alignment2]
        
        segmentid1 = p1['segmentid']
        transcriptid1 = p1['transcriptid']
        locusid1 = p1['locusid']
        crlid1 = p1['crlid']
        tx_pos_start1 = p1['tx_pos_start']
        tx_pos_end1 = p1['tx_pos_end']
        tx_pos_strand1 = p1['tx_pos_strand']
        cigar1 = p1['cigar']
        genomic_pos1 = p1['genomic_pos']
        locuspos1 = p1['locuspos']
        locusshare1 = p1['locusshare']
        # prob1, tpm1 not needed - we use prob1_float, tpm1_float instead
        
        segmentid2 = p2['segmentid']
        transcriptid2 = p2['transcriptid']
        locusid2 = p2['locusid']
        crlid2 = p2['crlid']
        tx_pos_start2 = p2['tx_pos_start']
        tx_pos_end2 = p2['tx_pos_end']
        tx_pos_strand2 = p2['tx_pos_strand']
        cigar2 = p2['cigar']
        genomic_pos2 = p2['genomic_pos']
        locuspos2 = p2['locuspos']
        locusshare2 = p2['locusshare']
        # prob2, tpm2 not needed - we use prob2_float, tpm2_float instead
        
        # these are multimappings of the same segment
        if segmentid1 == segmentid2 or locuspos1 == locuspos2 or crlid1 == crlid2:
            continue
        # check these are chimeric arms
        if not chira_utilities.is_chimeric(cigar1, cigar2, tx_pos_strand1 == "-",
                                           tx_pos_strand2 == "-", chimeric_overlap):
            continue
        switch_alignments = False
        # if the second reference is provided then it is assumed to be a split reference
        if len(d_ref_lengths2) != 0:
            # check if both transcripts are from the same reference database
            if transcriptid1 in d_ref_lengths1 and transcriptid2 in d_ref_lengths1 or \
                    transcriptid1 in d_ref_lengths2 and transcriptid2 in d_ref_lengths2:
                continue
            # switch orientation of the alignments based on references
            if transcriptid1 in d_ref_lengths2 or transcriptid2 in d_ref_lengths1:
                switch_alignments = True
        # not a split reference
        else:
            # then switch alignments based on the reference id
            if transcriptid2 > transcriptid1:
                switch_alignments = True
            elif locuspos2 > locuspos1:
                switch_alignments = True
        if switch_alignments:
            # OPTIMIZATION: Swap parsed structures instead of re-parsing
            # This avoids re-splitting the alignment strings
            p1, p2 = p2, p1
            segmentid1 = p1['segmentid']
            transcriptid1 = p1['transcriptid']
            locusid1 = p1['locusid']
            crlid1 = p1['crlid']
            tx_pos_start1 = p1['tx_pos_start']
            tx_pos_end1 = p1['tx_pos_end']
            tx_pos_strand1 = p1['tx_pos_strand']
            cigar1 = p1['cigar']
            genomic_pos1 = p1['genomic_pos']
            locuspos1 = p1['locuspos']
            locusshare1 = p1['locusshare']
            # prob1, tpm1 not needed - we use prob1_float, tpm1_float instead
            
            segmentid2 = p2['segmentid']
            transcriptid2 = p2['transcriptid']
            locusid2 = p2['locusid']
            crlid2 = p2['crlid']
            tx_pos_start2 = p2['tx_pos_start']
            tx_pos_end2 = p2['tx_pos_end']
            tx_pos_strand2 = p2['tx_pos_strand']
            cigar2 = p2['cigar']
            genomic_pos2 = p2['genomic_pos']
            locuspos2 = p2['locuspos']
            locusshare2 = p2['locusshare']
            # prob2, tpm2 not needed - we use prob2_float, tpm2_float instead
        # OPTIMIZATION: Use cached float values instead of repeated conversions
        # Parse once in parse_alignment_line(), reuse cached values here
        # Match original behavior: float("{:.2f}".format(float(prob))) 
        # Using round() produces equivalent results but is more efficient
        # Note: prob_float is alias for read_crl_fraction_float (EM-optimized probability)
        prob1_float = p1.get('prob_float', float(p1['prob']))  # read_crl_fraction for alignment 1
        prob2_float = p2.get('prob_float', float(p2['prob']))  # read_crl_fraction for alignment 2
        # Safety check: ensure read_crl_fraction != 0 (should already be filtered by pre-filtering, but defensive check)
        # Original code checked prob == 0 here; pre-filtering optimization moved it earlier for performance
        # This check ensures we don't process pairs with read_crl_fraction == 0 even if pre-filtering missed them
        if prob1_float == 0 or prob2_float == 0:
            continue
        # For chimeras, use read_crl_fraction directly as the score (not multiplied by locus_share)
        # Original: first_locus_score = float("{:.2f}".format(float(prob1)))
        # Equivalent: round to 2 decimal places (produces same result)
        first_locus_score = round(prob1_float, 2)  # read_crl_fraction for alignment 1
        second_locus_score = round(prob2_float, 2)  # read_crl_fraction for alignment 2
        combined_score = first_locus_score * second_locus_score
        
        # OPTIMIZATION: Use cached float values and f-strings for formatting
        tpm1_float = p1.get('tpm_float', float(p1['tpm']))
        tpm2_float = p2.get('tpm_float', float(p2['tpm']))
        tpm1 = f"{tpm1_float:.2f}"
        tpm2 = f"{tpm2_float:.2f}"

        geneid1, name1, region1, tx_len1 = extract_annotations(transcriptid1,
                                                               genomic_pos1,
                                                               d_regions,
                                                               f_gtf)
        geneid2, name2, region2, tx_len2 = extract_annotations(transcriptid2,
                                                               genomic_pos2,
                                                               d_regions,
                                                               f_gtf)

        arm1_start, arm1_end = chira_utilities.match_positions(cigar1, tx_pos_strand1 == "-")
        arm2_start, arm2_end = chira_utilities.match_positions(cigar2, tx_pos_strand2 == "-")
        read_length = chira_utilities.query_length(cigar1, tx_pos_strand1 == "-")

        if tx_len1 == "NA" or tx_len2 == "NA":
            if len(d_ref_lengths2) == 0:
                tx_len1 = d_ref_lengths1[transcriptid1]
                tx_len2 = d_ref_lengths1[transcriptid2]
            else:
                try:
                    tx_len1 = d_ref_lengths1[transcriptid1]
                    tx_len2 = d_ref_lengths2[transcriptid2]
                except KeyError:
                    tx_len1 = d_ref_lengths2[transcriptid1]
                    tx_len2 = d_ref_lengths1[transcriptid2]

        # Determine if miRNA is first or last in the read
        # Check if region1 or region2 is miRNA
        is_mirna1 = region1 in MIRNA_REGION_TYPES
        is_mirna2 = region2 in MIRNA_REGION_TYPES
        
        # Determine read orientation based on arm positions
        # arm1 is at 5' end if arm1_start < arm2_start, otherwise arm2 is at 5' end
        if (arm1_start < arm2_start and is_mirna1) or (arm1_start >= arm2_start and is_mirna2):
            mirna_position = "miRNA_first"
        elif (arm1_start < arm2_start and is_mirna2) or (arm1_start >= arm2_start and is_mirna1):
            mirna_position = "miRNA_last"
        else:
            mirna_position = "NA"  # Neither is miRNA (shouldn't happen in typical split reference)
        
        # OPTIMIZATION: Use f-strings for faster string formatting
        # Pre-format numeric values and use f-strings instead of str() and join()
        # All fields are already strings from parsed dict, except numeric values which we format
        alignment_info = f"{arm1_start},{arm1_end},{arm2_start},{arm2_end},{read_length}"
        chimera = [readid, transcriptid1, transcriptid2,
                   geneid1, geneid2, name1, name2, region1, region2,
                   tx_pos_start1, tx_pos_end1, tx_pos_strand1, str(tx_len1),
                   tx_pos_start2, tx_pos_end2, tx_pos_strand2, str(tx_len2),
                   alignment_info,
                   genomic_pos1, genomic_pos2,
                   locuspos1, locuspos2,
                   crlid1, crlid2,
                   tpm1, tpm2,
                   f"{first_locus_score:.2f}", f"{second_locus_score:.2f}", f"{combined_score:.2f}",
                   "NA",  # sequences
                   "NA",  # hybrid
                   "NA",  # hybrid_pos
                   "NA",  # mfe
                   mirna_position]
        chimera_found = True
        l_best_chimeras.append(chimera)

    # if there are no pairs, then it is a singleton
    if not chimera_found:
        # singleton read
        # OPTIMIZATION: Use cached parsed structures for singletons
        # All alignments are already parsed upfront, so just iterate through parsed_alignments
        for p in parsed_alignments:
            
            segmentid = p['segmentid']
            transcriptid = p['transcriptid']
            locusid = p['locusid']
            crlid = p['crlid']
            tx_pos_start = p['tx_pos_start']
            tx_pos_end = p['tx_pos_end']
            tx_pos_strand = p['tx_pos_strand']
            cigar = p['cigar']
            genomic_pos = p['genomic_pos']
            locuspos = p['locuspos']
            locusshare = p['locusshare']
            # prob, tpm not needed - we use prob_float, tpm_float instead

            geneid, name, region, tx_len = extract_annotations(transcriptid,
                                                               genomic_pos,
                                                               d_regions,
                                                               f_gtf)

            # OPTIMIZATION: Use cached float values instead of repeated conversions
            # Note: prob_float is alias for read_crl_fraction_float, locusshare_float is alias for locus_share_float
            prob_float = p.get('prob_float', float(p['prob']))  # read_crl_fraction (EM-optimized probability)
            locusshare_float = p.get('locusshare_float', float(locusshare))  # locus_share
            tpm_float = p.get('tpm_float', float(p['tpm']))
            # For singletons, use read_crl_fraction * locus_share as the alignment score
            locus_score = f"{(prob_float * locusshare_float):.2f}"
            tpm = f"{tpm_float:.2f}"

            arm_start, arm_end = chira_utilities.match_positions(cigar, tx_pos_strand == "-")
            read_length = chira_utilities.query_length(cigar, tx_pos_strand == "-")

            # OPTIMIZATION: Use f-strings for faster string formatting
            alignment_info = f"{arm_start},{arm_end},{read_length}"
            singleton = [readid, transcriptid, geneid, name, region,
                         tx_pos_start, tx_pos_end, tx_pos_strand, str(tx_len),
                         alignment_info,
                         genomic_pos,
                         locuspos,
                         crlid,
                         tpm,
                         locus_score]
            l_best_singletons.append(singleton)

    # OPTIMIZATION: Batch writing for better I/O performance
    # Collect lines in buffer, write in batches (e.g., every 10,000 lines) using writelines()
    # This reduces system calls and improves performance for large datasets
    BATCH_SIZE = 10000
    
    if len(l_best_chimeras) > 0:
        l_best_chimeras = update_best_hits(l_best_chimeras, "chimera")
        output_buffer = []
        for a in l_best_chimeras:
            if hybridize:
                add_locus_to_set(a[CHIMERA_IDX_LOCUS1], l_loci_bed)
                add_locus_to_set(a[CHIMERA_IDX_LOCUS2], l_loci_bed)
            # OPTIMIZATION: Use f-string for faster string formatting
            output_buffer.append("\t".join(str(x) for x in a) + "\n")
            if len(output_buffer) >= BATCH_SIZE:
                fh_chimeras.writelines(output_buffer)
                output_buffer = []
        # Write remaining lines
        if output_buffer:
            fh_chimeras.writelines(output_buffer)
    else:
        l_best_singletons = update_best_hits(l_best_singletons, "singleton")
        output_buffer = []
        for b in l_best_singletons:
            # OPTIMIZATION: Use f-string for faster string formatting
            output_buffer.append("\t".join(str(x) for x in b) + "\n")
            if len(output_buffer) >= BATCH_SIZE:
                fh_singletons.writelines(output_buffer)
                output_buffer = []
        # Write remaining lines
        if output_buffer:
            fh_singletons.writelines(output_buffer)




# --- IntaRNA / hybridization ---

def parse_intarna_csv(csv_path):
    """
    Parse IntaRNA CSV (--outCsvCols ... E,seq1,seq2). Per IntaRNA README: index 1 = target, 2 = query;
    id1 = target id, id2 = query id; seq1/subseq1/start1 = target, seq2/subseq2/start2 = query;
    hybridDPfull = target_part & query_part (id1 & id2). We store dotbracket/pos in id1&id2 order.
    Returns dict (id1, id2) -> (dotbracket, pos, energy_str, subseq1, subseq2, seq_for_id1, seq_for_id2).
    Keeps one row per pair (lowest energy).
    """
    d_hybrids = {}
    if not os.path.exists(csv_path):
        sys.stderr.write(f"Warning: IntaRNA result file not found: {csv_path}\n")
        sys.stderr.write("  Ensure batchtools LSF jobs completed successfully and the path is accessible from this host.\n")
        return d_hybrids
    with open(csv_path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    if len(lines) < 2:
        sys.stderr.write(f"Warning: IntaRNA result file empty or header-only (got {len(lines)} line(s)): {csv_path}\n")
        return d_hybrids

    data_lines = [ln for ln in lines[1:] if ln.strip()]
    for line in data_lines:
        line = line.strip()
        if not line:
            continue
        parts = line.replace(",", ";").split(";")
        if len(parts) < 8:
            continue
        id1 = parts[0].strip()
        id2 = parts[1].strip()
        start1, subseq1 = parts[2].strip(), parts[3].strip()
        start2, subseq2 = parts[4].strip(), parts[5].strip()
        dotbracket = parts[6].strip()
        energy_str = parts[7].strip()
        # IntaRNA: seq1 = target (id1), seq2 = query (id2). Store (seq_for_id1, seq_for_id2).
        seq_for_id1 = parts[8].strip() if len(parts) > 8 else None
        seq_for_id2 = parts[9].strip() if len(parts) > 9 else None
        if seq_for_id1 == "":
            seq_for_id1 = None
        if seq_for_id2 == "":
            seq_for_id2 = None
        if "&" not in dotbracket:
            continue
        # hybridDPfull is target&query (id1&id2) per IntaRNA README; keep as-is (id1='(', id2=')')
        pos = start1 + "&" + start2  # id1 & id2
        val = (dotbracket, pos, energy_str, subseq1, subseq2, seq_for_id1, seq_for_id2)
        try:
            energy_float = float(energy_str)
        except (ValueError, TypeError):
            energy_float = 0.0
        key = (id1, id2)
        if key not in d_hybrids:
            d_hybrids[key] = (val, energy_float)
        else:
            _, existing_energy = d_hybrids[key]
            if energy_float < existing_energy:
                d_hybrids[key] = (val, energy_float)
    d_hybrids = {k: v[0] for k, v in d_hybrids.items()}

    if not d_hybrids and data_lines:
        sys.stderr.write(f"Warning: No valid IntaRNA data rows in {csv_path}. First data line: {data_lines[0][:80]!r}\n")
    return d_hybrids



def hybridize_with_intarna(query_fa_path, target_fa_path, intarna_params, out_csv_path):
    """
    Run IntaRNA once on multi-FASTA query and target files (--outPairwise: one result per
    same-index pair). Used for local hybridization when --use_batchtools is not set.
    Writes CSV to out_csv_path. Returns dict (id1, id2) -> (dotbracket, pos, energy_str, subseq1, subseq2, seq_for_id1, seq_for_id2)
    via parse_intarna_csv (same format as batchtools).
    """
    cmd = ["IntaRNA"] + (intarna_params.split() if isinstance(intarna_params, str) else list(intarna_params))
    cmd += ["-q", query_fa_path, "-t", target_fa_path, "--out", out_csv_path]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if result.returncode != 0:
        return {}
    return parse_intarna_csv(out_csv_path)



def _load_loci_seqs_from_fasta(outdir, n, buffer_size=BUFFER_SIZE):
    """Load d_loci_seqs from loci.fa.<n> (shared helper)."""
    d_loci_seqs = {}
    fasta_file = os.path.join(outdir, "loci.fa.") + str(n)
    with open(fasta_file, 'r', buffering=buffer_size) as fh:
        current_id = None
        current_seq = []
        for line in fh:
            line_stripped = line.rstrip('\n\r')
            if line_stripped.startswith('>'):
                if current_id is not None:
                    locus_id = current_id[:-3] if len(current_id) >= 3 else current_id
                    sequence = ''.join(current_seq).upper().replace('T', 'U')
                    d_loci_seqs[locus_id] = sequence
                current_id = line_stripped[1:].strip()
                current_seq = []
            else:
                if current_id is not None:
                    current_seq.append(line_stripped)
        if current_id is not None:
            locus_id = current_id[:-3] if len(current_id) >= 3 else current_id
            sequence = ''.join(current_seq).upper().replace('T', 'U')
            d_loci_seqs[locus_id] = sequence
    return d_loci_seqs



def _collect_unique_locus_pairs(file_chimeras, d_loci_seqs, buffer_size, collect_lines=False):
    """
    One pass over a chimera chunk file: collect unique (locus1, locus2) pairs in first-seen order
    (only pairs where both loci are in d_loci_seqs). Optionally also collect every line with its
    locus ids for the batchtools finish phase.
    Returns (unique_pairs, lines_with_loci). If collect_lines is False, lines_with_loci is None.
    """
    unique_pairs = []
    seen = set()
    lines_with_loci = [] if collect_lines else None
    with open(file_chimeras, "r", buffering=buffer_size) as fh:
        for line in fh:
            a = line.rstrip("\n").split("\t")
            locuspos1 = a[CHIMERA_IDX_LOCUS1]
            locuspos2 = a[CHIMERA_IDX_LOCUS2]
            if collect_lines:
                lines_with_loci.append((line.rstrip("\n"), locuspos1, locuspos2))
            if locuspos1 not in d_loci_seqs or locuspos2 not in d_loci_seqs:
                continue
            key = (locuspos1, locuspos2)
            if key not in seen:
                seen.add(key)
                unique_pairs.append((locuspos1, locuspos2))
    return (unique_pairs, lines_with_loci)



def _write_paired_fasta(unique_pairs, d_loci_seqs, query_fa_path, target_fa_path, buffer_size):
    """
    Write multi-FASTA query.fa and target.fa for IntaRNA --outPairwise (one entry per pair;
    same index = same pair). We call IntaRNA with -q query.fa -t target.fa. IntaRNA uses
    index 1 = target, index 2 = query (README), so CSV id1=target_locus, id2=query_locus.
    We assign shorter locus to query, longer to target (IntaRNA recommendation).
    """
    with open(query_fa_path, "w", buffering=buffer_size) as fq, \
         open(target_fa_path, "w", buffering=buffer_size) as ft:
        for (lp1, lp2) in unique_pairs:
            s1, s2 = d_loci_seqs[lp1], d_loci_seqs[lp2]
            # Shorter as query (IntaRNA recommendation; typically miRNA side)
            if len(s1) <= len(s2):
                query_locus, target_locus = lp1, lp2
            else:
                query_locus, target_locus = lp2, lp1
            fq.write(">" + query_locus + "\n" + d_loci_seqs[query_locus] + "\n")
            ft.write(">" + target_locus + "\n" + d_loci_seqs[target_locus] + "\n")



def _merge_hybrid_into_chimeras(file_chimeras, file_out, d_loci_seqs, d_hybrids, buffer_size, remove_input=False):
    """
    Read chimera lines from file_chimeras, fill hybrid columns from d_hybrids (and optionally d_loci_seqs).
    Chimera columns use locus1 & locus2 order. d_hybrids keys are (id1, id2) with IntaRNA convention
    (id1=target, id2=query); we map to locus order and swap when key is (lp2, lp1).
    When remove_input is True (e.g. --remove_intermediate), deletes file_chimeras.
    """
    with open(file_chimeras, "r", buffering=buffer_size) as fh_in, \
         open(file_out, "w", buffering=buffer_size) as fh_out:
        for line in fh_in:
            a = line.rstrip("\n").split("\t")
            seq1 = seq2 = dotbracket = startpos = energy = hybrid_subseqs = "NA"
            locuspos1 = a[CHIMERA_IDX_LOCUS1]
            locuspos2 = a[CHIMERA_IDX_LOCUS2]
            lp1 = locuspos1.strip() if locuspos1 else locuspos1
            lp2 = locuspos2.strip() if locuspos2 else locuspos2
            # Resolve full sequences from d_loci_seqs (full locus FASTA); d_hybrids only provides structure/energy
            if (lp1, lp2) in d_hybrids:
                tup = d_hybrids[(lp1, lp2)]
                dotbracket, startpos, energy, subseq1, subseq2 = tup[0], tup[1], tup[2], tup[3], tup[4]
                hybrid_subseqs = subseq1 + "&" + subseq2
                if d_loci_seqs and lp1 in d_loci_seqs and lp2 in d_loci_seqs:
                    seq1, seq2 = d_loci_seqs[lp1], d_loci_seqs[lp2]
            elif (lp2, lp1) in d_hybrids:
                # Stored key (id1, id2) = (lp2, lp1): dotbracket/pos/subseqs are lp2&lp1; we need lp1&lp2
                tup = d_hybrids[(lp2, lp1)]
                db, p, energy, subseq_id1, subseq_id2 = tup[0], tup[1], tup[2], tup[3], tup[4]
                part_id1, part_id2 = db.split("&", 1)  # lp2 part, lp1 part (each has one bracket type)
                dotbracket = part_id2.replace(")", "(") + "&" + part_id1.replace("(", ")")  # lp1 & lp2
                startpos = p.split("&", 1)[1] + "&" + p.split("&", 1)[0]  # lp1 & lp2
                hybrid_subseqs = subseq_id2 + "&" + subseq_id1  # id2=lp1, id1=lp2 → lp1 & lp2
                if d_loci_seqs and lp1 in d_loci_seqs and lp2 in d_loci_seqs:
                    seq1, seq2 = d_loci_seqs[lp1], d_loci_seqs[lp2]
            elif d_loci_seqs and lp1 in d_loci_seqs and lp2 in d_loci_seqs:
                seq1 = d_loci_seqs[lp1]
                seq2 = d_loci_seqs[lp2]
            while len(a) <= CHIMERA_IDX_HYBRID_SUBSEQS:
                a.append("NA")
            a[CHIMERA_IDX_SEQUENCES] = seq1 + "&" + seq2
            a[CHIMERA_IDX_HYBRID] = dotbracket
            a[CHIMERA_IDX_HYBRID_POS] = startpos
            a[CHIMERA_IDX_MFE] = energy
            a[CHIMERA_IDX_HYBRID_SUBSEQS] = hybrid_subseqs
            fh_out.write("\t".join(a) + "\n")
    if remove_input and os.path.exists(file_chimeras):
        os.remove(file_chimeras)



def hybridize_and_write(outdir, intarna_params, n, sample_name, remove_intermediate=False):
    """
    Local hybridization: build multi-FASTA query/target from unique locus pairs, run IntaRNA
    once with --outPairwise, write chimeras-r.<n>. When remove_intermediate: deletes chimeras.<n> and local_intarna_<n>.
    Used when --hybridize and not --use_batchtools.
    """
    d_loci_seqs = _load_loci_seqs_from_fasta(outdir, n)
    file_chimeras = os.path.join(outdir, sample_name + ".chimeras." + str(n))
    file_out = os.path.join(outdir, sample_name + ".chimeras-r." + str(n))

    unique_pairs, _ = _collect_unique_locus_pairs(file_chimeras, d_loci_seqs, BUFFER_SIZE, collect_lines=False)

    d_hybrids = {}
    if unique_pairs:
        work_dir = os.path.join(outdir, "local_intarna_" + str(n))
        os.makedirs(work_dir, exist_ok=True)
        try:
            query_fa = os.path.join(work_dir, "query.fa")
            target_fa = os.path.join(work_dir, "target.fa")
            _write_paired_fasta(unique_pairs, d_loci_seqs, query_fa, target_fa, BUFFER_SIZE)
            out_csv = os.path.join(work_dir, "intarna_out" + str(n) + ".csv")
            d_hybrids = hybridize_with_intarna(query_fa, target_fa, intarna_params, out_csv)
        finally:
            if remove_intermediate:
                try:
                    shutil.rmtree(work_dir, ignore_errors=True)
                except OSError:
                    pass

    _merge_hybrid_into_chimeras(file_chimeras, file_out, d_loci_seqs, d_hybrids, BUFFER_SIZE, remove_input=remove_intermediate)



def prepare_hybridization_batch(outdir, intarna_params, n, sample_name, batchtools_work_dir):
    """
    Phase 1 for batchtools: load chimera chunk (*.chimeras.n), write multi-FASTA query.fa and
    target.fa for one IntaRNA --outPairwise run per chunk. Does not run IntaRNA (worker does).
    Saves d_loci_seqs to chunk dir so finish phase can retrieve full locus sequences even when
    IntaRNA seq1/seq2 are NA for a pair.
    """
    d_loci_seqs = _load_loci_seqs_from_fasta(outdir, n, BUFFER_SIZE)
    file_chimeras = os.path.join(outdir, sample_name + ".chimeras." + str(n))
    chunk_dir = os.path.join(batchtools_work_dir, str(n))
    os.makedirs(chunk_dir, exist_ok=True)

    unique_pairs, _ = _collect_unique_locus_pairs(file_chimeras, d_loci_seqs, BUFFER_SIZE, collect_lines=False)

    query_fa = os.path.join(chunk_dir, "query.fa")
    target_fa = os.path.join(chunk_dir, "target.fa")
    _write_paired_fasta(unique_pairs, d_loci_seqs, query_fa, target_fa, BUFFER_SIZE)
    # Persist full locus sequences for finish phase (avoids re-reading loci.fa.<n> and covers NA seq1/seq2)
    with open(os.path.join(chunk_dir, "loci_seqs.pkl"), "wb") as f:
        pickle.dump(d_loci_seqs, f)



def finish_hybridization_write(outdir, n, sample_name, batchtools_work_dir, compress=False, remove_intermediate=False):
    """
    Phase 3 for batchtools: read *.chimeras.n and IntaRNA result CSV, write chimeras-r.<n>.
    Uses IntaRNA CSV (seq1, seq2) when available and falls back to d_loci_seqs loaded from
    loci_seqs.pkl in the chunk dir. When remove_intermediate, deletes *.chimeras.n.
    """
    chunk_dir = os.path.join(batchtools_work_dir, str(n))
    file_chimeras = os.path.join(outdir, sample_name + ".chimeras." + str(n))
    output_file = os.path.join(outdir, sample_name + ".chimeras-r." + str(n))
    loci_pkl = os.path.join(chunk_dir, "loci_seqs.pkl")
    d_loci_seqs = {}
    if os.path.exists(loci_pkl):
        with open(loci_pkl, "rb") as f:
            d_loci_seqs = pickle.load(f)
    result_csv = os.path.join(chunk_dir, "result.csv")
    d_hybrids = parse_intarna_csv(result_csv) if os.path.exists(result_csv) else {}
    _merge_hybrid_into_chimeras(file_chimeras, output_file, d_loci_seqs or None, d_hybrids, BUFFER_SIZE, remove_input=remove_intermediate)



def run_batchtools_r_script(args, batchtools_work_dir, batchtools_registry, intarna_params):
    """
    Submit IntaRNA jobs via R batchtools (same strategy as chira_map.py). One LSF job per chunk;
    each job runs process_intarna_chunk_batchtools.py, which runs IntaRNA once per chunk with
    multi-FASTA query.fa and target.fa (--outPairwise: one result per same-index pair).
    Writes config.json and chunks.json to registry dir, then calls submit_intarna_batchtools.R.
    """
    # Build chunks list: one entry per chunk dir that has query.fa and target.fa (had pairs)
    chunks_data = []
    for n in sorted([d for d in os.listdir(batchtools_work_dir) if os.path.isdir(os.path.join(batchtools_work_dir, d)) and d.isdigit()], key=int):
        chunk_dir = os.path.join(batchtools_work_dir, n)
        if os.path.exists(os.path.join(chunk_dir, "query.fa")) and os.path.exists(os.path.join(chunk_dir, "target.fa")):
            chunks_data.append({"chunk_idx": int(n)})
    if not chunks_data:
        return

    reg_dir = os.path.abspath(batchtools_registry)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(reg_dir, exist_ok=True)

    # Resource defaults (match chira_map.py)
    cores_per_job = getattr(args, 'batchtools_cores', None) or 8
    if not getattr(args, 'batchtools_cores', None):
        print("WARNING: --batchtools_cores not specified. Using default: 8 cores per job.", file=sys.stderr)
    if getattr(args, 'batchtools_memory', None):
        mem_match = re.match(r'(\d+)(GB|G|MB|M)', (getattr(args, 'batchtools_memory') or '').upper())
        if mem_match:
            mem_value = int(mem_match.group(1))
            mem_unit = mem_match.group(2)
            total_gb = mem_value / 1024 if mem_unit in ('MB', 'M') else mem_value
            memory_per_job = f"{total_gb / cores_per_job:.1f}GB"
        else:
            memory_per_job = f"{max(0.5, 4.0 / cores_per_job):.1f}GB"
    else:
        memory_per_job = f"{max(0.5, 4.0 / cores_per_job):.1f}GB"
    walltime = getattr(args, 'batchtools_walltime', None) or "48:00"
    queue = getattr(args, 'batchtools_queue', None) or "long"
    conda_env = getattr(args, 'batchtools_conda_env', None) or os.environ.get('CONDA_DEFAULT_ENV', '')
    if conda_env:
        conda_env = os.path.expanduser(conda_env)
    max_parallel = getattr(args, 'batchtools_max_parallel', None)
    job_name_prefix = f"chira_intarna_{os.path.basename(args.outdir)}_{int(time.time())}"

    # Resolve template (same as chira_map.py)
    template_file = getattr(args, 'batchtools_template', None) or ''
    if template_file == 'lsf-simple':
        pass
    elif template_file:
        if not os.path.isabs(template_file):
            template_file = os.path.join(script_dir, template_file)
        template_file = os.path.abspath(template_file)
        if not os.path.exists(template_file):
            print(f"WARNING: Template not found: {template_file}. Using 'lsf-simple'.", file=sys.stderr)
            template_file = "lsf-simple"
    else:
        default_tmpl = os.path.join(script_dir, "lsf_custom.tmpl")
        template_file = default_tmpl if os.path.exists(default_tmpl) else "lsf-simple"
    if template_file != "lsf-simple":
        template_file = os.path.abspath(template_file)

    python_script = os.path.abspath(os.path.join(script_dir, "process_intarna_chunk_batchtools.py"))
    if not os.path.exists(python_script):
        raise FileNotFoundError(f"Python worker script not found: {python_script}")

    config = {
        "reg_dir": reg_dir,
        "queue": queue,
        "cores_per_job": cores_per_job,
        "memory_per_job": memory_per_job,
        "walltime": walltime,
        "conda_env": conda_env,
        "intarna_params": intarna_params,
        "job_name_prefix": job_name_prefix,
        "max_parallel": max_parallel,
        "template_file": template_file,
        "batchtools_work_dir": os.path.abspath(batchtools_work_dir),
        "python_script": python_script,
    }

    config_file = os.path.join(reg_dir, "config.json")
    chunks_file = os.path.join(reg_dir, "chunks.json")
    with open(config_file, 'w', encoding='utf-8') as f:
        json.dump(config, f, indent=2, ensure_ascii=False)
    with open(chunks_file, 'w', encoding='utf-8') as f:
        json.dump(chunks_data, f, indent=2, ensure_ascii=False)

    r_script = os.path.join(script_dir, "submit_intarna_batchtools.R")
    if not os.path.exists(r_script):
        raise FileNotFoundError(f"R script not found: {r_script}")

    cmd = ["Rscript", r_script, os.path.abspath(config_file), os.path.abspath(chunks_file)]
    print("Submitting IntaRNA jobs via batchtools (chira_map pattern): " + " ".join(cmd), file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=86400, cwd=script_dir)
    if result.returncode != 0:
        sys.stderr.write(result.stderr or "")
        raise RuntimeError(f"batchtools R script failed with code {result.returncode}. stderr: {result.stderr!r}")



def run_batchtools_hybridization(args, common_intarna_params, batchtools_work_dir, batchtools_registry):
    """Run hybridization using batchtools: prepare (MPIRE) -> R batchtools -> finish (MPIRE)."""
    # Phase 1: prepare FASTA and manifests per chunk
    prep_args = [(
        args.outdir,
        common_intarna_params,
        str(k),
        args.sample_name,
        batchtools_work_dir
    ) for k in range(args.processes)]
    with WorkerPool(n_jobs=args.processes) as pool:
        pool.map(prepare_hybridization_batch, prep_args, progress_bar=False)
    # Phase 2: run R batchtools (submits to cluster, waits for completion; same pattern as chira_map)
    run_batchtools_r_script(args, batchtools_work_dir, batchtools_registry, common_intarna_params)
    # Phase 3: write final output from result CSVs
    finish_args = [(
        args.outdir,
        str(k),
        args.sample_name,
        batchtools_work_dir,
        args.compress,
        getattr(args, 'remove_intermediate', False)
    ) for k in range(args.processes)]
    with WorkerPool(n_jobs=args.processes) as pool:
        pool.map(finish_hybridization_write, finish_args, progress_bar=False)



# --- Chunk writing ---

def write_chimeras(chunk_start, chunk_end, total_read_count, d_ref_lengths1, d_ref_lengths2, hybridize,
                   chimeric_overlap, f_gtf, outdir, crl_file, tpm_threshold, score_cutoff, n, sample_name, compress=False, shared_objects=None):
    """
    Write chimeric reads and singletons for a chunk of reads.
    
    OPTIMIZATION: If shared_objects is provided, use shared memory for d_ref_lengths1 and d_ref_lengths2.
    This reduces memory overhead by 50-90% when using MPIRE.
    
    Args:
        shared_objects: Optional dict containing shared objects (d_ref_lengths1, d_ref_lengths2)
                       If provided, d_ref_lengths1 and d_ref_lengths2 arguments are ignored
    """
    # OPTIMIZATION: Use shared objects if available (MPIRE), otherwise use passed arguments
    if shared_objects:
        d_ref_lengths1 = shared_objects['d_ref_lengths1']
        d_ref_lengths2 = shared_objects['d_ref_lengths2']
    d_regions = {}
    l_loci_bed = set()
    file_chimeras = os.path.join(outdir, sample_name + ".chimeras." + n)
    file_singletons = os.path.join(outdir, sample_name + ".singletons." + n)
    # Intermediate files are NOT compressed - only final merged files are compressed
    # This avoids CPU overhead during parallel processing and merge operations
    
    # Open files (intermediate files are always uncompressed)
    open_func = open
    open_mode = "w"
    
    # OPTIMIZATION: Use indexed file access to jump directly to chunk start
    # This avoids scanning the entire file sequentially (critical for large files)
    index_file = crl_file + '.idx'
    use_index = os.path.exists(index_file)
    
    if use_index:
        # Load index
        try:
            with open(index_file, 'rb') as f:
                read_index_dict, read_list = pickle.load(f)
            
            # Get read IDs for this chunk
            # Original logic analysis:
            # - Processes reads where: chunk_start <= read_count + 1 AND read_count <= chunk_end
            # - This translates to: read_count from max(0, chunk_start-1) to chunk_end (inclusive)
            # - Since read_count directly maps to read_list indices (0-indexed), we need:
            #   read_list[start_idx:end_idx+1] where start_idx = max(0, chunk_start-1) and end_idx = chunk_end
            start_idx = max(0, chunk_start - 1)
            end_idx = min(len(read_list) - 1, chunk_end)
            
            if start_idx <= end_idx and end_idx < len(read_list):
                chunk_read_ids = read_list[start_idx:end_idx+1]  # Include end
                start_pos = read_index_dict[chunk_read_ids[0]][0] if chunk_read_ids else 0
                end_read_id = chunk_read_ids[-1] if chunk_read_ids else None
            else:
                # Fallback to sequential if index doesn't cover range
                use_index = False
                chunk_read_ids = None
        except (IOError, KeyError, IndexError) as e:
            print(f"Warning: Could not use index file ({e}), falling back to sequential access", file=sys.stderr)
            use_index = False
            chunk_read_ids = None
    
    # make bed entry for extracing locus sequence and hybridizing
    with open_func(file_chimeras, open_mode, buffering=BUFFER_SIZE) as fh_chimeras, \
         open_func(file_singletons, open_mode, buffering=BUFFER_SIZE) as fh_singletons:
        
        if use_index and chunk_read_ids:
            # OPTIMIZATION: Jump directly to chunk start using index
            # This avoids scanning the entire file sequentially (8-10x speedup for large files)
            chunk_read_ids_set = set(chunk_read_ids)
            with open(crl_file, "rb", buffering=BUFFER_SIZE) as fh_crl:
                fh_crl.seek(start_pos)
                # Ensure we're at the start of a line (in case seek landed in middle of line)
                # Read one byte to check, but this should already be at line start from index
                
                prev_readid = None
                l_readlines = []
                last_read_processed = False
                
                for line_bytes in fh_crl:
                    line_str = line_bytes.decode('utf-8', errors='replace')
                    f = line_str.rstrip('\n').split('\t')
                    if len(f) < 1:
                        continue
                    
                    # last field after | represents the segment id, rest of the string before is read id
                    readid = '|'.join(f[0].split("|")[:-1])
                    
                    if readid != prev_readid:
                        # Process previous read if it's in our chunk
                        if prev_readid is not None and prev_readid in chunk_read_ids_set:
                            l_read_alignments = filter_alignments(l_readlines, tpm_threshold, score_cutoff)
                            extract_and_write(prev_readid, l_read_alignments, l_loci_bed, d_ref_lengths1, d_ref_lengths2, f_gtf,
                                              d_regions, chimeric_overlap, fh_chimeras, fh_singletons, hybridize)
                            
                            # Check if this was the last read in our chunk
                            if prev_readid == end_read_id:
                                last_read_processed = True
                                break
                        
                        prev_readid = readid
                        l_readlines = []
                        
                        # If we've moved past our chunk (new readid is not in our chunk), stop reading
                        # But only if we've already processed all reads in our chunk
                        if readid not in chunk_read_ids_set:
                            # We've moved past our chunk, stop reading
                            break
                    
                    # Accumulate lines for reads in our chunk
                    if readid in chunk_read_ids_set:
                        l_readlines.append(line_str)
                
                # Process last read if we have lines and haven't processed it yet
                if not last_read_processed and prev_readid is not None and prev_readid in chunk_read_ids_set and l_readlines:
                    l_read_alignments = filter_alignments(l_readlines, tpm_threshold, score_cutoff)
                    extract_and_write(prev_readid, l_read_alignments, l_loci_bed, d_ref_lengths1, d_ref_lengths2, f_gtf,
                                      d_regions, chimeric_overlap, fh_chimeras, fh_singletons, hybridize)
        else:
            # Fallback: Sequential access (original method)
            with open(crl_file, "r", buffering=BUFFER_SIZE) as fh_crl:
                prev_readid = None
                l_readlines = []
                read_count = 0
                for line in fh_crl:
                    f = line.rstrip('\n').split('\t')
                    # last field after | represents the segment id, rest of the string before is read id
                    readid = '|'.join(f[0].split("|")[:-1])

                    if prev_readid != readid:
                        if chunk_start > read_count + 1:
                            read_count += 1
                            prev_readid = readid
                            l_readlines = []
                            continue
                        if read_count > chunk_end:
                            break
                        if prev_readid is not None:
                            l_read_alignments = filter_alignments(l_readlines, tpm_threshold, score_cutoff)
                            extract_and_write(prev_readid, l_read_alignments, l_loci_bed, d_ref_lengths1, d_ref_lengths2, f_gtf,
                                              d_regions, chimeric_overlap, fh_chimeras, fh_singletons, hybridize)

                        read_count += 1
                        prev_readid = readid
                        l_readlines = []
                    l_readlines.append(line)

                # write last read info if it's within our chunk
                # At end of file, the last read is still in l_readlines and hasn't been processed
                # Process it if it's within our chunk range
                if l_readlines and chunk_start <= read_count + 1 <= chunk_end:
                    l_read_alignments = filter_alignments(l_readlines, tpm_threshold, score_cutoff)
                    extract_and_write(readid, l_read_alignments, l_loci_bed, d_ref_lengths1, d_ref_lengths2, f_gtf,
                                      d_regions, chimeric_overlap, fh_chimeras, fh_singletons, hybridize)

    if hybridize:
        # loci sequences are neeeded to hybridize
        with open(os.path.join(outdir, "loci.bed.") + str(n), "w", buffering=BUFFER_SIZE) as fh_bed:
            for bed_line in l_loci_bed:
                fh_bed.write(bed_line + "\n")




# --- Counts and merge ---

def parse_counts_file(crl_file, tpm_cutoff, rebuild_index=False):
    d_crl_tpm = defaultdict(float)
    l_loci_bed = set()
    prev_readid = None
    read_count = 0
    index_file = crl_file + '.idx'
    use_existing_index = (not rebuild_index) and os.path.exists(index_file)
    if use_existing_index:
        print(f"Using existing index: {index_file}", file=sys.stderr)

    # Create read index while parsing (unless using existing index)
    read_index_dict = {}
    read_list = []
    current_pos = 0
    line_count_for_read = 0

    with open(crl_file, "rb", buffering=BUFFER_SIZE) as fh_crl_file:
        for line_bytes in fh_crl_file:
            line_str = line_bytes.decode('utf-8', errors='replace')
            f = line_str.rstrip('\n').split('\t')
            if len(f) < 13:
                current_pos = fh_crl_file.tell()
                continue

            readid = '|'.join(f[0].split("|")[:-1])

            if readid != prev_readid:
                if not use_existing_index:
                    if prev_readid is not None:
                        read_index_dict[prev_readid] = (read_index_dict[prev_readid][0], line_count_for_read)
                        line_count_for_read = 0
                    if readid not in read_index_dict:
                        read_index_dict[readid] = (current_pos, 0)
                        read_list.append(readid)

                read_count += 1
                prev_readid = readid
                if use_existing_index:
                    line_count_for_read = 0

            if not use_existing_index:
                line_count_for_read += 1

            crlid = f[3]
            crl_tpm = f[12]
            d_crl_tpm[crlid] = float(crl_tpm)
            b = f[9].split(":")
            locus_bed_entry = "\t".join([":".join(b[0:-3]), b[-3], b[-2], f[9], "1", b[-1]])
            if locus_bed_entry not in l_loci_bed:
                l_loci_bed.add(locus_bed_entry)

            current_pos = fh_crl_file.tell()

        if not use_existing_index and prev_readid is not None:
            read_index_dict[prev_readid] = (read_index_dict[prev_readid][0], line_count_for_read)

    if use_existing_index:
        # Verify existing index matches current CRL (same number of reads); remove if stale
        try:
            with open(index_file, 'rb') as f:
                _, existing_read_list = pickle.load(f)
            if len(existing_read_list) != read_count:
                print(f"Warning: Index has {len(existing_read_list):,} reads but CRL has {read_count:,}; removing stale index", file=sys.stderr)
                try:
                    os.remove(index_file)
                except OSError as err:
                    print(f"Warning: Could not remove stale index: {err}", file=sys.stderr)
        except (IOError, Exception) as e:
            print(f"Warning: Could not verify index ({e}); leaving as-is", file=sys.stderr)
    else:
        # Save index for use by write_chimeras
        print(f"Saving read index to {index_file}...", file=sys.stderr)
        with open(index_file, 'wb') as f:
            pickle.dump((read_index_dict, read_list), f)
        print(f"Created index for {len(read_list):,} reads", file=sys.stderr)

    uniq_tpms = sorted(list(set(d_crl_tpm.values())))
    
    # Clamp tpm_cutoff to [0, 1) to prevent index out of bounds
    if tpm_cutoff < 0:
        tpm_cutoff = 0.0
    elif tpm_cutoff >= 1.0:
        tpm_cutoff = 0.999999  # Just below 1.0 to ensure valid index
    
    # Handle empty list case
    if len(uniq_tpms) == 0:
        tpm_threshold = 0.0
    else:
        # Calculate index: ensure it's in valid range [0, len-1]
        index = int(tpm_cutoff * len(uniq_tpms))
        # Clamp index to valid range (shouldn't be needed with tpm_cutoff < 1.0, but defensive)
        index = min(index, len(uniq_tpms) - 1)
        tpm_threshold = uniq_tpms[index]

    return read_count, tpm_threshold



def merge_files(inprefix, outfile, header, r, compress=False, remove_intermediate=False):
    # Intermediate files are always uncompressed; only final output is compressed when requested
    temp_files = []
    for i in range(r):
        temp_file = inprefix + "." + str(i)
        if os.path.exists(temp_file):
            temp_files.append(temp_file)
    if not temp_files:
        return
    
    # Write header first, then append sorted unique lines
    # Note: gzip.open doesn't support buffering parameter; use BUFFER_SIZE for uncompressed output only
    open_func = gzip.open if compress else open
    open_mode = "wt" if compress else "w"
    
    # Only apply buffering for uncompressed files (gzip.open handles buffering internally)
    if compress:
        fh_out = open_func(outfile, open_mode)
    else:
        fh_out = open_func(outfile, open_mode, buffering=BUFFER_SIZE)
    
    with fh_out:
        # Write header first
        fh_out.write(header + "\n")
        
        # Use shell commands to merge and sort files (more efficient for large files)
        # Intermediate files are uncompressed, so we can use cat directly
        # OPTIMIZATION: Use parallel sort if available (GNU sort >= 8.6); use more threads to shorten merge walltime
        sort_cmd = "sort -u"
        try:
            result = subprocess.run(["sort", "--version"], capture_output=True, text=True, timeout=2)
            if "GNU coreutils" in result.stdout:
                n_sort = max(1, min(r * 2, (os.cpu_count() or 16)))  # more threads when available to speed merge
                sort_cmd = f"sort --parallel={n_sort} -u"
        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
            pass
        
        escaped_files = " ".join([f"'{f}'" for f in temp_files])
        cmd = f"cat {escaped_files} | {sort_cmd}"
        # Execute command and write output to file
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            sys.stderr.write(f"Warning: merge_files command failed: {stderr}\n")
        else:
            # Write sorted unique lines (header already written)
            # Compression happens automatically if open_func is gzip.open
            fh_out.write(stdout)

    if remove_intermediate:
        for i in range(r):
            temp_file = inprefix + "." + str(i)
            if os.path.exists(temp_file):
                os.remove(temp_file)


def write_interaction_summary(outdir, sample_name, compress=False, num_threads=4):
    d_interactions = defaultdict(lambda: defaultdict(list))
    chimeras_file = os.path.join(outdir, sample_name + ".chimeras.txt")
    if compress:
        chimeras_file += ".gz"

    # Note: gzip.open doesn't support buffering parameter; use BUFFER_SIZE for uncompressed input only
    open_func = gzip.open if compress else open
    open_mode = "rt" if compress else "r"

    # Only apply buffering for uncompressed files (gzip.open handles buffering internally)
    if compress:
        fh_in = open_func(chimeras_file, open_mode)
    else:
        fh_in = open_func(chimeras_file, open_mode, buffering=BUFFER_SIZE)

    with fh_in:
        next(fh_in)
        for line in fh_in:
            f = line.rstrip("\n").split("\t")
            (readid, ref1, ref2, region1, region2, 
            locus1, locus2, tpm1, tpm2, score1, score2) = (f[0], f[1], f[2],
                                                            f[7], f[8], f[20],
                                                            f[21], f[24], f[25],
                                                            f[26], f[27])
            gene_name1 = f[5] 
            gene_name2 = f[6] 
            tpm = str(float(tpm1) + float(tpm2))
            score = str(float(score1) * float(score2))
            sequence1 = sequence2 = "NA"
            if f[29] != "NA":
                [sequence1, sequence2] = f[29].split("&")
            dotbracket = f[30]
            hybrid_start_pos = f[31]
            mfe = f[32]
            interaction = "\t".join(locus1.split(":")) + "\t" + "\t".join(locus2.split(":"))
            hybridization_pos = interaction
            [refid1, ref_start1, re1, ref_strand1, refid2, ref_start2, re2, ref_strand2] = interaction.split("\t")
            hybridized_sequences = "NA\tNA"
            interaction_otherway = "\t".join(locus2.split(":")) + "\t" + "\t".join(locus1.split(":"))

            if interaction_otherway in d_interactions:
                interaction = interaction_otherway
                hybridization_pos = interaction
                (ref2, ref1, region2, region1, tpm2, tpm1, score2, score1) = (f[1], f[2], f[7], f[8], f[24], f[25],
                                                                                f[26], f[27])
                gene_name1, gene_name2 = gene_name2, gene_name1
            if dotbracket != "NA":
                hybrid_subseqs = f[34] if len(f) > 34 else "NA"
                hybridization_pos1 = hybridization_pos2 = ""
                hybridized_sequence1 = hybridized_sequence2 = "NA"
                if hybrid_subseqs != "NA" and "&" in hybrid_subseqs:
                    # Column 34 and 31 are locus1&locus2: subseq1&subseq2, start1&start2
                    hybridized_sequence1, hybridized_sequence2 = hybrid_subseqs.split("&", 1)
                    start1, start2 = hybrid_start_pos.split("&", 1)
                    start1, start2 = int(start1), int(start2)  # 1-based
                    end1 = start1 + len(hybridized_sequence1) - 1
                    end2 = start2 + len(hybridized_sequence2) - 1
                    hybridization_pos1 = f'{refid1}:{int(ref_start1) + start1 - 1}-{int(ref_start1) + end1}:{ref_strand1}'
                    hybridization_pos2 = f'{refid2}:{int(ref_start2) + start2 - 1}-{int(ref_start2) + end2}:{ref_strand2}'
                    hybridization_pos = hybridization_pos1 + '&' + hybridization_pos2
                    hybridized_sequences = hybridized_sequence1 + '&' + hybridized_sequence2

                if interaction_otherway in d_interactions:
                    hybrid_start_pos = f[31].split('&')[1] + '&' + f[31].split('&')[0]
                    hybridization_pos = hybridization_pos2 + '&' + hybridization_pos1
                    hybridized_sequences = hybridized_sequence2 + '&' + hybridized_sequence1
                    target_dotbracket = dotbracket.split("&")[0].replace("(", ")")
                    query_dotbracket = dotbracket.split("&")[1].replace(")", "(")
                    dotbracket = query_dotbracket + "&" + target_dotbracket
                    [sequence2, sequence1] = f[29].split("&")

            d_interactions[interaction]["readid"].append(readid)
            d_interactions[interaction]["ref1"].append(ref1)
            d_interactions[interaction]["ref2"].append(ref2)
            d_interactions[interaction]["region1"].append(region1)
            d_interactions[interaction]["region2"].append(region2)
            d_interactions[interaction]["gene_name1"].append(gene_name1)
            d_interactions[interaction]["gene_name2"].append(gene_name2)
            common_info = "\t".join([sequence1, sequence2, dotbracket, mfe, hybridized_sequences,
                                        hybrid_start_pos, hybridization_pos,
                                        tpm1, tpm2, tpm, score1, score2, score])
            d_interactions[interaction]["common"] = [common_info]

    with open(os.path.join(outdir, "interactions.temp"), "w", buffering=BUFFER_SIZE) as fh_out:
        for interaction in d_interactions.keys():
            fh_out.write("\t".join([str(len(set(d_interactions[interaction]["readid"]))),
                                    interaction,
                                    d_interactions[interaction]["common"][0],
                                    ";".join(sorted(set(d_interactions[interaction]["region1"]))),
                                    ";".join(sorted(set(d_interactions[interaction]["region2"]))),
                                    ";".join(sorted(set(d_interactions[interaction]["ref1"]))),
                                    ";".join(sorted(set(d_interactions[interaction]["ref2"]))),
                                    ";".join(sorted(set(d_interactions[interaction]["gene_name1"]))),
                                    ";".join(sorted(set(d_interactions[interaction]["gene_name2"])))]) + "\n")

    header_interactions = "\t".join([
        "supporting_read_count",
        "locus_1_chromosome", "locus_1_start", "locus_1_end", "locus_1_strand",
        "locus_2_chromosome", "locus_2_start", "locus_2_end", "locus_2_strand",
        "locus_1_sequence", "locus_2_sequence", "hybridization_structure_dotbracket",
        "hybridization_mfe_kcal_mol", "hybridized_sequence_segments",
        "hybridization_start_positions","hybridization_genomic_coordinates",
        "tpm_locus_1", "tpm_locus_2", "tpm_combined",
        "alignment_score_locus_1", "alignment_score_locus_2", "combined_alignment_score",
        "annotation_region_locus_1", "annotation_region_locus_2",
        "reference_transcript_id_1", "reference_transcript_id_2",
        "gene_name_1", "gene_name_2"
    ])
    # OPTIMIZATION: Use parallel sort if available (GNU sort supports --parallel option)
    # For systems with GNU sort >= 8.6, this can significantly speed up sorting large interaction files
    sort_cmd = "sort -k 1nr,1"
    try:
        # Check if GNU sort with --parallel is available (GNU sort >= 8.6)
        result = subprocess.run(["sort", "--version"], capture_output=True, text=True, timeout=2)
        if "GNU coreutils" in result.stdout:
            # Use parallel sort with specified number of threads
            # Default to 4 threads, but can be overridden for better performance
            sort_cmd = f"sort --parallel={max(1, num_threads)} -k 1nr,1"
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
        # Fall back to standard sort if check fails
        pass

    os.system(sort_cmd + " " + os.path.join(outdir, "interactions.temp") + " > " + os.path.join(outdir, "interactions.sorted"))
    with open(os.path.join(outdir, sample_name + ".interactions.txt"), "w", buffering=BUFFER_SIZE) as fh_out:
        fh_out.write("# Note: To identify which locus is miRNA vs target, check the region1 and region2 fields.\n")
        fh_out.write("# miRNA annotations include: miRNA, 3p_mature_mir, 5p_mature_mir, mature_mir\n")
        fh_out.write(header_interactions + "\n")
        with open(os.path.join(outdir, "interactions.sorted"), "r", buffering=BUFFER_SIZE) as fh_in:
            for line in fh_in:
                fh_out.write(line)
    os.remove(os.path.join(outdir, "interactions.temp"))
    os.remove(os.path.join(outdir, "interactions.sorted"))




# --- CLI and main flow ---

def validate_arguments(args):
    """Validate command-line arguments."""
    if args.hybridize and args.f_gtf and not args.f_ref:
        sys.stderr.write("Need the reference fasta file to hybridize. Make sure to provide the genomic fasta file"
                         " in case you already provided a GTF file.\n")
        sys.exit(1)

    if args.temperature < 0 or args.temperature > 100:
        sys.stderr.write("IntaRNA tempertature must be between 0 and 100!\n")
        sys.exit(1)



def print_configuration(args):
    """Print configuration parameters."""
    print('CRL file                             : ' + args.crl_file)
    print('Output directory                     : ' + args.outdir)
    if args.f_gtf:
        print('Annotation file                      : ' + args.f_gtf)
    print('Number of processes                  : ' + str(args.processes))
    print('TPM cutoff                           : ' + str(args.tpm_cutoff))
    print('Score cutoff                         : ' + str(args.score_cutoff))
    print('Chimeric overlap                     : ' + str(args.chimeric_overlap))
    print('Hybridize chimeric loci?             : ' + str(args.hybridize))
    print('Do not enforce seed interaction      : ' + str(args.no_seed))
    print('1st priority reference fasta file    : ' + args.ref_fasta1)
    if args.ref_fasta2:
        print('2nd priority reference fasta file    : ' + args.ref_fasta2)
    if args.f_ref:
        print('Reference genomic fasta file         : ' + args.f_ref)
    print('Summarize interactions at loci level : ' + str(args.summarize))
    print('Sample name                            : ' + args.sample_name)
    print('Compress output files with gzip        : ' + str(args.compress))
    print("===================================================================")



def setup_references(args):
    """Extract reference lengths from fasta files."""
    d_reflen1 = {}
    d_reflen2 = {}
    chira_utilities.extract_reflengths(args.ref_fasta1, d_reflen1)
    if args.ref_fasta2:
        chira_utilities.extract_reflengths(args.ref_fasta2, d_reflen2)
    return d_reflen1, d_reflen2



def run_chimera_extraction(args, d_reflen1, d_reflen2, tpm_cutoff_value, no_of_reads):
    """
    Run multiprocessing for chimera extraction.
    
    OPTIMIZATION: Uses MPIRE with shared objects if available for better performance.
    Falls back to multiprocessing.Process if MPIRE is not installed.
    """
    print(str(datetime.datetime.now()), " START: multiprocessing")
    
    # Prepare chunk arguments as tuples (not dictionaries) to avoid MPIRE unpacking issues
    # MPIRE will unpack dictionaries as keyword arguments, causing TypeError
    # Using tuples ensures arguments are passed as positional arguments
    chunk_args = []
    for k in range(args.processes):
        s = k * math.ceil(no_of_reads / args.processes)
        e = min(s + math.floor(no_of_reads / args.processes), no_of_reads)
        print(k, s, e, no_of_reads)
        chunk_args.append((
            s,  # chunk_start
            e,  # chunk_end
            no_of_reads,  # total_read_count
            args.hybridize,  # hybridize
            args.chimeric_overlap,  # chimeric_overlap
            args.f_gtf,  # f_gtf
            args.outdir,  # outdir
            args.crl_file,  # crl_file
            tpm_cutoff_value,  # tpm_threshold
            args.score_cutoff,  # score_cutoff
            str(k),  # n
            args.sample_name,  # sample_name
            args.compress  # compress
        ))
    
    # OPTIMIZATION: Use MPIRE with shared objects for better performance
    # Shared objects reduce memory overhead by 50-90% for large datasets
    # d_ref_lengths1 and d_ref_lengths2 are shared in memory instead of pickled per process
    # Pass shared objects directly to WorkerPool (MPIRE 2.10.1+)
    # Shared objects are passed as first argument to worker functions
    # Use 'fork' start method for copy-on-write (most memory efficient) when available
    # On Windows, fall back to default start method (spawn/forkserver)
    shared_objects_dict = {
        'd_ref_lengths1': d_reflen1,
        'd_ref_lengths2': d_reflen2
    }
    
    def process_chunk(shared_objects, chunk_start, chunk_end, total_read_count, hybridize,
                     chimeric_overlap, f_gtf, outdir, crl_file, tpm_threshold, score_cutoff,
                     n, sample_name, compress):
        """Wrapper function to process a chunk with shared objects.
        
        Args:
            shared_objects: Dictionary of shared objects (MPIRE) - passed as first arg by MPIRE
                - d_ref_lengths1: Dictionary of reference lengths for first priority reference
                - d_ref_lengths2: Dictionary of reference lengths for second priority reference
            chunk_start, chunk_end, total_read_count, hybridize, chimeric_overlap, f_gtf,
            outdir, crl_file, tpm_threshold, score_cutoff, n, sample_name, compress:
                Chunk processing parameters passed as positional arguments
        """
        return write_chimeras(
            chunk_start,
            chunk_end,
            total_read_count,
            None,  # d_ref_lengths1 - will use shared object
            None,  # d_ref_lengths2 - will use shared object
            hybridize,
            chimeric_overlap,
            f_gtf,
            outdir,
            crl_file,
            tpm_threshold,
            score_cutoff,
            n,
            sample_name,
            compress,
            shared_objects=shared_objects
        )
    
    # Determine start method: use 'fork' for copy-on-write on Unix systems, default on Windows
    # 'fork' provides copy-on-write semantics: shared objects only copied when modified
    # This is the most memory-efficient option when available (50-90% memory reduction)
    # Reference: https://sybrenjansen.github.io/mpire/master/usage/workerpool/shared_objects.html

    # Use 'fork' if available (Unix/Linux), otherwise use default (spawn/forkserver on Windows)
    # fork is not available on Windows, so we check platform and start method
    try:
        # Check if fork is available (Unix/Linux systems)
        if sys.platform != 'win32':
            # Try to get current start method - if it's already set to spawn, don't override
            current_method = multiprocessing.get_start_method(allow_none=True)
            if current_method != 'spawn':
                start_method = 'fork'
            else:
                start_method = None  # Use default (spawn)
        else:
            start_method = None  # Windows doesn't support fork
    except (AttributeError, ValueError):
        # Fallback if get_start_method is not available or raises error
        start_method = 'fork' if sys.platform != 'win32' else None
    
    with WorkerPool(n_jobs=args.processes, shared_objects=shared_objects_dict, 
                    start_method=start_method) as pool:
        pool.map(process_chunk, chunk_args, progress_bar=False)



def prepare_reference_file(args):
    """Prepare reference fasta file for hybridization."""
    if args.f_ref:
        return args.f_ref
    else:
        f_reference = os.path.join(args.outdir, 'merged_reference.fa')
        l_ref = [args.ref_fasta1]
        if args.ref_fasta2:
            l_ref.append(args.ref_fasta2)
        with open(f_reference, 'w', buffering=BUFFER_SIZE) as fh_out_ref:
            for fname in l_ref:
                with open(fname, 'r', buffering=BUFFER_SIZE) as infile:
                    for lin in infile:
                        fh_out_ref.write(lin)
        return f_reference



def extract_loci_sequences(args, f_reference):
    """Extract FASTA sequences for loci using bedtools getfasta."""
    for k in range(args.processes):
        bed = os.path.join(args.outdir, "loci.bed.") + str(k)
        fa = os.path.join(args.outdir, "loci.fa.") + str(k)
        getfasta_cmd = chira_utilities.get_bedtools_command('getfasta')
        process = subprocess.Popen(getfasta_cmd + ['-s', '-nameOnly',
                                    '-fi', f_reference,
                                    '-bed', bed,
                                    '-fo', fa],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT,
                                   universal_newlines=True)
        for l in sorted(set(process.stdout.readlines())):
            print(l, end="")



def build_intarna_params(args):
    """Build common IntaRNA parameters string."""
    noseed_param = ""
    if args.no_seed:
        noseed_param = "--noSeed"
    parts = ["--outMode C", "--outCsvCols id1,id2,start1,subseq1,start2,subseq2,hybridDPfull,E",
             noseed_param, "-m", args.intarna_mode, "--acc", args.accessibility, "--outPairwise=1",
             "--outNumber=1",
             "--temperature", str(args.temperature), "--seedBP", str(args.seed_bp),
             "--seedMinPu", str(args.seed_min_pu), "--accW", str(args.acc_width)]
    # Per-pair runs are small (1×1); use one thread so cluster can run more jobs in parallel.
    parts.append("--threads")
    parts.append("0")
    return " ".join(parts)



def run_hybridization(args):
    """Run hybridization: local IntaRNA if --use_batchtools is not set, else batchtools (LSF)."""
    use_batchtools = getattr(args, 'use_batchtools', False)
    chimeras_prefix = os.path.join(args.outdir, args.sample_name + ".chimeras-r")
    f_reference = prepare_reference_file(args)

    # Extract FASTA sequences for loci (needed for both local and batchtools)
    extract_loci_sequences(args, f_reference)

    common_intarna_params = build_intarna_params(args)

    if use_batchtools:
        batchtools_work_dir = os.path.abspath(os.path.join(args.outdir, "batchtools_work"))
        batchtools_registry = getattr(args, 'batchtools_registry', None) or os.path.abspath(os.path.join(batchtools_work_dir, "registry"))
        os.makedirs(batchtools_work_dir, exist_ok=True)
        run_batchtools_hybridization(args, common_intarna_params, batchtools_work_dir, batchtools_registry)
        if getattr(args, 'remove_intermediate', False):
            if os.path.exists(batchtools_work_dir):
                shutil.rmtree(batchtools_work_dir, ignore_errors=True)
    else:
        # Local hybridization: run IntaRNA per chunk via hybridize_and_write (MPIRE)
        remove_intermediate = getattr(args, 'remove_intermediate', False)
        local_hybrid_args = [(args.outdir, common_intarna_params, str(k), args.sample_name, remove_intermediate) for k in range(args.processes)]
        with WorkerPool(n_jobs=args.processes) as pool:
            pool.map(lambda t: hybridize_and_write(*t), local_hybrid_args, progress_bar=False)

    # Cleanup intermediate files when --remove_intermediate
    if getattr(args, 'remove_intermediate', False):
        for k in range(args.processes):
            for suf in ("loci.fa.", "loci.bed."):
                p = os.path.join(args.outdir, suf + str(k))
                if os.path.exists(p):
                    os.remove(p)
        if not args.f_ref:
            if os.path.exists(f_reference):
                os.remove(f_reference)
            if os.path.exists(f_reference + ".fai"):
                os.remove(f_reference + ".fai")

    return chimeras_prefix



def merge_output_files(args, chimeras_prefix):
    """Merge output files from multiple processes."""
    # File name prefixes
    if not args.hybridize:
        chimeras_prefix = os.path.join(args.outdir, args.sample_name + ".chimeras")
    
    chimeras_file = os.path.join(args.outdir, args.sample_name + ".chimeras.txt")
    singletons_file = os.path.join(args.outdir, args.sample_name + ".singletons.txt")
    singletons_prefix = os.path.join(args.outdir, args.sample_name + ".singletons")
    
    # Add .gz extension if compression is enabled
    if args.compress:
        chimeras_file += ".gz"
        singletons_file += ".gz"
    
    # Header fields
    header_chimeras = "\t".join(["read_id",
                                 "transcript_id_1",
                                 "transcript_id_2",
                                 "gene_id_1",
                                 "gene_id_2",
                                 "gene_symbol_1",
                                 "gene_symbol_2",
                                 "annotation_region_1",
                                 "annotation_region_2",
                                 "transcript_start_1",
                                 "transcript_end_1",
                                 "transcript_strand_1",
                                 "transcript_length_1",
                                 "transcript_start_2",
                                 "transcript_end_2",
                                 "transcript_strand_2",
                                 "transcript_length_2",
                                 "read_alignment_info",
                                 "genomic_coordinates_1",
                                 "genomic_coordinates_2",
                                 "locus_id_1",
                                 "locus_id_2",
                                 "crl_group_id_1",
                                 "crl_group_id_2",
                                 "tpm_1",
                                 "tpm_2",
                                 "alignment_score_1",
                                 "alignment_score_2",
                                 "combined_alignment_score",
                                 "hybridized_sequences",
                                 "hybridization_structure",
                                 "hybridization_positions",
                                 "hybridization_mfe_kcal_mol",
                                 "mirna_read_position",
                                 "hybridized_subsequences"])
    merge_files(chimeras_prefix, chimeras_file, header_chimeras, args.processes, args.compress, getattr(args, 'remove_intermediate', False))
    
    header_singletons = "\t".join(["read_id",
                                   "transcript_id",
                                   "gene_id",
                                   "gene_symbol",
                                   "annotation_region",
                                   "transcript_start",
                                   "transcript_end",
                                   "transcript_strand",
                                   "transcript_length",
                                   "read_alignment_info",
                                   "genomic_coordinates",
                                   "locus_id",
                                   "crl_group_id",
                                   "tpm",
                                   "alignment_score"])
    merge_files(singletons_prefix, singletons_file, header_singletons, args.processes, args.compress, getattr(args, 'remove_intermediate', False))



def parse_arguments():
    """Parse command-line arguments. Order: required I/O, references, parallelism, filtering, output, hybridization, batchtools, version."""
    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: extract chimeras',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required I/O
    parser.add_argument('-l', '--loci', action='store', dest='crl_file', required=True,
                        metavar='', help='Input BED file with alignments (e.g. loci.counts from chira_quantify)')
    parser.add_argument('-o', '--out', action='store', dest='outdir', required=True,
                        metavar='', help='Path to output directory')
    parser.add_argument('-n', '--sample_name', action='store', dest='sample_name', required=True,
                        metavar='', help='Sample name prefix for output files')
    
    # Reference FASTA
    parser.add_argument('-f1', '--ref_fasta1', action='store', dest='ref_fasta1', required=True,
                        metavar='', help='First priority reference FASTA file')
    parser.add_argument('-g', '--gtf', action='store', dest='f_gtf', required=False,
                        metavar='', help='Annotation GTF file')
    parser.add_argument('-f2', '--ref_fasta2', action='store', dest='ref_fasta2', required=False,
                        metavar='', help='Second priority reference FASTA file (e.g. miRNA)')
    parser.add_argument('-f', '--ref', action='store', dest='f_ref', required=False,
                        metavar='', help='Reference genomic FASTA (for IntaRNA accessibility)')

    # main job parallelism
    parser.add_argument('-p', '--processes', action='store', type=int, default=8, metavar='',
                        dest='processes',
                        help="""Number of processes to use for chimera extraction, hybridization prep/finish, and merging.""")

    # Filtering
    parser.add_argument('-tc', '--tpm_cutoff', action='store', type=chira_utilities.score_float, default=0, metavar='',
                        dest='tpm_cutoff',
                        help="""Transcripts with less than this percentile TPMs will be discarded in
the final output. [0, 1) - values >= 1.0 will be clamped to just below 1.0.""")
    parser.add_argument('-sc', '--score_cutoff', action='store', type=chira_utilities.score_float, default=0.0, metavar='',
                        dest='score_cutoff',
                        help='Hybrids with less than this score will be discarded in the final output. [0-1.0)')
    parser.add_argument('-co', '--chimeric_overlap', action='store', type=int, default=2, metavar='',
                        dest='chimeric_overlap',
                        help='Maximum number of bases allowed between the chimeric segments of a read')

    # Output options
    parser.add_argument("-s", '--summarize', action='store_true', dest='summarize',
                        help="Summarize interactions at loci level")
    parser.add_argument('-z', '--gzip', action='store_true', dest='compress',
                        help='Compress output files (chimeras and singletons) with gzip')

    # CRL index (.idx)
    parser.add_argument('--rebuild_index', action='store_true', dest='rebuild_index',
                        help='Rebuild CRL index even if .idx exists (default: use existing .idx when present)')
    parser.add_argument('--remove_index', action='store_true', dest='remove_index',
                        help='Remove the CRL index file (.idx) after processing (default: keep it)')

    # Hybridization (IntaRNA)
    parser.add_argument("-r", '--hybridize', action='store_true', dest='hybridize',
                        help="Run IntaRNA to hybridize the predicted chimeras")
    parser.add_argument("-ns", '--no_seed', action='store_true', dest='no_seed',
                        help="Do not enforce seed interactions")
    parser.add_argument("-acc", '--accessibility', type=str, choices=["C", "N"], default='N', required=False,
                        dest='accessibility', metavar='', help='IntaRNA accessibility: C (compute) or N (not)')
    parser.add_argument("-m", '--intarna_mode', type=str, choices=["H", "M", "S"], default='H', required=False,
                        dest='intarna_mode', metavar='', help='IntaRNA mode: H (heuristic), M (exact), S (seed-only)')
    parser.add_argument('-t', '--temperature', action='store', type=float, default=37, metavar='',
                        dest='temperature',
                        help='IntaRNA temperature parameter in Celsius to setup the VRNA energy parameters')
    parser.add_argument('-sbp', '--seed_bp', action='store', type=int, default=5, metavar='',
                        dest='seed_bp', choices=range(2, 20),
                        help='IntaRNA --seedBP parameter: number of inter-molecular base pairs within the seed region')
    parser.add_argument('-smpu', '--seed_min_pu', action='store', type=chira_utilities.score_float, default=0,
                        metavar='', dest='seed_min_pu',
                        help="""IntaRNA --seedMinPu parameter: minimal unpaired probability
(per sequence) a seed region may have.""")
    parser.add_argument('-accw', '--acc_width', action='store', type=int, default=150, metavar='',
                        dest='acc_width', choices=range(0, 99999),
                        help='IntaRNA --accW parameter: sliding window size for accessibility computation')

    # Batchtools (HPC cluster for IntaRNA)
    parser.add_argument('--use_batchtools', action='store_true', dest='use_batchtools',
                        help='Use R batchtools to submit IntaRNA jobs to HPC cluster. Requires --hybridize, R with batchtools and IntaRNA on cluster PATH.')
    parser.add_argument('--remove_intermediate', action='store_true', dest='remove_intermediate',
                        help='Remove intermediate files after success: loci.fa.<n>, loci.bed.<n>, merged reference (if not -f), work dirs (batchtools_work/ or local_intarna_<n>/), *.chimeras.<n>, and merge inputs (chimeras-r.<n>, singletons.<n>). Default: keep them.')
    parser.add_argument('--batchtools_registry', action='store', dest='batchtools_registry', default=None, metavar='',
                        help='Directory for batchtools registry (default: <outdir>/batchtools_work/registry)')
    parser.add_argument('--batchtools_template', action='store', dest='batchtools_template', default='', metavar='',
                        help='Path to batchtools LSF template (or "lsf-simple"). Default: lsf_custom.tmpl if present.')
    parser.add_argument('--batchtools_queue', action='store', dest='batchtools_queue', default='long', metavar='',
                        help='LSF queue name for batchtools jobs (default: long)')
    parser.add_argument('--batchtools_cores', action='store', type=int, default=1, metavar='',
                        dest='batchtools_cores',
                        help='Cores per LSF job (default: 1). IntaRNA runs with 1 thread; use 1 to run more jobs in parallel.')
    parser.add_argument('--batchtools_memory', action='store', dest='batchtools_memory', default=None, metavar='',
                        help='Total memory per job (e.g. 8GB, 64GB). Converted to per-core for LSF.')
    parser.add_argument('--batchtools_walltime', action='store', dest='batchtools_walltime', default='48:00', metavar='',
                        help='Walltime per job (e.g. 48:00 or 240:00)')
    parser.add_argument('--batchtools_conda_env', action='store', dest='batchtools_conda_env', default=None, metavar='',
                        help='Conda environment path for cluster jobs (optional)')
    parser.add_argument('--batchtools_max_parallel', action='store', type=int, default=None, metavar='',
                        dest='batchtools_max_parallel',
                        help='Max concurrent batchtools jobs (default: None = all chunks at once for minimum walltime)')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {chira_utilities.__version__}')

    return parser.parse_args()



def main():
    """Main function to orchestrate the chimera extraction workflow."""
    args = parse_arguments()
    # Use absolute paths so batchtools jobs (which write to absolute paths in jobs.json) and
    # the finish phase find the same result.csv files, regardless of the directory chira_extract was run from.
    args.outdir = os.path.abspath(args.outdir)
    args.crl_file = os.path.abspath(args.crl_file)
    print_configuration(args)
    validate_arguments(args)

    # Parse annotations if GTF file provided
    if args.f_gtf:
        print("Parsing the annotation file")
        parse_annotations(args.f_gtf)

    # Setup references
    d_reflen1, d_reflen2 = setup_references(args)

    # Parse CRLs file
    print("Parsing CRLs file")
    no_of_reads, tpm_cutoff_value = parse_counts_file(args.crl_file, args.tpm_cutoff, rebuild_index=args.rebuild_index)
    print("Done")

    # Run chimera extraction
    run_chimera_extraction(args, d_reflen1, d_reflen2, tpm_cutoff_value, no_of_reads)

    # Run hybridization if enabled (returns chimeras prefix)
    chimeras_prefix = None
    if args.hybridize:
        chimeras_prefix = run_hybridization(args)

    # Merge output files
    merge_output_files(args, chimeras_prefix)
    
    print(str(datetime.datetime.now()), " END: multiprocessing")
    
    # Write interaction summary if requested
    if args.summarize:
        write_interaction_summary(args.outdir, args.sample_name, args.compress, args.processes)
    
    # Cleanup: Remove index file only if requested
    if args.remove_index:
        index_file = args.crl_file + '.idx'
        if os.path.exists(index_file):
            try:
                os.remove(index_file)
                print(f"Removed index file: {index_file}", file=sys.stderr)
            except OSError as e:
                print(f"Warning: Could not remove index file {index_file}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()


