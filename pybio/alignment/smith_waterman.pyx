""" Functions for performing Smith-Waterman """

cimport numpy as cnp
import numpy as np
cimport ssw
from collections import namedtuple
from cStringIO import StringIO

mat50 = np.array([
   # A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
     5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1, -5, # A
    -2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1,  0, -1, -5, # R
    -1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  5,  0, -1, -5, # N
    -2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  6,  1, -1, -5, # D
    -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -3, -1, -5, # C
    -1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0,  4, -1, -5, # Q
    -1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1,  5, -1, -5, # E
     0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -2, -1, -5, # G
    -2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0,  0, -1, -5, # H
    -1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4, -3, -1, -5, # I
    -2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4, -3, -1, -5, # L
    -1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0,  1, -1, -5, # K
    -1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3, -1, -1, -5, # M
    -3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4, -4, -1, -5, # F
    -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -1, -1, -5, # P
     1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0,  0, -1, -5, # S
     0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -5, # T
    -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -1, -5, # W
    -2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -2, -1, -5, # Y
     0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -3, -3, -1, -5, # V
    -2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -3,  6,  1, -1, -5, # B
    -1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  1,  5, -1, -5, # Z
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5, # X
    -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1  # *
    ], dtype=np.int8)

dna_match = np.array([
#    A   C   G   T   N
     2, -3, -3, -3,  0, # A
    -3,  2, -3, -3,  0, # C
    -3, -3,  2, -3,  0, # G
    -3, -3, -3,  2,  0, # T
     0,  0,  0,  0,  0  # N
    ], dtype=np.int8)

sam_op_table = ['M', 'I', 'D']

amino_acid_sequence_order = 'ARNDCQEGHILKMFPSTWYVBZX*'

amino_acid_conversion_table = np.array([
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23,  0, 20,  4,  3,  6, 13,  7,  8,  9, 23, 11, 10, 12,  2, 23,
    14,  5,  1, 15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
    23,  0, 20,  4,  3,  6, 13,  7,  8,  9, 23, 11, 10, 12,  2, 23,
    14,  5,  1, 15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23])

dna_sequence_order = 'ACGTN'

nucleotide_conversion_table = np.array([
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  0,  4,  1,  4,  4,  4,  2,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  3,  0,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  0,  4,  1,  4,  4,  4,  2,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  3,  0,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4])

cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] convert_sequence(sequence, conversion_table=nucleotide_conversion_table):
    """ Converts the sequence into a int8 array according to a conversion table """
    cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] seq = np.array([conversion_table[ord(char)] for char in sequence], dtype=np.int8)
    return seq

def build_match_matrix(match_score=2, mismatch_score=-3, sequence_order=dna_sequence_order):
    """ Create a matrix with zeros for the outside column and row, and a diagonal matrix of mismatch score with match score on the diagonals for the rest. Return flattened. """
    n = len(sequence_order)
    match_matrix = np.zeros((n,n), dtype=np.int8)
    match_matrix[:n-1,:n-1].fill(mismatch_score)
    np.fill_diagonal(match_matrix[:n-1,:n-1], match_score)
    return match_matrix.flatten()

def smith_waterman_proteins(query, target, match_matrix=mat50, cnp.int8_t score_size=2, cnp.uint8_t gap_open_penalty=5, cnp.uint8_t gap_extend_penalty=2, cnp.uint16_t score_filter=0, cnp.uint32_t distance_filter=0, sequence_order=amino_acid_sequence_order, conversion_table=amino_acid_conversion_table):
    """
    Performs a Smith-Waterman sequence alignment on amino acid sequences
    
    Example usage:
    pybio.alignment.smith_waterman_proteins('AGRLQLQVIICVATRTAGYTEYG','GRLQLQVILQVIICVICVAAATRTAGYTE')
    """
    smith_waterman(query, target, match_matrix, score_size, gap_open_penalty, gap_extend_penalty, score_filter, distance_filter, sequence_order, conversion_table)

def smith_waterman(query_sequence, target_sequence, match_matrix=dna_match, cnp.int8_t score_size=2, cnp.uint8_t gap_open_penalty=5, cnp.uint8_t gap_extend_penalty=2, cnp.uint16_t score_filter=0, cnp.uint32_t distance_filter=0, sequence_order=dna_sequence_order, conversion_table=nucleotide_conversion_table):
    """
    Performs a Smith-Waterman sequence alignment.

    Example usage:
    pybio.alignment.smith_waterman('ACTAGAATGGCT','CCATACTGAACTGACTAAC', gap_open_penalty=2)

    Parameters
    ----------
    query_sequence : str
        sequence to query
    target_sequence : str
        sequence to align the query against
    match_matrix : list
        a flattened matrix of alignment scores
    sequence_order : str
        a string listing the symbols in the query/target alphabet in an order matching the alignment matrix
    """

    # Create profile for query sequence

    cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] query = convert_sequence(query_sequence, conversion_table)
    cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] matrix = match_matrix

    cdef ssw.s_profile* profile
    profile = ssw.ssw_init(<cnp.int8_t*> query.data, len(query),
                           <cnp.int8_t*> matrix.data, len(sequence_order), score_size)
 
    # Do the match

    cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] reference
    reference = convert_sequence(target_sequence, conversion_table)

    cdef cnp.uint8_t bit_flag
    bit_flag = 0
    if score_filter != 0:
        bit_flag = bit_flag | 0x2
    if distance_filter != 0:
        bit_flag = bit_flag | 0x4
    if bit_flag == 0 or bit_flag == 8:
        bit_flag = bit_flag | 0x1

    cdef cnp.uint32_t mask_length = max(len(query_sequence)/2, 15)

    cdef ssw.s_align* align
    align = ssw.ssw_align(profile, <cnp.int8_t*> reference.data, len(target_sequence), 
                          gap_open_penalty, gap_extend_penalty, bit_flag,
                          score_filter, distance_filter, mask_length)

    cigar = get_cigar_segments(align)

    print 'cigar',cigar_to_str(cigar)

    print alignment_as_str(query_sequence, target_sequence, cigar, align.ref_begin1, align.read_begin1)

    ssw.init_destroy(profile)
    ssw.align_destroy(align)

cdef get_cigar_segments(ssw.s_align* align):
    CigarSegment = namedtuple('CigarSegment', 'length op')
    cigar_segments = []
    for i in range(align.cigarLen):
        segment_length = align.cigar[i] >> 4
        sam_op = sam_op_table[align.cigar[i] & 0xf]
        cigar_segments.append(CigarSegment(segment_length, sam_op))
    return cigar_segments

def cigar_to_str(cigar):
    return "".join([str(segment.length) + segment.op for segment in cigar])

def get_alignment_str(sequence, cigar, begin, gap_op):
    sequence_io = StringIO(sequence)
    alignment_str = ' ' * begin + sequence_io.read(max(0,-begin))
    for segment in cigar:
        alignment_str += '-' * segment.length if segment.op == gap_op else sequence_io.read(segment.length)
    alignment_str += sequence_io.read()
    return alignment_str

def alignment_as_str(query, target, cigar, begin_query, begin_target):
    query_str = get_alignment_str(query, cigar, begin_query - begin_target, 'D')
    target_str = get_alignment_str(target, cigar, begin_target - begin_query, 'I')
    return "\n".join([target_str, query_str])


