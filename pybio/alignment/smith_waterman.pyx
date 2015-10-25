""" Functions for performing Smith-Waterman """

#from __future__ import print_function

cimport numpy as cnp
cimport ssw
from collections import namedtuple
from pybio import DnaSequence, RnaSequence, ProteinSequence
from matrix import blastn_default_matrix, blosum62
from alignment import Alignment, CigarSegment, Cigar

__all__ = ['smith_waterman']

_default_args = {
        DnaSequence : {'match_matrix':blastn_default_matrix, 'gap_open_penalty':5, 'gap_extend_penalty':2},
        RnaSequence : {'match_matrix':blastn_default_matrix, 'gap_open_penalty':5, 'gap_extend_penalty':2},
        ProteinSequence : {'match_matrix':blosum62, 'gap_open_penalty':11, 'gap_extend_penalty':1},
        }

def smith_waterman(query, target, **kwargs):
    """
    Perform a Smith-Waterman alignment on two sequences.

    Parameters
    ----------
    query : Sequence
        The query sequence.
    target : Sequence
        The target sequence of matching type with the query Sequence.
    match_matrix : numpy array
        A square numpy array with dimensions equal to the size of the Sequence alphabet.
        This matrix gives the score for every match or mismatch.
        Defaults to an identity matrix; ones for match, zero for no match.
    gap_open_penalty : np.uint8
        Penalty for creating a gap in the alignment.
    gap_extend_penalty : np.uint8
        Penalty for extending a gap in the alignment.
    score_filter : np.uint16
        Alignment score must be greater than the score_filter to return an alignment.
        Defaults to 0.
    distance_filter : np.uint32
        The target and query sequence must have a length greater than the distance_filter
        to return an alignment.
        Defaults to 0.
    score_size : np.unint8 in [0,1,2]
        If your estimated best alignment score is < 255 this should be 0.
        If your estimated best alignment score is >= 255, this should be 1.
        If you don't know, this should be 2.
        Defaults to 2.

    Match Defaults
    --------------
    DnaSequence and RnaSequence defaults match BLASTN defaults:
        match_matrix = blastn_default_matrix
        gap_open_penalty = 5
        gap_extend_penalty = 2
    ProteinSequence defaults match BLASTP defaults:
        match_matrix = blosum62
        gap_open_penalty = 11
        gap_extension_penalty = 1
    General defaults:
        match_matrix = identity
        gap_open_penalty = 5
        gap_extend_penalty = 2

    Notes
    -----
    This implementation uses the SSW Library, written by Mengyao Zhao and Wan-Ping Lee. SSW is 
    a striped implementation that uses Single-Instruction Multiple-Data (SIMD) to parallelize 
    on a single processor. The implementation is ~50 times faster than a naive Smith-Waterman.

    References
    ----------
    Zhao, Mengyao, Wan-Ping Lee, Erik P. Garrison, & Gabor T. Marth. 
    "SSW Library: An SIMD Smith-Waterman C/C++ Library for Applications". PLOS ONE (2013). Web. 11 July 2014.
    http://www.plosone.org/article/info:doi/10.1371/journal.pone.0082138
    """

    smith_waterman_args = _default_args[type(query)].copy() if type(query) in _default_args else {}
    for arg, value in kwargs.iteritems():
        smith_waterman_args[arg] = value

    if smith_waterman_args['match_matrix'] is None:
        smith_waterman_args = cnp.identity(len(query.alphabet))

    return _smith_waterman(query, target, **smith_waterman_args)
    
def _smith_waterman(query, target, cnp.ndarray match_matrix, cnp.uint8_t gap_open_penalty=5, cnp.uint8_t gap_extend_penalty=2, cnp.uint16_t score_filter=0, cnp.uint32_t distance_filter=0, cnp.int8_t score_size=2, ):
    """
    Perform a Smith-Waterman alignment on two sequences.

    Parameters
    ----------
    query : Sequence
        the query sequence
    target : Sequence
        the target sequence of matching type with the query Sequence
    match_matrix : numpy array
        a flat numpy array with dimensions equal to the size of the Sequence alphabet
        This matrix gives the score for every match or mismatch.
    gap_open_penalty : np.uint8
        penalty for creating a gap in the alignment
    gap_extend_penalty : np.uint8
        penalty for extending a gap in the alignment
    score_filter : np.uint16
        alignment score must be greater than the score_filter to return an alignment
    distance_filter : np.uint32
        the target and query sequence must have a length greater than the distance_filter
        to return an alignment
    score_size : np.unint8 in [0,1,2]
        If your estimated best alignment score is < 255 this should be 0.
        If your estimated best alignment score is >= 255, this should be 1.
        If you don't know, this should be 2.
    """

    cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] query_sequence = query._sequence
    cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] target_sequence = target._sequence
    cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] matrix = match_matrix.flatten()

    cdef cnp.uint8_t bit_flag = 0
    bit_flag = bit_flag | 0x2

    #if score_filter != 0:
    #    bit_flag = bit_flag | 0x2
    #if distance_filter != 0:
    #    bit_flag = bit_flag | 0x4
    #if bit_flag == 0 or bit_flag == 8:
    #    bit_flag = bit_flag | 0x1

    cdef cnp.uint32_t mask_length = max(len(query_sequence)/2, 15)

    cdef ssw.s_profile* profile = ssw.ssw_init(<cnp.int8_t*> query_sequence.data, len(query_sequence), <cnp.int8_t*> matrix.data, len(query.alphabet), score_size)
    cdef ssw.s_align* align = ssw.ssw_align(profile, <cnp.int8_t*> target_sequence.data, len(target_sequence), gap_open_penalty, gap_extend_penalty, bit_flag, score_filter, distance_filter, mask_length)

    cigar = get_cigar(align)

    alignment = Alignment(query, target, align.score1, align.score2, cigar, align.ref_begin1, align.ref_end1, align.read_begin1, align.read_end1)

    ssw.init_destroy(profile)
    ssw.align_destroy(align)

    return alignment

cdef get_cigar(ssw.s_align* align):
    sam_op_table = ['M', 'I', 'D']
    cigar = []
    for i in range(align.cigarLen):
        segment_length = align.cigar[i] >> 4
        sam_op = sam_op_table[align.cigar[i] & 0xf]
        cigar.append(CigarSegment(segment_length, sam_op))
    return cigar

