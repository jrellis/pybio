from collections import namedtuple
from StringIO import StringIO

__all__ = ['Alignment', 'CigarSegment']

class Alignment(object):
    """
    An alignment of two sequences.
    """
    def __init__(self, query, target, cigar, begin_query, begin_target):
        """
        Parameters
        ----------
        query : Sequence
        target : Sequence
        cigar : list of CigarSegments
            A CIGAR, as defined in the SAM Format.
        begin_query : int
            the position at which the query sequence begins
        begin_target : int
            the position at which the target sequence begins
        """
        self.query = query
        self.target = target
        self.cigar = Cigar(cigar)
        self.begin_query = begin_query
        self.begin_target = begin_target

    def __str__(self):
        query_str = self._get_alignment_str(self.query, self.cigar, self.begin_query - self.begin_target, 'D')
        target_str = self._get_alignment_str(self.target, self.cigar, self.begin_target - self.begin_query, 'I')
        return "\n".join([target_str, query_str])
    
    def __repr__(self):
        return self.__class__.__name__ + '\n' + self.__str__()

    @staticmethod
    def _get_alignment_str(sequence, cigar, begin, gap_op):
        if type(sequence) is str:
            sequence_io = StringIO(sequence)
        else:
            sequence_io = StringIO(sequence.to_string())
        alignment_str = ' ' * begin + sequence_io.read(max(0,-begin))
        for segment in cigar:
            alignment_str += '-' * segment.length if segment.op == gap_op else sequence_io.read(segment.length)
        alignment_str += sequence_io.read()
        return alignment_str

class Cigar(list):
    """
    A list of CigarSegments, named tuples with a length and op code.

    The cigar is defined as part of the Sequence Alignment/Map (SAM) Format
    specification.

    A cigar describes the an alignment in pairs, a length and an op code.
    The op code describes the alignment state of a query segment with a
    target segment.

    M   alignment match
    I   insertion to the target
    D   deletion from the target
    N   skipped region from the target
    S   soft clipping
    H   hard clipping
    P   padding
    =   sequence match
    X   sequence mismatch

    References
    ----------
    1. Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer,Gabor Marth, Goncalo Abecasis, Richard Durbin, and 1000 Genome Project Data Processing Subgroup
    "The Sequence Alignment/Map format and SAMtools". Bioinformatics 2009 Aug 15; 25(16): 2078-2079
    2. The SAM/BAM Format Specification Working Group 
    "Sequence Alignment/Map Format Specification". https://samtools.github.io/hts-specs/SAMv1.pdf
    """
    def to_string(self):
        return "".join([str(segment.length) + segment.op for segment in self])
    def __str__(self):
        return self.to_string()
    def __repr__(self):
        return self.__class__.__name__ + '(' + self.__str__() + ')'

class CigarSegment(namedtuple('CigarSegment',['length', 'op'])):
    """
    A named tuple to hold segments of a cigar, which describes an alignment.

    The cigar is defined as part of the Sequence Alignment/Map (SAM) Format
    specification.

    A cigar describes the an alignment in pairs, a length and an op code.
    The op code describes the alignment state of a query segment with a
    target segment.

    M   alignment match
    I   insertion to the target
    D   deletion from the target
    N   skipped region from the target
    S   soft clipping
    H   hard clipping
    P   padding
    =   sequence match
    X   sequence mismatch

    References
    ----------
    1. Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer,Gabor Marth, Goncalo Abecasis, Richard Durbin, and 1000 Genome Project Data Processing Subgroup
    "The Sequence Alignment/Map format and SAMtools". Bioinformatics 2009 Aug 15; 25(16): 2078-2079
    2. The SAM/BAM Format Specification Working Group 
    "Sequence Alignment/Map Format Specification". https://samtools.github.io/hts-specs/SAMv1.pdf
    """
    pass
