from collections import namedtuple
from StringIO import StringIO
from itertools import izip

__all__ = ['Alignment', 'CigarSegment']

class Alignment(object):
    """
    An alignment of two sequences.
    """
    def __init__(self, query, target, optimal_alignment_score, suboptimal_alignment_score, cigar, begin_query, end_query, begin_target, end_target):
        """
        Parameters
        ----------
        query : Sequence
        target : Sequence
        cigar : list of CigarSegments
            A CIGAR, as defined in the SAM Format.
        optimal_alignment_score : int
            the optimal alignment score of the alignment
        suboptimal_alignment_score : int
            the suboptimal alignment score of the alignment
        begin_query : int
            the position at which the query sequence begins
        end_query : int
            the position at which the query sequence ends
        begin_target : int
            the position at which the target sequence begins
        end_target : int
            the position at which the target sequence ends
        """
        self.query = query
        self.target = target
        self.optimal_alignment_score = optimal_alignment_score
        self.suboptimal_alignment_score = suboptimal_alignment_score
        self.cigar = Cigar(cigar)
        self.begin_query = begin_query
        self.end_query = end_query
        self.begin_target = begin_target
        self.end_target = end_target
        self.query_str = self._seq_cig(self.query[self.begin_query:self.end_query+1], self.cigar, 'D')
        self.target_str = self._seq_cig(self.target[self.begin_target:self.end_target+1], self.cigar, 'I')
        self.alignment_str = self._calc_alignment_str()
    
    def _calc_alignment_str(self):
        def alignment_str_char(query_char, target_char):
            if query_char=='-' or target_char=='-':
                return ' '
            if query_char != target_char:
                return '*'
            return '|'
        return ''.join([alignment_str_char(query_char, target_char) for query_char, target_char in izip(self.query_str, self.target_str)])
            
    def __getitem__(self, slc):

        pass

    def _get_line(self, slc, digits):
        #query_line = 'query ' + str(self.begin_query + slc.start).rjust(digits) + '  ' + self.query_str[slc] + '  ' + str(self.end_query + slc.stop
        query_start = str(self.begin_query + slc.start).rjust(digits)
        query_stop = str(self.begin_query + slc.stop).rjust(digits)
        query_str = self.query_str[slc]
        alignment_str = self.alignment_str[slc]
        target_str = self.target_str[slc]
        space_digits = digits * ' '
        target_start = str(self.begin_target + slc.start).rjust(digits)
        target_stop = str(self.begin_target + slc.stop).rjust(digits)
        query_line = 'query  {query_start}  {query_str}  {query_stop}'
        align_line = '       {space_digits}  {alignment_str}'
        target_line = 'target {target_start}  {target_str}  {target_stop}'
        return '\n'.join([query_line, align_line, target_line]).format(**locals())

    def to_string(self):
        identifier_str = "target: {target}\n query: {query}\n".format(target=self.target.identifier, query=self.query.identifier) if self.target.identifier or self.query.identifier else ""
    
        #target_str = self._seq_cig(self.target[self.begin_target:self.end_target+1], self.cigar, 'I')
        #query_str = self._seq_cig(self.query[self.begin_query:self.end_query+1], self.cigar, 'D')

        score_str = "\noptimal alignment score: {optimal_alignment_score}, suboptimal alignment score: {suboptimal_alignment_score}".format(optimal_alignment_score=self.optimal_alignment_score, suboptimal_alignment_score=self.suboptimal_alignment_score)
        return "\n".join([identifier_str, self.target_str, self.alignment_str, self.query_str, score_str])

    def _seq_cig(self, seq, cigar, gap_op):
        sequence_io = StringIO(seq.to_string())
        alignment_str = ''
        for segment in cigar:
            alignment_str += '-' * segment.length if segment.op == gap_op else sequence_io.read(segment.length)
        return alignment_str

    def __str__(self):
        identifier_str = "target: {target}\n query: {query}\n".format(target=self.target.identifier, query=self.query.identifier) if self.target.identifier or self.query.identifier else ""
        #target_str = self._get_alignment_str(self.target, self.cigar, self.begin_target, self.begin_target - self.begin_query, 'I')
        #query_str = self._get_alignment_str(self.query, self.cigar, self.begin_query, self.begin_query - self.begin_target, 'D')
        #target_str = self._get_alignment_str(self.target, self.cigar, self.begin_query, self.begin_query - self.begin_target, 'I')
        #query_str = self._get_alignment_str(self.query, self.cigar, self.begin_target, self.begin_target - self.begin_query, 'D')
        score_str = "\noptimal alignment score: {optimal_alignment_score}, suboptimal alignment score: {suboptimal_alignment_score}".format(optimal_alignment_score=self.optimal_alignment_score, suboptimal_alignment_score=self.suboptimal_alignment_score)
        return "\n".join([identifier_str, self.target_str, self.query_str, score_str])
    
    def __repr__(self):
        return self.__class__.__name__ + '\n\n' + self.__str__()

    @staticmethod
    def _get_alignment_str(sequence, cigar, begin, offset, gap_op):
        if type(sequence) is str:
            sequence_io = StringIO(sequence)
        else:
            sequence_io = StringIO(sequence.to_string())
        alignment_str = ' ' * max(0,offset) + sequence_io.read(begin-offset)
        for segment in cigar:
            alignment_str += '-' * segment.length if segment.op == gap_op else sequence_io.read(segment.length)
        alignment_str += sequence_io.read()
        return alignment_str

class Cigar(list):
    """
    A list of CigarSegments, namedtuples with a length and op code.

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
    A named tuple to hold segments of a Cigar, which describes an alignment.
    The cigar is defined as part of the Sequence Alignment/Map (SAM) Format
    specification.

    The CigarSegment describes a segment of the alignment with a lenght and
    an op code defing the match state of that segment.

    M   alignment match
    I   insertion to the target
    D   deletion from the target
    N   skipped region from the target
    S   soft clipping
    H   hard clipping
    P   padding
    =   sequence match
    X   sequence mismatch
    """
    pass
