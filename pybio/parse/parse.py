"""Classes and methods for parsing files"""

from __future__ import print_function
import logging
from collections import namedtuple
from itertools import izip_longest

__all__ = ['Fasta', 'Fastq']
logger = logging.getLogger(__name__)

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

class Fasta(namedtuple('Fasta', ['identifier', 'description', 'sequence'])):
    """
    A representation of the data in a parsed FASTA entry.

    identifier - an identifier for a sequence often following the NCBI FASTA defline format
    description - a description of the sequence
    sequence - the amino acid or nucleoacid sequence in the FASTA format convention

    The Fasta class has methods for parsing and writing to FASTA files.

    Class methods
    -------------
    parse_string
        parses a Fasta object from a FASTA-format string
    parse_iterator
        parses a FASTA file and returns an iterator of Fasta objects
    """
    @staticmethod
    def _parse_description_line(description_line):
        description_line = description_line[1:].strip().split(' ', 1)
        if len(description_line) is 2:
            return description_line
        elif description_line:
            return description_line, None
        else:
            return None, None
        
    @classmethod
    def parse_string(cls, fasta_str):
        if not fasta_str.startswith('>'):
            raise ValueError('Entry is not in FASTA format.')
        description_line, sequence = fasta_str.split('\n', 1)
        identifier, description = cls._parse_description_line(description_line)
        sequence = ''.join(sequence.split())
        return cls(identifier, description, sequence)

    @classmethod
    def parse_iterator(cls, fasta_file):
        with open(fasta_file, 'r') as fasta:
            sequence = ''
            for line in fasta:
                if not line.strip():
                    continue
                if line.startswith('>'):
                    if sequence:
                        yield Fasta(identifier, description, sequence)
                    sequence = ''
                    identifier, description = cls._parse_description_line(line)
                else:
                    sequence += ''.join(line.split())
            if sequence:
                yield Fasta(identifier, description, sequence)
    
class Fastq(namedtuple('Fastq', ['identifier', 'description', 'sequence', 'quality_scores'])):
    """
    A representation of the data in a parsed FASTQ entry.

    The Fastq class has methods for parsing and writing to FASTQ files.

    identifier - an identifier for a sequence often following the NCBI FASTA defline format
    description - a description of the sequence
    sequence - the amino acid or nucleoacid sequence in the FASTA format convention
    quality_scores - Phred quality scores for each item in the sequence. The sequence is encoded as
    ASCII characters from '!' to '~'.

    Phred quality scores give the probability of the base call being incorrect, according to
    the formula P = 10^(-Q/10), where Q is the Phred quality score.

    Class methods
    -------------
    parse_string
        parses a Fastq object from a FASTQ-format string
    parse_iterator
        parses a FASTQ file and returns an iterator of Fastq objects
    """
    @staticmethod
    def _parse_description_line(description_line):
        description_line = description_line[1:].strip().split(' ', 1)
        if len(description_line) is 2:
            return description_line
        elif description_line:
            return description_line, None
        else:
            return None, None
    
    @classmethod
    def parse_entry(cls, entry):
        description_line, sequence, optional_description_line, quality_scores = entry
        identifier, description, = cls._parse_description_line(description_line)
        return cls(identifier, description, sequence.strip(), quality_scores.strip())

    @classmethod
    def parse_string(cls, fastq_str):
        if not fastq_str.startswith('@'):
            raise ValueError('Entry is not in FASTQ format.')
        description_line, sequence, optional_description_line, quality_scores = fastq_str.split('\n', 3)
        identifier, description, = cls._parse_description_line(description_line)
        return cls(identifier, description, sequence.strip(), quality_scores.strip())

    @classmethod
    def parse_iterator(cls, fastq_file):
        with open(fastq_file, 'r') as fastq:
            for entry in grouper(fastq, 4):
                yield cls.parse_entry(entry)

