""" A Sequence class for biosequences """

from __future__ import print_function
import logging
import types
import numpy as np
from itertools import izip_longest

__all__ = ['Sequence', 'DnaSequence', 'RnaSequence', 'ProteinSequence']
logger = logging.getLogger(__name__)

class Sequence(object):
    """
    A biological sequence.
    """ 
    def __init__(self, sequence, alphabet=None):
        """
        Parameters
        ----------
        sequence : str or pybio Sequence
            string with the characters in the Sequence
        alphabet : list of chars, default None
            the characters in this Sequence
        """
        self.identifier = None
        self.description = None
        self._max_output_chars = 40
        
        initializers = {
                types.NoneType : self._init_empty,
                str : self._init_from_str,
                list : self._init_from_str,
                Sequence : self._init_from_sequence
                }

        for cls in type(sequence).mro():
            if cls in initializers:
                initializers[cls](sequence, alphabet)
                break

        if not hasattr(self, '_sequence'):
            raise ValueError('Sequence cannot be initialized from {sequence_type}'.format(sequence_type=type(sequence).__name__))

    @classmethod
    def _from_ndarray(cls, array, alphabet):
        sequence = cls(None)
        sequence._sequence = array
        sequence._alphabet = alphabet
        return sequence

    @property
    def alphabet(self):
        """Get the sequence alphabet."""
        return self._alphabet

    def to_list(self):
        """Returns the sequence as a list of characters."""
        return [self._alphabet[value] for value in self._sequence]
        
    def to_string(self):
        """Returns the sequence as a string."""
        return ''.join(self.to_list())

    def _init_empty(self, sequence, alphabet):
        self._sequence = None
        self._alphabet = list(alphabet)

    def _init_from_str(self, sequence, alphabet):
        self._alphabet = list(alphabet) if alphabet else list(set(sequence))
        self._validate_str_sequence(sequence, self._alphabet)
        inverse_alphabet = {letter:value for value,letter in enumerate(self._alphabet)}
        self._sequence = np.array([inverse_alphabet[letter] for letter in sequence], dtype=np.int8)

    def _init_from_sequence(self, sequence, alphabet):
        if alphabet != sequence.alphabet:
            raise ValueError('The alphabet argument does not match the sequence alphabet.')
        self._alphabet = sequence.alphabet
        self._sequence = sequence._sequence

    def _validate_str_sequence(self, sequence, alphabet):
        for index, char in enumerate(sequence):
            if char not in alphabet:
                raise ValueError('Sequence contains char {char} at position {index} that is not in the alphabet'.format(char=char,index=index))

    def __getitem__(self, slc):
        return self._from_ndarray(self._sequence[slc], self._alphabet)

    def __str__(self):
        if len(self._sequence) < self._max_output_chars:
            return self.to_string()
        else:
            return self[0:self._max_output_chars/2].to_string() + '...' + self[-self._max_output_chars/2:].to_string()

    def __repr__(self):
        return self.__class__.__name__ + '(' + self.__str__() + ')'

    def __len__(self):
        return len(self._sequence)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.to_list() == other.to_list()
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

class DnaSequence(Sequence):
    """
    A DNA sequence.
    """
    _dna_alphabet = 'ACGTN'
    _transcription = {'T':'U'}
    _complement = {'A':'T','T':'A','G':'C','C':'G'}

    def __init__(self, sequence):
        """
        Parameters
        ----------
        sequence : str or pybio Sequence
            the characters in the DNA Sequence, must be in the DNA alphabet
        """
        super(DnaSequence, self).__init__(sequence, DnaSequence._dna_alphabet)

    def transcribe(self):
        """Returns the transcription to an RNA sequence."""
        rna_alphabet = _convert_alphabet(self._alphabet, DnaSequence._transcription)
        return RnaSequence._from_ndarray(self._sequence, rna_alphabet)

    def complement(self):
        """Returns the complement of the sequence."""
        complement_sequence = _convert_sequence(self._sequence, self._alphabet, self._complement)
        return DnaSequence._from_ndarray(complement_sequence, self._alphabet)

    def reverse_complement(self):
        """Returns the reverse complement of the sequence."""
        return self.complement()[::-1]

    def translate(self):
        """Returns the translated protein sequence."""
        rna = self.transcribe()
        return rna.translate()

class RnaSequence(Sequence):
    """
    An RNA sequence.
    """
    _rna_alphabet = 'ACGUN'
    _reverse_transcription = {'U':'T'}
    _complement = {'A':'U','U':'A','G':'C','C':'G'}  
    _translation = {'UUU':'F','UUC':'F','UUA':'L','UUG':'L',
                    'CUU':'L','CUC':'L','CUA':'L','CUG':'L',
                    'AUU':'I','AUC':'I','AUA':'I','AUG':'M',
                    'GUU':'V','GUC':'V','GUA':'V','GUG':'V',
                    'UCU':'S','UCC':'S','UCA':'S','UCG':'S',
                    'CCU':'P','CCC':'P','CCA':'P','CCG':'P',
                    'ACU':'T','ACC':'T','ACA':'T','ACG':'T',
                    'GCU':'A','GCC':'A','GCA':'A','GCG':'A',
                    'UAU':'Y','UAC':'Y','UAA':'Stop','UAG':'Stop',
                    'CAU':'H','CAC':'H','CAA':'Q','CAG':'Q',
                    'AAU':'N','AAC':'N','AAA':'K','AAG':'K',
                    'GAU':'D','GAC':'D','GAA':'E','GAG':'E',
                    'UGU':'C','UGC':'C','UGA':'Stop','UGG':'W',
                    'CGU':'R','CGC':'R','CGA':'R','CGG':'R',
                    'AGU':'S','AGC':'S','AGA':'R','AGG':'R',
                    'GGU':'G','GGC':'G','GGA':'G','GGG':'G'}
    
    def __init__(self, sequence):
        """
        Parameters
        ----------
        sequence : str or pybio Sequence
            the characters in the DNA Sequence, must be in the RNA alphabet
        """
        super(RnaSequence, self).__init__(sequence, RnaSequence._rna_alphabet)

    def reverse_transcribe(self):
        """Returns the reverse transcription to a DNA sequence."""
        rna_alphabet = _convert_alphabet(self._alphabet, RnaSequence._reverse_transcription)
        return DnaSequence._from_ndarray(self._sequence, rna_alphabet)

    def complement(self):
        """Returns the complement of the sequence."""
        complement_sequence = _convert_sequence(self._sequence, self._alphabet, self._complement)
        return RnaSequence._from_ndarray(complement_sequence, self._alphabet)

    def reverse_complement(self):
        """Returns the reverse complement of the sequence."""
        return self.complement()[::-1]

    def translate(self):
        """Returns the translated protein sequence, up to the first stop codon"""
        #TODO: Use Numba and test speed
        offset = len(self._sequence) % 3
        if offset % 3 != 0:
            logger.warning('The length of this sequence is not a multiple of three...')

        def build_codon(triplet):
            return ''.join([self._rna_alphabet[base] for base in triplet])

        protein_sequence = []
        for triplet in _grouper(self._sequence[:len(self._sequence)-offset], 3):
            amino_acid = self._translation[build_codon(triplet)]
            if amino_acid == 'Stop':
                break
            protein_sequence.append(amino_acid)

        return ProteinSequence(protein_sequence)

class ProteinSequence(Sequence):
    """
    A protein sequence.
    """
    _protein_alphabet = 'ARNDCQEGHILKMFPSTWYVBZX*'

    def __init__(self, sequence):
        """
        Parameters
        ----------
        sequence : str or pybio Sequence
            the characters in the protein sequence, must be in the protein amino acid alphabet
        """
        super(ProteinSequence, self).__init__(sequence, ProteinSequence._protein_alphabet)

def _convert_sequence(sequence, alphabet, conversion):
    #TODO: Documentation and Numba
    letter_to_index = {letter:index for index, letter in enumerate(alphabet)}
    index_list = np.array([letter_to_index[conversion[letter]] if letter in conversion else letter_to_index[letter] for letter in alphabet], dtype=np.int8)
    return index_list[sequence]

def _convert_alphabet(alphabet, conversion):
    """
    Returns a new alphabet with items converted according to the conversion dictionary.

    Parameters
    ----------
    alphabet : str or list
        the alphabet to convert
    conversion : dict
        a dictionary to convert the alphabet
    """
    return [conversion[letter] if letter in conversion else letter for letter in alphabet]

def _grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks, from the itertools recipe book"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)
