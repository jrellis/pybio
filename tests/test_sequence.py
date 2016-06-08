import pybio
import pytest
    
seq_string = 'mwklrxinoidoqxlakbrp'

dna_string = 'ACGTTGCAC'
dna_bad = 'ACGTTGCACQ'
dna_complement = 'TGCAACGTG'
dna_reverse_complement = 'GTGCAACGT'
dna_transcribe = 'ACGUUGCAC'
dna_translate = 'TLH'

rna_string = 'ACGUUGCAC'
rna_bad = 'ACGTTGCAC'
rna_complement = 'UGCAACGUG'
rna_reverse_complement = 'GUGCAACGU'
rna_reverse_transcribe = 'ACGTTGCAC'
rna_translate = 'TLH'

protein_string = 'TLH'
protein_bad = '.'

class TestSequence:

    def test_init_random_string(self):
        pybio.Sequence(seq_string)

    def test_alphabet_property(self):
        seq = pybio.Sequence(seq_string)
        assert sorted(seq.alphabet) == sorted(set(seq_string))

    def test_string_property(self):
        seq = pybio.Sequence(seq_string)
        assert seq.string == seq_string

    def test_to_string(self):
        seq = pybio.Sequence(seq_string)
        assert seq.to_string() == seq_string

    def test_to_list(self):
        seq = pybio.Sequence(seq_string)
        assert seq.to_list() == list(seq_string)

    def test_slice(self):
        seq = pybio.Sequence(seq_string)
        assert seq[3:6].string == seq_string[3:6]

    def test_len(self):
        seq = pybio.Sequence(seq_string)
        assert len(seq) == len(seq_string)

    def test_iter(self):
        seq = pybio.Sequence(seq_string)
        assert [x for x in seq] == [x for x in seq_string]

class TestDnaSequence:
    
    def test_init(self):
        pybio.DnaSequence(dna_string)

    def test_bad_character(self):
        with pytest.raises(ValueError):
            pybio.DnaSequence(dna_bad)

    def test_complement(self):
         assert pybio.DnaSequence(dna_string).complement() == pybio.DnaSequence(dna_complement)

    def test_reverse_complement(self):
        assert pybio.DnaSequence(dna_string).reverse_complement() == pybio.DnaSequence(dna_reverse_complement)

    def test_transcribe(self):
        assert pybio.DnaSequence(dna_string).transcribe() == pybio.RnaSequence(dna_transcribe)

    def test_translate(self):
        assert pybio.DnaSequence(dna_string).translate() == pybio.ProteinSequence(dna_translate)

class TestRnaSequence:
    def test_init(self):
        pybio.RnaSequence(rna_string)
    
    def test_bad_character(self):
        with pytest.raises(ValueError):
            pybio.RnaSequence(rna_bad)

    def test_complement(self):
         assert pybio.RnaSequence(rna_string).complement() == pybio.RnaSequence(rna_complement)

    def test_reverse_complement(self):
        assert pybio.RnaSequence(rna_string).reverse_complement() == pybio.RnaSequence(rna_reverse_complement)

    def test_reverse_transcribe(self):
        assert pybio.RnaSequence(rna_string).reverse_transcribe() == pybio.DnaSequence(rna_reverse_transcribe)

    def test_translate(self):
        assert pybio.RnaSequence(rna_string).translate() == pybio.ProteinSequence(rna_translate)

class TestProteinSequence:

    def test_init(self):
        pybio.ProteinSequence(protein_string)

    def test_bad_character(self):
        with pytest.raises(ValueError):
            pybio.ProteinSequence(protein_bad)

