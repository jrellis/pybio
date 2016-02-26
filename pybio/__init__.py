"""
PyBio - a bioinformatics data analysis library for Python

See https://github.com/jrellis/pybio for more documentation, or see the
docstrings of classes, modules, and functions in the pybio namespace:

Sequence
DnaSequence
RnaSequence
ProteinSequence
entrez
alignment
parse
"""

from sequence import *
import alignment, parse, entrez
import logging

logging.getLogger('pybio').addHandler(logging.NullHandler())

def _run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

if _run_from_ipython():
    logging.basicConfig(level=logging.INFO)
