from sequence import *
import alignment
import logging

logging.getLogger('pybio').addHandler(logging.NullHandler())

def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

if run_from_ipython():
    logging.basicConfig(level=logging.INFO)
