from sequence import *
import alignment, parse
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
