""" Functions for accessing the NCBI Entrez databases """

import logging
import requests
import time

__all__ = ['set_email', 'set_tool', 'einfo', 'esearch', 'epost', 'esummary', 'efetch', 'elink', 'egquery', 'espell', 'ecitmatch']

logger = logging.getLogger(__name__)

default_email=None
default_tool='pybio'
last_access_time=0
ncbi_delay=0.3334
eutils_url_template='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{eutil_tool}.fcgi'

def set_email(email):
    """Set the Entrez email to this email address."""
    global default_email
    default_email = email

def set_tool(tool):
    """Set the Entrez tool to this name."""
    global default_tool 
    default_tool = tool

def set_ncbi_delay(delay):
    """
    Set the NCBI delay time.

    This is set to ensure that the maximum number of queries
    per second is not exceeded.
    """
    global ncbi_delay 
    ncbi_delay = delay

def _add_email_and_tool(kwargs):
    if not 'email' in kwargs:
        if not default_email:
            logger.warning('No email set for Entrez. Please use pybio.entrez.set_email() to set an email address.')
        else:
            kwargs['email'] = default_email
    if not 'tool' in kwargs:
        kwargs['tool'] = default_tool

def _ncbi_delay():
    time_since_last_call = time.time() - last_access_time
    if time_since_last_call < ncbi_delay:
        time.sleep(time_since_last_call)

def _entrez_post(eutil_tool, **kwargs):
    _ncbi_delay()
    _add_email_and_tool(kwargs)
    response = requests.post(eutils_url_template.format(**locals()), params=kwargs)
    last_access = time.time()
    return response

def _entrez_get(eutil_tool, **kwargs):
    _ncbi_delay()
    _add_email_and_tool(kwargs)
    response = requests.get(eutils_url_template.format(**locals()), params=kwargs)
    last_access = time.time()
    return response

def einfo(**kwargs):
    """
    Direct access to the Entrez EInfo utility
    
    Provides information on Entrez databases.

    Parameters
    ----------
    db
        the Entrez database. If not specified, a list of all valid
        Entrez databases is returned.
    version
        the version of the EInfo XML returned
    retmode
        the format of the output: 'xml' or 'json'
    """
    return _entrez_get('einfo', **kwargs)

def esearch(**kwargs):
    """
    Direct access to the Entrez ESearch utility
    
    Returns UIDs matching a query.
    """
    return _entrez_post('esearch', **kwargs)

def epost(**kwargs):
    """
    Direct access to the Entrez EPost utility
    
    Stores a set of UIDs from a given database on the History Server.
    Returns the query key and web environment for the dataset.
    """
    return _entrez_post('epost', **kwargs)

def esummary(**kwargs):
    """
    Direct access to the Entrez ESummary utility
    
    Returns document summaries matching a set of UIDs.
    """
    return _entrez_post('esummary', **kwargs)

def efetch(**kwargs):
    """
    Direct access to the Entrez EFetch utility
    
    Returns data records matching a set of UIDs.
    """
    return _entrez_post('efetch', **kwargs)

def elink(**kwargs):
    """
    Direct access to the Entrez ELink utility
    
    Returns UIDs linked to an input set of UIDs.
    """
    return _entrez_post('efetch', **kwargs)

def egquery(**kwargs):
    """
    Direct access to the Entrez EGQuery utility
    
    Returns the number of records matching a query.
    """
    return _entrez_post('egquery', **kwargs)

def espell(**kwargs):
    """
    Direct access to the Entrez ESpell utility
    
    Returns spelling suggestions for a query in a database.
    """
    return _entrez_post('espell', **kwargs)

def ecitmatch(**kwargs):
    """
    Direct access to the Entrez ECitMatch utility
    
    Returns PubMed IDs.
    """
    return _entrez_post('ecitmatch', **kwargs)


