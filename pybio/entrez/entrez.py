""" Functions for accessing the NCBI Entrez databases """

import logging
import requests
import time
from threading import Thread
import Queue
from xml.etree import cElementTree
from pybio import ProteinSequence, DnaSequence

__all__ = ['set_email', 'set_tool', 'einfo', 'esearch', 'epost', 'esummary', 'efetch', 'elink', 'egquery', 'espell', 'ecitmatch', 'esearch_id', 'esearch_esummary', 'esearch_efetch', 'search_proteins', 'search_nucleotide_sequences', 'get_entry', 'get_nucleotide_sequence', 'get_protein']

logger = logging.getLogger(__name__)

default_email=None
default_tool='pybio'
last_access_time=0
ncbi_delay=0.3334
eutils_url_template='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{eutil_tool}.fcgi'
ecitmatch_url_template='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/{eutil_tool}.cgi'

# Helper functions for Entrez

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

def _entrez_post(eutil_tool, template=eutils_url_template, **kwargs):
    _ncbi_delay()
    _add_email_and_tool(kwargs)
    response = requests.post(template.format(**locals()), params=kwargs)
    last_access = time.time()
    return response

def _entrez_get(eutil_tool, **kwargs):
    _ncbi_delay()
    _add_email_and_tool(kwargs)
    response = requests.get(eutils_url_template.format(**locals()), params=kwargs)
    last_access = time.time()
    return response

def _call_check(call, **kwargs):
    result = call(**kwargs)
    if not result.ok:
        raise RuntimeError('{call} failed with status code {status_code}: {reason}'.format(call=call, status_code=search_result.status_code, reason=search_result.reason))
    return result

# Direct wrappers over Entrez tools

def einfo(**kwargs):
    """
    Direct access to the Entrez EInfo utility
    
    Provides information on Entrez databases.
    """
    return _entrez_get('einfo', **kwargs)

def esearch(db, term, **kwargs):
    """
    Direct access to the Entrez ESearch utility
    
    Returns UIDs matching a query.

    Parameters
    ----------
    db
        Entrez database.
    term
        Entrez text query.
    """
    return _entrez_post('esearch', db=db, term=term, **kwargs)

def epost(db, id, **kwargs):
    """
    Direct access to the Entrez EPost utility
    
    Stores a set of UIDs from a given database on the History Server.
    Returns the query key and web environment for the dataset.

    Parameters
    ----------
    db
        Entrez database.
    id
        Comma separated Entrez UID list. 
    """
    return _entrez_post('epost', db=db, id=id, **kwargs)

def esummary(db, **kwargs):
    """
    Direct access to the Entrez ESummary utility
    
    Returns document summaries matching a set of UIDs.

    Parameters
    ----------
    db
        Entrez database.
    """
    return _entrez_post('esummary', db=db, **kwargs)

def efetch(db, **kwargs):
    """
    Direct access to the Entrez EFetch utility
    
    Returns data records matching a set of UIDs.

    Parameters
    ----------
    db
        Entrez database.
    """
    return _entrez_post('efetch', db=db, **kwargs)

def elink(db, dbfrom, cmd, **kwargs):
    """
    Direct access to the Entrez ELink utility
    
    Returns UIDs linked to an input set of UIDs.

    Parameters
    ----------
    db
        Entrez database from which to retrieve UIDs.
    dbfrom
        Entrez database with the input UIDs
    cmd
        ELink command mode.
    """
    return _entrez_post('efetch', db=db, dbfrom=dbfrom, cmd=cmd, **kwargs)

def egquery(term, **kwargs):
    """
    Direct access to the Entrez EGQuery utility
    
    Returns the number of records matching a query.

    Parameters
    ----------
    term
        Entrez text query.
    """
    return _entrez_post('egquery', term=term, **kwargs)

def espell(db, term, **kwargs):
    """
    Direct access to the Entrez ESpell utility
    
    Returns spelling suggestions for a query in a database.

    Parameters
    ----------
    db
        Entrez database.
    term
        Entrez text query.
    """
    return _entrez_post('espell', db=db, term=term, **kwargs)

def ecitmatch(bdata, db='pubmed', rettype='xml', **kwargs):
    """
    Direct access to the Entrez ECitMatch utility
    
    Returns PubMed IDs.

    Parameters
    ----------
    bdata
        Citation string in this format: journal_title|year|volume|first_page|author_name|your_key|
    db
        Entrez database. Only 'pubmed' is supported.
    rettype
        Retrieval type. Only 'xml' is supported.
    """
    return _entrez_post('ecitmatch', template=ecitmatch_url_template, db=db, rettype=rettype, bdata=bdata, **kwargs)

# Entrez Pipelines

def esearch_id(**kwargs):
    """
    Make an esearch Entrez query and return a list of ids.
    """
    search_result = _call_check(esearch, **kwargs)
    element_tree = cElementTree.fromstring(search_result.text)
    return [id.text for id in element_tree.find('IdList').findall('Id')]
    
def esearch_history(**kwargs):
    """
    Return the QueryKey and WebEnv matching an Entrez esearch query.
    """
    kwargs['usehistory'] = 'y'
    search_result = _call_check(esearch, **kwargs)
    element_tree = cElementTree.fromstring(search_result.text)
    return element_tree.find('QueryKey').text, element_tree.find('WebEnv').text

def esearch_esummary(db, **kwargs):
    """
    Return summaries of data records matching an Entrez esearch query.
    """
    query_key, web_env = esearch_history(db=db, **kwargs)
    kwargs['query_key'] = query_key
    kwargs['WebEnv'] = web_env
    return _call_check(esummary, db=db, **kwargs)

def esearch_efetch(db, **kwargs):
    """
    Return summaries of data records matching an Entrez esearch query.
    """
    query_key, web_env = esearch_history(db=db, **kwargs)
    kwargs['query_key'] = query_key
    kwargs['WebEnv'] = web_env
    return _call_check(efetch, db=db, **kwargs)

def get_entry(db, accession, return_type='fasta'):
    """
    Get an Entrez entry by accession number.

    Parameters
    ----------
    db
        Entrez database
    accession
        accession number
    return_type
        return type from Entrez, corresponding to rettype in Entrez calls
    """
    term = '{accession}[accession]'.format(**locals())
    return esearch_efetch(db=db, rettype=return_type, term=term)

def get_nucleotide_sequence(accession):
    """
    Get an Entrez nucleotide entry by accession number.

    Parameters
    ----------
    accession
        accession number
    """
    result = get_entry('nucleotide', accession).text
    if '<ERROR>' in result:
        error = cElementTree.fromstring(result).find('ERROR').text
        raise ValueError('Error accessing {accession}: {error}'.format(**locals()))
    return DnaSequence.from_fasta(result)

def get_protein(accession):
    """
    Get an Entrez protein entry by accession number.

    Parameters
    ----------
    accession
        accession number
    """
    result = get_entry('protein', accession).text
    if '<ERROR>' in result:
        error = cElementTree.fromstring(result).find('ERROR').text
        raise ValueError('Error accessing {accession}: {error}'.format(**locals()))
    return ProteinSequence.from_fasta(result)

def _make_term(**kwargs):
    terms = []
    for key,value in kwargs.iteritems():
        field = key.replace('_', ' ')
        terms.append('{value}[{field}]'.format(**locals()))
    return ' '.join(terms)

def threaded_entrez_function_iter(entrez_function, db, rettype='fasta', max_count=None, max_per_call=100, queue_size=3, **kwargs):
    """
    Returns an iterator over results from an entrez function taking a database and a term (like esearch).
    """
    term = _make_term(**kwargs)
    queue = Queue.Queue(queue_size)
    
    if max_count:
        def not_completed(retstart):
            return max_count - retstart > 0
        def next_retmax(retstart):
            return min(max_per_call, max_count - retstart)
    else:
        def not_completed(retstart):
            return True
        def next_retmax(retstart):
            return max_per_call

    def do_entrez_query():
        try:
            retstart = 0
            while not_completed(retstart):
                retmax = next_retmax(retstart)
                result = entrez_function(db, rettype=rettype, retmax=retmax, retstart=retstart, term=term)
                if not result.ok or '<ERROR>' in result.text:
                    break
                queue.put(result)
                retstart = retstart + max_per_call
        except Exception as exception:
            logger.exception(exception)
        queue.put(None)
        return
    
    thread = Thread(target=do_entrez_query)
    thread.daemon = True
    thread.start()

    while True: 
        result = queue.get()
        if not result:
            raise StopIteration()
        yield result

def search_proteins(max_count=None, max_per_call=100, queue_size=3, **kwargs):
    """
    Return from query terms an iterator of ProteinSequences from Entrez proteins.

    Parameters
    ----------
    max_count
        maximum number of proteins returned
    max_per_call
        maximum number of entries returned per call to Entrez
    queue_size
        the size of the queue holding asyncronous results from Entrez queries

    Search Parameters
    -----------------
    all_fields(ALL)  
        All terms from all searchable fields
    uid(UID)     
        Unique number assigned to each sequence
    filter(FILT)     
        Limits the records
    text_word(WORD)   
        Free text associated with record
    title(TITL)     
        Words in definition line
    keyword(KYWD)     
        Nonstandardized terms provided by submitter
    author(AUTH)      
        Author(s) of publication
    journal(JOUR)     
        Journal abbreviation of publication
    volume(VOL)       
        Volume number of publication
    issue(ISS)        
        Issue number of publication
    page_number(PAGE)
        Page number(s) of publication
    organism(ORGN)
        Scientific and common names of organism, and all higher levels of taxonomy
    accession(ACCN)   
        Accession number of sequence
    primary_accession(PACC)
        Does not include retired secondary accessions
    gene_name(GENE)
        Name of gene associated with sequence
    protein_name(PROT)
        Name of protein associated with sequence
    ECNO
        EC number for enzyme or CAS registry number
    publication_date(PDAT)
        Date sequence added to GenBank
    modification_date(MDAT)
        Date of last update
    substance_name(SUBS)
        CAS chemical name or MEDLINE Substance Name
    properties(PROP)
        Classification by source qualifiers and molecule type
    SeqID(SQID)
        SeqID String
    bioproject(GPRJ)
        BioProject
    sequence_length(SLEN)
        Length of sequence
    molecular_weight(MLWT)
        Molecular Weight
    feature_key(FKEY)
        Feature annotated on sequence
    primary_organism(PORG)
        Scientific and common names of primary organism, and all higher levels of taxonomy
    assembly(ASSM)
        Assembly
    division(DIV)
        Division
    strain(STRN)
        Strain
    isolate(ISOL)
        Isolate
    cultivar(CULT)
        Cultivar
    breed(BRD)
        Breed
    """
    for entry in threaded_entrez_function_iter(esearch_efetch, 'protein', max_count=max_count, max_per_call=max_per_call, queue_size=queue_size, **kwargs):
        for sequence in ProteinSequence.sequences_from_fasta(entry.text):
            yield sequence

def search_nucleotide_sequences(max_count=None, max_per_call=100, queue_size=3, **kwargs):
    """
    Return from query terms an iterator of DnaSequences from Entrez nucleotides.

    Parameters
    ----------
    max_count
        maximum number of nucleotides returned
    max_per_call
        maximum number of entries returned per call to Entrez
    queue_size
        the size of the queue holding asyncronous results from Entrez queries

    Search Parameters
    -----------------
    all_fields(ALL)  
        All terms from all searchable fields
    uid(UID)     
        Unique number assigned to each sequence
    filter(FILT)     
        Limits the records
    text_word(WORD)   
        Free text associated with record
    title(TITL)     
        Words in definition line
    keyword(KYWD)     
        Nonstandardized terms provided by submitter
    author(AUTH)      
        Author(s) of publication
    journal(JOUR)     
        Journal abbreviation of publication
    volume(VOL)       
        Volume number of publication
    issue(ISS)        
        Issue number of publication
    page_number(PAGE)
        Page number(s) of publication
    organism(ORGN)
        Scientific and common names of organism, and all higher levels of taxonomy
    accession(ACCN)   
        Accession number of sequence
    primary_accession(PACC)
        Does not include retired secondary accessions
    gene_name(GENE)
        Name of gene associated with sequence
    protein_name(PROT)
        Name of protein associated with sequence
    ECNO
        EC number for enzyme or CAS registry number
    publication_date(PDAT)
        Date sequence added to GenBank
    modification_date(MDAT)
        Date of last update
    substance_name(SUBS)
        CAS chemical name or MEDLINE Substance Name
    properties(PROP)
        Classification by source qualifiers and molecule type
    SeqID(SQID)
        SeqID String
    bioproject(GPRJ)
        BioProject
    sequence_length(SLEN)
        Length of sequence
    feature_key(FKEY)
        Feature annotated on sequence
    primary_organism(PORG)
        Scientific and common names of primary organism, and all higher levels of taxonomy
    component_accession
        Component accessions for an assembly
    assembly(ASSM)
        Assembly
    division(DIV)
        Division
    strain(STRN)
        Strain
    isolate(ISOL)
        Isolate
    cultivar(CULT)
        Cultivar
    breed(BRD)
        Breed
    BioSample(BIOS)
        BioSample   
    """
    for entry in threaded_entrez_function_iter(esearch_efetch, 'nucleotide', max_count=max_count, max_per_call=max_per_call, queue_size=queue_size, **kwargs):
        for sequence in DnaSequence.sequences_from_fasta(entry.text):
            yield sequence
    

