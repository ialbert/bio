'''
A simplified reimplementation of Entrez Direct interfaces

http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4

'''
import requests
from biorun import utils
from biorun.libs import xmltodict
import json

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

logger = utils.logger

# The keys that are valid environment keys.
ENV_KEYS = [ "webenv", "use_history", "query_key" ]

def efetch(db, env={}, **kwds):
    """
    Perform an efetch query.
    Returns a dictionary based data structure.
    """
    try:

        params = dict(db=db)

        # Fill in with environment.
        for key, value in env.items():
            if key in ENV_KEYS:
                params[key] = value

        params.update(kwds)
        r = requests.get(EFETCH_URL, params=params)

        logger.info(r.url)

        data = xmltodict.parse(r.text)

        # Roundtrip to get rid of OrderedDicts
        data = json.loads(json.dumps(data))

        return data

    except Exception as exc:
        logger.error(exc)



def esearch(db, **kwds):
    """
    Performs an esearch query.
    Returns a dictionary based data structure.
    """
    try:
        params = dict(db=db)
        params.update(kwds)

        r = requests.get(ESEARCH_URL, params=params)

        logger.info(r.url)

        data = xmltodict.parse(r.text)

        # Roundtrip to get rid of OrderedDicts
        data = json.loads(json.dumps(data))

        # Parse the web environment
        field = data['eSearchResult']

        env = dict(
            webenv=field['WebEnv'],
            count=field['Count'],
            query_key=field['QueryKey'],
            use_history="y"
        )

        return env

    except Exception as exc:
        logger.error(exc)
