import requests, sys, json

URL = 'https://rest.ensembl.org'

ENA_SEARCH_URL = f"{URL}/lookup/id"

def search(acc, url=ENA_SEARCH_URL):
    """
    Connect to the lookup service.
    """
    conn = f"{url}/{acc}"
    r = requests.get(conn, headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    data = r.json()
    try:
        text = json.dumps(data, indent=4)
        print (text)
    except Exception as exc:
        print(f"# JSON decoding error: {exc}")
        sys.exit(1)

if __name__ == '__main__':

    elems = [
           'ENSG00000157764',
            #'ENST00000361390', 'ENSP00000354686',
    ]

    for acc in elems:
        search(acc)
