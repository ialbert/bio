"""
Wrappers to interact with NCBI datasets API.
"""

import requests

from biorun import utils


DATASETS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"


def check_accession(value):
    """Validate and normalize an NCBI assembly accession."""
    acc = value.upper()
    url = f"{DATASETS_API}/{acc}/check"

    try:
        response = requests.get(url, headers={"Accept": "application/json"}, timeout=60)
        response.raise_for_status()
        data = response.json()
    except (requests.RequestException, ValueError) as exc:
        utils.error(f"Error checking {acc}: {exc}")

    valid = {item.upper() for item in data.get("valid_assemblies", [])}
    if acc not in valid:
        utils.error(f"invalid accession number: {acc}")

    return acc

