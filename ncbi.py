import requests


def fetch_sequence(accession):

    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        f"efetch.fcgi?db=nuccore&id={accession}"
        f"&rettype=fasta&retmode=text"
    )

    response = requests.get(url)

    if response.status_code == 200:
        return response.text

    return None
