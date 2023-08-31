from Bio import Entrez


def name_to_accessions(taxon_name: str) -> list:
    Entrez.email = "clr96@nau.edu"

    handle = Entrez.esearch(
        db="nuccore",
        retMax=100,
        term=f"{taxon_name}[ORGANISM] AND 5000:1000000000[SLEN]",
        idtype="acc",
    )
    record = Entrez.read(handle)

    return [acc for acc in record["IdList"]]
