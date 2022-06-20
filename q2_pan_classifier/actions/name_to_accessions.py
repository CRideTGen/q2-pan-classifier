
from qiime2 import Metadata
from qiime2.plugin import Str

from q2_pan_classifier.format_types import NCBIAccFile


def name_to_accessions(taxon_name: str) -> list:

    with open("/scratch/cridenour/Projects/Qiime2Plugins/q2-pan-classifier/q2_pan_classifier/data/test_data.txt", 'r') as rr:
        tt = rr.readlines()
    return tt
