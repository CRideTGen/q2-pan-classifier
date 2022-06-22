#   Copyright 2021 Chase Ridenour

#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import time
import pkg_resources

import skbio
from Bio import Entrez
import pandas as pd
import qiime2.plugin.model as model
from q2_types.feature_data import DNAFASTAFormat
from qiime2.plugin import SemanticType, TextFileFormat
from contextlib import contextmanager

DNAFastaNCBI = SemanticType('DNAFastaNCBI')
NCBIAccFile = SemanticType("NCBIAccFile")


class NCBIAccFileFormat(model.TextFileFormat):
    """

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _validate_(self, level):
        pass

    @contextmanager
    def _open_acc_file_(self):
        acc_file = self.open()
        try:
            yield acc_file
        finally:
            acc_file.close()

    def get_accession_numbers(self) -> list:
        with self._open_acc_file_() as acc_file:
            return acc_file.readlines()


class EntrezFetch:
    def __init__(self, email, db, ids, rettype, retmode):
        Entrez.email = email
        self.handle = Entrez.efetch(db=db,
                                    id=ids,
                                    rettype=rettype,
                                    retmode=retmode)

    def __enter__(self):
        records = Entrez.read(self.handle)
        return records

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.handle.close()


def _accession_number_split_(accession_numbers: list, sub_list_size: int = 100) -> list:
    out_list = list()
    ac_len = len(accession_numbers)

    total_count = 0
    sublist = list()

    while total_count < ac_len:

        index = 0
        while index < sub_list_size and total_count < ac_len:
            sublist.append(accession_numbers[total_count])
            index += 1
            total_count += 1

        out_list.append(sublist.copy())
        sublist.clear()

    return out_list


def _find_tax_(tax_df_input: pd.DataFrame, taxon: str) -> int:
    genus = tax_df_input["Genus"]
    if genus.str.contains(taxon).any():
        index = genus.index[genus.str.contains(taxon) == True][0]
    else:
        index = -1

    return index


def _get_taxonomy_(ac_numbers_subset: list, tax_df: pd.DataFrame, email: str) -> list:
    names = ["Unranked", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
    taxonomy = []

    with EntrezFetch(email=email, db="nuccore", ids=ac_numbers_subset, rettype="gb", retmode="xml") as records:

        for rec in records:
            tax_tmp = rec['GBSeq_taxonomy']
            tax_split = tax_tmp.split(';')
            scientific_name = rec['GBSeq_organism']

            for tax_level in reversed(tax_split):
                index = _find_tax_(tax_df, tax_level.strip())

                if not index == -1:
                    tax_raw = tax_df.loc[index, names]
                    tax_list = tax_raw.to_list()

                    tax_str = ";".join(tax_list)

                    taxonomy.append(tax_str + ';' + scientific_name)
                    break

    return taxonomy


class DNAFastaNCBIFormat(DNAFASTAFormat):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.PIPE = '|'
        self.accession_numbers = []
        self.names = []
        self.taxonomy = []
        self.email = "clr96@nau.edu"

    def get_accession_numbers(self) -> None:
        samps = skbio.read(str(self), format="fasta")

        for samp in samps:

            name = samp.metadata['id']

            self.names.append(name.strip())

            if "gb:" in name and self.PIPE in name:
                a_n_tmp = name.split(":")[1].split(self.PIPE)[0]
                self.accession_numbers.append(a_n_tmp)
            elif self.PIPE in name:
                a_n_tmp = name.split(self.PIPE)[1]
                self.accession_numbers.append(a_n_tmp)
            else:
                a_n_tmp = name.split("-")[0]
                self.accession_numbers.append(a_n_tmp)

    def get_taxonomy(self):

        if not self.accession_numbers:
            raise ValueError("Missing accession numbers")

        taxonomy_reference_csv = pkg_resources.resource_filename('q2_pan_classifier', 'data/ViralTaxonomy.csv')

        tax_df = pd.read_csv(taxonomy_reference_csv)

        accession_number_subsets = _accession_number_split_(self.accession_numbers)

        for subset in accession_number_subsets:
            taxonomy_out = _get_taxonomy_(subset, tax_df, self.email)
            self.taxonomy = taxonomy_out
            time.sleep(1)

    def _validate_(self, level):
        super()._validate_(level=level)

        # fasta file check
        # get it from DNAFASTA
        # Check if sequence names have pipes
        # Check if all accession numbers are valid
        # if not raise error
        # use model.ValidationError
        # TODO: make validatoin function to check if ncbi taxonomic names are there


DNAFastaNCBIDirFormat = model.SingleFileDirectoryFormat('DNAFastaNCBIDirFormat', 'taxonomy.tsv', DNAFastaNCBIFormat)
NCBIAccFileDirectoryFormat = model.SingleFileDirectoryFormat('NCBIAccFileDirectoryFormat', 'acc.txt', NCBIAccFileFormat)
