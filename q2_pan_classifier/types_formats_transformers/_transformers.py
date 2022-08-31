#   Copyright 2021 Chase Ridenour
#
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
import pandas as pd
from q2_pan_classifier.types_formats_transformers._format import DNAFastaNCBIFormat, NCBIAccFileFormat
from q2_pan_classifier.plugin_setup import plugin

from q2_types.feature_data import TSVTaxonomyFormat, DNAFASTAFormat



@plugin.register_transformer
def _2(ref_seqs: list) -> TSVTaxonomyFormat:

    tax_out = TSVTaxonomyFormat()

    with open(tax_out.path, 'w') as ff:
        ff.write('\t'.join(tax_out.HEADER) + '\n')
        for name in ref_seqs:
            ff.write('\t'.join([name, 'virus']) + '\n')

    return tax_out


@plugin.register_transformer
def _3(ref_seqs: DNAFASTAFormat) -> TSVTaxonomyFormat:

    seq_names = [name.metadata['id'] for name in ref_seqs.view(pd.Series)]

    tax_out = TSVTaxonomyFormat()

    with open(tax_out.path, 'w') as ff:
        ff.write('\t'.join(tax_out.HEADER) + '\n')
        for name in seq_names:
            ff.write('\t'.join([name, 'virus']) + '\n')

    return tax_out


@plugin.register_transformer
def _4(ff: DNAFastaNCBIFormat) -> list:
    return ff.accession_numbers

@plugin.register_transformer
def _5(ref_seqs: DNAFastaNCBIFormat) -> TSVTaxonomyFormat:
    ref_seqs.get_accession_numbers()
    ref_seqs.get_taxonomy()
    seq_names = ref_seqs.names
    seq_taxonomy = ref_seqs.taxonomy

    tax_out = TSVTaxonomyFormat()

    with open(tax_out.path, 'w') as ff:
        ff.write('\t'.join(tax_out.HEADER) + '\n')
        for name, tax in zip(seq_names, seq_taxonomy):
            ff.write('\t'.join([name, tax]) + '\n')

    return tax_out

@plugin.register_transformer
def _6(input_list: list) -> NCBIAccFileFormat:

    out_path = NCBIAccFileFormat()

    with open(out_path.path, 'w') as ff:
        for i in input_list:
            ff.write(f"{i.strip()}\n")

    return out_path

@plugin.register_transformer
def _7(input_acc: NCBIAccFileFormat) -> list:
    return input_acc.get_accession_numbers()


#TODO Need to write tests for the transformer
#test plugin base from qiime2.plugin
#tests embedded in module