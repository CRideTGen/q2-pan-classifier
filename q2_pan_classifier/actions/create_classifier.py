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
import qiime2

def generate_taxonomy(ref_seqs: pd.Series) -> list:

    seq_names = [name.metadata['id'] for name in ref_seqs]

    return seq_names

def create_classifier(ctx,
                      ref_seqs_file,
                      ref_tax=None,
                      ref_tax_file=None,
                      f_primer=None,
                      r_primer=None,
                      min_len=None,
                      max_len=None):

    #importing external plugins to be used later
    extract_refs = ctx.get_action('feature_classifier', 'extract_reads')
    train_classifier = ctx.get_action('feature_classifier', 'fit_classifier_naive_bayes')

    results = []

    #importing reference sequences and taxonomy
    ref_seqs = qiime2.Artifact.import_data(type='FeatureData[Sequence]',
                                           view=ref_seqs_file,
                                           view_type=None)

    if ref_tax and ref_tax_file:
        raise ValueError("Please only provide ref_tax OR ref_tax_file not both")

    if ref_tax:
        ref_tax_out = ref_tax
    elif ref_tax_file:
        ref_tax_out = qiime2.Artifact.import_data(type='FeatureData[Taxonomy]',
                                                  view=ref_tax_file,
                                                  view_type=None)
    else:
        ref_tax_out = qiime2.Artifact.import_data(type='FeatureData[Taxonomy]',
                                                  view=ref_seqs_file,
                                                  view_type='DNAFastaNCBIFormat')

    # using imported plugins to extract reference and train classifier
    if f_primer and r_primer:
        if min_len and max_len:
            trimmed_refs = extract_refs(sequences=ref_seqs,
                                        f_primer=f_primer,
                                        r_primer=r_primer,
                                        min_length=min_len,
                                        max_length=max_len)
        else:
            trimmed_refs = extract_refs(sequences=ref_seqs,
                                        f_primer=f_primer,
                                        r_primer=r_primer)

        trained_class = train_classifier(reference_reads=trimmed_refs.reads,
                                         reference_taxonomy=ref_tax_out)
        results += trimmed_refs
    else:
        trained_class = train_classifier(reference_reads=ref_seqs,
                                         reference_taxonomy=ref_tax_out)
        results += [ref_seqs]

    results += [ref_tax_out]
    results += trained_class

    return tuple(results)