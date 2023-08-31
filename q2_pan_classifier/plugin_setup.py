#   Copyright 2021 Evan Bolyen
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
import importlib

from q2_dada2._stats import DADA2Stats
from q2_feature_classifier.classifier import TaxonomicClassifier, Taxonomy
from q2_types.feature_data import FeatureData, Sequence
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.per_sample_sequences import PairedEndSequencesWithQuality, SequencesWithQuality
from q2_types.sample_data import SampleData
from qiime2.plugin import Plugin, Visualization
from qiime2.plugin import Str, Int, Range

import q2_pan_classifier.actions as actions
from q2_pan_classifier.format_types import (DNAFastaNCBI, DNAFastaNCBIFormat, DNAFastaNCBIDirFormat, NCBIAccFile,
                                            NCBIAccFileFormat, NCBIAccFileDirectoryFormat)

# This is the plugin object. It is what the framework will load and what an
# interface will interact with. Basically every registration we perform will
# involve this object in some way.
plugin = Plugin(name="pan-classifier",
                version="0.0.1.dev",
                package="q2_pan_classifier",
                website="https://github.com/ebolyen/q2-reveal")

plugin.register_semantic_types(DNAFastaNCBI, NCBIAccFile)
plugin.register_semantic_type_to_format(DNAFastaNCBI, DNAFastaNCBIDirFormat)
plugin.register_semantic_type_to_format(NCBIAccFile, NCBIAccFileDirectoryFormat)
plugin.register_formats(DNAFastaNCBIFormat, DNAFastaNCBIDirFormat)
plugin.register_formats(NCBIAccFileFormat, NCBIAccFileDirectoryFormat)

importlib.import_module("q2_pan_classifier.transformers")

plugin.methods.register_function(

    function=actions.name_to_accessions,
    inputs={},
    outputs=[('acc_file', NCBIAccFile)],
    parameters={'taxon_name': Str},
    input_descriptions=None,
    parameter_descriptions=None,
    name='Get accession numbers',
    description="Gets Accession Numbers"
)

# plugin.pipelines.register_function(
#     function=actions.get_taxonomies_by_taxon_level,
#     inputs={'accession_file': NCBIAccFile},
#     outputs=[('ref_seqs', FeatureData[Sequence]),
#              ('tax_ref', FeatureData[Taxonomy]),
#              ],
#     parameters={},
#     input_descriptions=None,
#     parameter_descriptions={
#     },
#     output_descriptions={'ref_seqs': 'Path where reference sequence Artifact will be written',
#                          'trained_classifier': 'Path where the trained classifier Artifact will be written'
#                          },
#     name='Get_Taxonomies_by_name',
#     description="test"
# )

plugin.methods.register_function(
    function=actions.generate_taxonomy,
    inputs={
        'ref_seqs': FeatureData[Sequence]
    },
    outputs=[('taxonomy', FeatureData[Taxonomy])],
    parameters={},
    input_descriptions=None,
    parameter_descriptions=None,
    name='Generate Taxonomy',
    description="Creates a taxonomy from reference fasta file"
)

plugin.pipelines.register_function(
    function=actions.create_classifier,
    inputs={'ref_tax': FeatureData[Taxonomy]},
    outputs=[('ref_seqs', FeatureData[Sequence]),
             ('tax_ref', FeatureData[Taxonomy]),
             ('trained_classifier', TaxonomicClassifier)],
    parameters={
        'ref_seqs_file': Str,
        'ref_tax_file': Str,
        'f_primer': Str,
        'r_primer': Str,
        'min_len': Int % Range(0, None),
        'max_len': Int % Range(0, None)
    },
    input_descriptions=None,
    parameter_descriptions={
        'ref_seqs_file': "File path to the reference sequence fasta file. ",
        'ref_tax_file': "File path to the reference taxonomy file. Must be a TSV file. ",
        'f_primer': "Forward primer to trim the reference sequences to be only the amplified region",
        'r_primer': "Reverse primer to trim the reference sequences to be only the amplified region",
        'min_len': "Minimum length of the trimmed sequences. Any sequences shorter than this will be removed from the analysis",
        'max_len': "Maximum length of the trimmed sequences. Any sequences longer than this will be removed from the analysis"
    },
    output_descriptions={'ref_seqs': 'Path where reference sequence Artifact will be written',
                         'trained_classifier': 'Path where the trained classifier Artifact will be written'
                         },
    name='Create Classifier',
    description="test"
)

plugin.pipelines.register_function(
    function=actions.prep_sequence_reads_single,
    inputs={'sequences': SampleData[SequencesWithQuality]},
    outputs=[
        ('trimmed_reads', SampleData[SequencesWithQuality]),
        ('table_viz', Visualization)
    ],
    parameters={
        'sequences_directory': Str,
        'metadata_template_dir': Str,
        'primer_f': Str,
        'primer_r': Str
    },
    input_descriptions=None,
    parameter_descriptions={
        'sequences_directory': 'Path to sequences directory',
        'primer_f': 'Sequence of an adapter ligated to the 3\' end. The '
                    'adapter and any subsequent bases are trimmed. If a `$` '
                    'is appended, the adapter is only found if it is at the '
                    'end of the read. Search in forward read. If your '
                    'sequence of interest is "framed" by a 5\' and a 3\' '
                    'adapter, use this parameter to define a "linked" primer '
                    '- see https://cutadapt.readthedocs.io for complete '
                    'details.',
        'primer_r': 'Sequence of an adapter ligated to the 3\' end. The '
                    'adapter and any subsequent bases are trimmed. If a `$` '
                    'is appended, the adapter is only found if it is at the '
                    'end of the read. Search in reverse read. If your '
                    'sequence of interest is "framed" by a 5\' and a 3\' '
                    'adapter, use this parameter to define a "linked" primer '
                    '- see https://cutadapt.readthedocs.io for complete '
                    'details.'
    },
    output_descriptions={
        'table_viz': 'Visualization of the demultiplexed sequences to assess read quality.'
                     'The Dada2 truncation inputs will be determined using this visualization.'
    },
    name='Read Quality Visualization',
    description="test"
)

plugin.pipelines.register_function(
    function=actions.prep_sequence_reads_paired,
    inputs={'sequences': SampleData[PairedEndSequencesWithQuality]},
    outputs=[
        ('trimmed_reads', SampleData[PairedEndSequencesWithQuality]),
        ('table_viz', Visualization)
    ],
    parameters={
        'sequences_directory': Str,
        'metadata_template_dir': Str,
        'primer_f': Str,
        'primer_r': Str
    },
    input_descriptions=None,
    parameter_descriptions={
        'sequences_directory': 'Path to sequences directory',
        'primer_f': 'Sequence of an adapter ligated to the 3\' end. The '
                    'adapter and any subsequent bases are trimmed. If a `$` '
                    'is appended, the adapter is only found if it is at the '
                    'end of the read. Search in forward read. If your '
                    'sequence of interest is "framed" by a 5\' and a 3\' '
                    'adapter, use this parameter to define a "linked" primer '
                    '- see https://cutadapt.readthedocs.io for complete '
                    'details.',
        'primer_r': 'Sequence of an adapter ligated to the 3\' end. The '
                    'adapter and any subsequent bases are trimmed. If a `$` '
                    'is appended, the adapter is only found if it is at the '
                    'end of the read. Search in reverse read. If your '
                    'sequence of interest is "framed" by a 5\' and a 3\' '
                    'adapter, use this parameter to define a "linked" primer '
                    '- see https://cutadapt.readthedocs.io for complete '
                    'details.'
    },
    output_descriptions={
        'table_viz': 'Visualization of the demultiplexed sequences to assess read quality.'
                     'The Dada2 truncation inputs will be determined using this visualization.'
    },
    name='Read Quality Visualization',
    description="test"
)

plugin.pipelines.register_function(
    function=actions.classify_reads,
    inputs={'samp_reads': SampleData[PairedEndSequencesWithQuality],
            'trained_classifier': TaxonomicClassifier,
            'dada2_table': FeatureTable[Frequency],
            'dada2_rep_seqs': FeatureData[Sequence],
            'dada2_stats': SampleData[DADA2Stats]
            },
    outputs=[('classified', FeatureData[Taxonomy]),
             ('barplot_taxonomy', Visualization),
             ('merged_table', Visualization),
             ('dada2_table_out', FeatureTable[Frequency]),
             ('dada2_rep_seqs_out', FeatureData[Sequence]),
             ('dada2_stats_out', SampleData[DADA2Stats])
             ],
    parameters={
        'trunc_len_f': Int % Range(0, None),
        'trunc_len_r': Int % Range(0, None),
    },
    input_descriptions={
        'samp_reads': 'Path to sample reads Artifact (.qza file)',
        'trained_classifier': 'Path to trained classifier Artifact (.qza file)',
    },
    parameter_descriptions={
        'trunc_len_f': ('Position at which forward read sequences should be '
                        'truncated due to decrease in quality. This truncates '
                        'the 3\' end of the of the input sequences, which '
                        'will be the bases that were sequenced in the last '
                        'cycles. Reads that are shorter than this value '
                        'will be discarded. After this parameter is applied '
                        'there must still be at least a 12 nucleotide overlap '
                        'between the forward and reverse reads. If 0 is '
                        'provided, no truncation or length filtering will be '
                        'performed'),
        'trunc_len_r': ('Position at which reverse read sequences should be '
                        'truncated due to decrease in quality. This truncates '
                        'the 3\' end of the of the input sequences, which '
                        'will be the bases that were sequenced in the last '
                        'cycles. Reads that are shorter than this value '
                        'will be discarded. After this parameter is applied '
                        'there must still be at least a 12 nucleotide overlap '
                        'between the forward and reverse reads. If 0 is '
                        'provided, no truncation or length filtering will be '
                        'performed')
    },
    output_descriptions={
        'dada2_table_out': 'The resulting feature table.',
        'dada2_rep_seqs_out': 'The resulting feature sequences. Each '
                              'feature in the feature table will be '
                              'represented by exactly one sequence.',
        'dada2_stats_out': 'Stats on Dada2 clustering and filtering',
        'classified': 'Resulting Taxonomic Artifact from the classification'
    },
    name='Classify Reads',
    description="Using a trained classifier to classify unknown reads"
)

plugin.pipelines.register_function(
    function=actions.classify_reads_single,
    inputs={'samp_reads': SampleData[SequencesWithQuality],
            'trained_classifier': TaxonomicClassifier,
            'dada2_table': FeatureTable[Frequency],
            'dada2_rep_seqs': FeatureData[Sequence],
            'dada2_stats': SampleData[DADA2Stats]
            },
    outputs=[('classified', FeatureData[Taxonomy]),
             ('barplot_taxonomy', Visualization),
             ('merged_table', Visualization),
             ('dada2_table_out', FeatureTable[Frequency]),
             ('dada2_rep_seqs_out', FeatureData[Sequence]),
             ('dada2_stats_out', SampleData[DADA2Stats])
             ],
    parameters={
        'trunc_len': Int % Range(0, None)
    },
    input_descriptions={
        'samp_reads': 'Path to sample reads Artifact (.qza file)',
        'trained_classifier': 'Path to trained classifier Artifact (.qza file)',
    },
    parameter_descriptions={
        'trunc_len': ('Position at which forward read sequences should be '
                      'truncated due to decrease in quality. This truncates '
                      'the 3\' end of the of the input sequences, which '
                      'will be the bases that were sequenced in the last '
                      'cycles. Reads that are shorter than this value '
                      'will be discarded. After this parameter is applied '
                      'there must still be at least a 12 nucleotide overlap '
                      'between the forward and reverse reads. If 0 is '
                      'provided, no truncation or length filtering will be '
                      'performed')
    },
    output_descriptions={
        'dada2_table_out': 'The resulting feature table.',
        'dada2_rep_seqs_out': 'The resulting feature sequences. Each '
                              'feature in the feature table will be '
                              'represented by exactly one sequence.',
        'dada2_stats_out': 'Stats on Dada2 clustering and filtering',
        'classified': 'Resulting Taxonomic Artifact from the classification'
    },
    name='Classify Reads',
    description="Using a trained classifier to classify unknown reads"
)

plugin.visualizers.register_function(
    function=actions.visualization_final,
    inputs={
        # 'feat_table': FeatureTable[Frequency]
    },
    parameters=None,
    input_descriptions={},
    parameter_descriptions={},
    name='Final Visualization',
    description=("Making things look like things")
)
