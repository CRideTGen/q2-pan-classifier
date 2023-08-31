import qiime2
import qiime2.plugins.dada2.actions as dada2_actions
from pathlib import Path
from q2_feature_classifier._taxonomic_classifier import TaxonomicClassifier
from q2_types.feature_data import FeatureData, Sequence, Taxonomy
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.per_sample_sequences import PairedEndSequencesWithQuality
from q2_types.sample_data import SampleData
from qiime2 import Metadata
from qiime2.plugins.cutadapt.methods import trim_paired
from qiime2.plugins.demux.visualizers import summarize as demux_summarize
from qiime2.plugins.feature_classifier.methods import classify_sklearn
from qiime2.plugins.feature_table.methods import transpose
from qiime2.plugins.feature_table.visualizers import summarize as table_summarize
from qiime2.plugins.feature_table.visualizers import tabulate_seqs
from qiime2.plugins.metadata.visualizers import tabulate
from qiime2.plugins.taxa.visualizers import barplot
from qiime2.plugins.feature_table.methods import merge

# first
def import_viral_sequences(path_to_manifest_file: Path):
    viral_seq_art = qiime2.Artifact.import_data(type='SampleData[PairedEndSequencesWithQuality]',
                                                view=str(path_to_manifest_file),
                                                view_type='PairedEndFastqManifestPhred33V2')
    return viral_seq_art

# second
def trim_sequences(demux_seqs: SampleData[PairedEndSequencesWithQuality], forward_primer: list, reverse_primer: list):
    trimmed_seqs = trim_paired(
        demultiplexed_sequences=demux_seqs,
        front_f=forward_primer,
        front_r=reverse_primer
    )

    return trimmed_seqs

# third
def visualize_sequence_quality(trimmed_demux_seqs: SampleData[PairedEndSequencesWithQuality]):
    visual_out = demux_summarize(
        data=trimmed_demux_seqs
    )

    return visual_out

# fourth
def cluster_sequences(demux_seqs: SampleData[PairedEndSequencesWithQuality], trunc_f: int, trunc_r: int):
    table, rep_seqs, stats = dada2_actions.denoise_paired(
        demultiplexed_seqs=demux_seqs,
        trim_left_f=0,
        trim_left_r=0,
        trunc_len_f=trunc_f,
        trunc_len_r=trunc_r
    )

    return table, rep_seqs, stats


# fifth
def visualize_feature_table(table: FeatureTable[Frequency], rep_seqs: FeatureData[Sequence], metadata: Metadata):
    visualization = table_summarize(
        table=table,
        sample_metadata=metadata
    )

    tabulate = tabulate_seqs(
        data=rep_seqs
    )
    return visualization, tabulate

# sixth
def classify_sequences(classifier: TaxonomicClassifier, rep_seqs: FeatureData[Sequence]):
    taxonomy = classify_sklearn(
        classifier=classifier,
        reads=rep_seqs
    )

    return taxonomy

# seventh
def create_taxa_bar_plot(table: FeatureTable[Frequency], taxonomy: FeatureData[Taxonomy], metadata: Metadata):
    visualization = barplot(
        table=table,
        taxonomy=taxonomy,
        metadata=metadata
    )

    return visualization

# eighth
def transpose_table(table: FeatureTable[Frequency]):
    transposed_table = transpose(
        table=table
    )

    return transposed_table

# ninth
def tabulate_table(taxonomy: FeatureData[Taxonomy], rep_seqs: FeatureData[Sequence],
                   transposed_table: FeatureTable[Frequency]):
    table_m = transposed_table.view(view_type=Metadata)
    tax_m = taxonomy.view(view_type=Metadata)
    rep_m = rep_seqs.view(view_type=Metadata)

    merged = tax_m.merge(rep_m, table_m)
    visualization = tabulate(
        input=merged
    )

    return visualization