import argparse
import subprocess
import yaml
import qiime2

from pathlib import Path
import qiime2.plugins.dada2.actions as dada2_actions
from q2_feature_classifier._taxonomic_classifier import TaxonomicClassifier
from q2_types.feature_data import FeatureData, Sequence, Taxonomy
from q2_types.feature_table import FeatureTable, Frequency
from qiime2.plugins.cutadapt.methods import trim_paired
from q2_types.per_sample_sequences import PairedEndSequencesWithQuality
from qiime2.plugins.demux.visualizers import summarize as demux_summarize
from q2_types.sample_data import SampleData
from qiime2.plugins.feature_table.visualizers import summarize as table_summarize
from qiime2 import Metadata
from qiime2.plugins.feature_table.visualizers import tabulate_seqs
from qiime2.plugins.feature_classifier.methods import classify_sklearn
from qiime2.plugins.taxa.visualizers import barplot
from qiime2.plugins.feature_table.methods import transpose
from qiime2.plugins.metadata.visualizers import tabulate



parser = argparse.ArgumentParser()

parser.add_argument("-c", '--config', required=True)

args = parser.parse_args()

def main():
    with open(args.config, 'r') as cf:
        config_loaded = yaml.safe_load(cf)
        manifest_file_path = config_loaded['manifest_path']
        metadata_path = config_loaded['metadata_path']
        metadata = Metadata.load(metadata_path)
        classifier_path = config_loaded['classifier_path']
        classifier = qiime2.Artifact.load(classifier_path)
        forward_primer = config_loaded['forward_primer']
        reverse_primer = config_loaded['reverse_primer']
        trunc_len_f = config_loaded['trunc_len_f']
        trunc_len_r = config_loaded['trunc_len_r']

        # make the output directories if they do not exist
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
            os.mkdir(visual_output_dir)

        viral_sequences = import_viral_sequences(manifest_file_path)
        #click.echo(viral_sequences)
        #viral_sequences.save(output_dir + "/visual_output/viral_sequences.qza")


        # TODO: log the ones that are trimmed
        # will go into the database
        trimmed_sequences, = trim_sequences(viral_sequences, [forward_primer], [reverse_primer])
        click.echo(trimmed_sequences)
        #trimmed_sequences.save(output_dir + "/visual_output/paired_end_demux_trimmed.qza")


        visual_out = visualize_sequence_quality(trimmed_sequences)
        click.echo(visual_out)
        #visual_out.save(output_dir + "/visual_output/demux_visual_out.qza")

        # check if the truncate length forward and reverse are given, run if they are.
        if(trunc_len_r != "" and trunc_len_f != ""):
            table, rep_seqs, stats = cluster_sequences(trimmed_sequences, int(trunc_len_f), int(trunc_len_r))
            click.echo(table)
            click.echo(rep_seqs)
            click.echo(stats)
            #table.save(output_dir + "/visual_output/table_dada2.qza")
            #rep_seqs.save(output_dir + "/visual_output/representative_sequences.qza")
            #stats.save(output_dir + "/visual_output/stats_dada2.qza")


            visualization, tabulate = visualize_feature_table(table, rep_seqs, metadata)
            click.echo(visualization)
            click.echo(tabulate)
            #visualization.save(output_dir + "/visual_output/table_summarize.qza")
            #tabulate.save(output_dir + "/visual_output/tabulate_rep_seqs.qza")


            taxonomy, = classify_sequences(classifier, rep_seqs)
            click.echo(taxonomy)
            #taxonomy.save(output_dir + "/visual_output/taxonomy.qza")

            visualization_barplot = create_taxa_bar_plot(table, taxonomy, metadata)
            click.echo(visualization_barplot)
            #visualization_barplot.save(output_dir + "/visual_output/taxa_bar_plot.qza")


            transposed_table, = transpose_table(table)
            click.echo(transposed_table)
            #transposed_table.save(output_dir + "/visual_output/transposed_table.qza")


            tabulated_table = tabulate_table(taxonomy, rep_seqs, transposed_table)
            click.echo(tabulated_table)
            #tabulated_table.save(output_dir + "/visual_output/tabulated_table.qza")


        else:
            click.echo("\nTruncate length forward and reverse were not given.\n"
                       "Analysis stopped at the clustering of sequences.\n")

        ### TODO: Notes
        # create a output directory:
        #   visualations_output
        #   flag to save other earlier output

if __name__ == "__main__":
    main()

