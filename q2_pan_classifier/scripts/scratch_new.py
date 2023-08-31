import argparse

import click
import qiime2
import yaml
from qiime2 import Metadata

from q2_pan_classifier.scripts.script_methods import (
    import_viral_sequences,
    trim_sequences,
    visualize_sequence_quality,
    cluster_sequences,
    visualize_feature_table,
    classify_sequences,
    create_taxa_bar_plot,
    transpose_table,
    tabulate_table,
)

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--config", required=True)

args = parser.parse_args()


def main():
    print("hello")
    with open(args.config, "r") as cf:
        config_loaded = yaml.safe_load(cf)
        manifest_file_path = config_loaded["manifest_path"]
        metadata_path = config_loaded["metadata_path"]
        metadata = Metadata.load(metadata_path)
        classifier_path = config_loaded["classifier_path"]
        classifier = qiime2.Artifact.load(classifier_path)
        forward_primer = config_loaded["forward_primer"]
        reverse_primer = config_loaded["reverse_primer"]
        trunc_len_f = config_loaded["trunc_len_f"]
        trunc_len_r = config_loaded["trunc_len_r"]

        # make the output directories if they do not exist
        # if not os.path.isdir(output_dir):
        #    os.mkdir(output_dir)
        #    os.mkdir(visual_output_dir)

        viral_sequences = import_viral_sequences(manifest_file_path)
        # click.echo(viral_sequences)
        viral_sequences.save("./output/viral_sequences")

        # TODO: log the ones that are trimmed
        # will go into the database
        (trimmed_sequences,) = trim_sequences(
            viral_sequences, [forward_primer], [reverse_primer]
        )
        # click.echo(trimmed_sequences)
        trimmed_sequences.save("./output/paired_end_demux_trimmed")

        visual_out = visualize_sequence_quality(trimmed_sequences)
        # click.echo(visual_out.visualization)
        visual_out.visualization.save("./output/demux_visual_out")

        # check if the truncate length forward and reverse are given, run if they are.
        if trunc_len_r != "" and trunc_len_f != "":
            table, rep_seqs, stats = cluster_sequences(
                trimmed_sequences, int(trunc_len_f), int(trunc_len_r)
            )
            # click.echo(table)
            # click.echo(rep_seqs)
            # click.echo(stats)
            table.save("./output/table_dada2")
            rep_seqs.save("./output/representative_sequences")
            stats.save("./output/stats_dada2")

            visualization, tabulate = visualize_feature_table(table, rep_seqs, metadata)
            # click.echo(visualization)
            # click.echo(tabulate)
            visualization.visualization.save("./output/table_summarize")
            tabulate.visualization.save("./output/tabulate_rep_seqs")

            (taxonomy,) = classify_sequences(classifier, rep_seqs)
            # click.echo(taxonomy)
            taxonomy.save("./output/taxonomy")
            #
            visualization_barplot = create_taxa_bar_plot(table, taxonomy, metadata)
            # #click.echo(visualization_barplot)
            visualization_barplot.visualization.save("./output/taxa_bar_plot")
            #
            #
            (transposed_table,) = transpose_table(table)
            # #click.echo(transposed_table)
            transposed_table.save("./output/transposed_table")
            #
            #
            tabulated_table = tabulate_table(taxonomy, rep_seqs, transposed_table)
            # #click.echo(tabulated_table)
            tabulated_table.visualization.save("./output/tabulated_table")

        else:
            click.echo(
                "\nTruncate length forward and reverse were not given.\n"
                "Analysis stopped at the clustering of sequences.\n"
            )

        ## TODO: Notes
        # create a output directory:
        #  visualations_output
        #  flag to save other earlier output


if __name__ == "__main__":
    main()
