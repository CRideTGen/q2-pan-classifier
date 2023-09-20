import os
import tempfile

import qiime2
from q2_pan_classifier.utils.uils import (
    _get_cpus,
    _generate_manifest_file_,
    _write_metadata_template_,
)


def prep_sequence_reads_paired(ctx, sequences=None, sequences_directory=None, metadata_template_dir=None, primer_f=None,
                               primer_r=None):
    results = []

    # importing external plugins to be used later

    create_table_viz = ctx.get_action('demux', 'summarize')
    cut_adapt = ctx.get_action('cutadapt', 'trim_paired')

    # importing sequences

    if sequences and sequences_directory:
        raise ValueError("Please only provide a sequence qza or sequence_directory path not both.")
    elif sequences:
        read_seqs = sequences

    elif sequences_directory:

        sequences_directory = os.path.abspath(sequences_directory)
        temp_dir = tempfile.TemporaryDirectory()

        manifest_file_path, names = _generate_manifest_file_(sequences_directory, temp_dir.name)

        read_seqs = qiime2.Artifact.import_data(type='SampleData[PairedEndSequencesWithQuality]',
                                                view=manifest_file_path,
                                                view_type='PairedEndFastqManifestPhred33V2')

        # write metadata template
        _write_metadata_template_(sample_names=names, output_dir=metadata_template_dir)

    else:
        raise ValueError("Please provide either a sequence qza or sequence directory path"
                         )

    # using plugins to trim reads and create reads visualization

    if primer_f and primer_r:
        trimmed_reads = cut_adapt(demultiplexed_sequences=read_seqs,
                                  front_f=[primer_f],
                                  front_r=[primer_r],
                                  cores=_get_cpus())

        table_viz = create_table_viz(data=trimmed_reads.trimmed_sequences)

        results += trimmed_reads

    else:
        table_viz = create_table_viz(data=read_seqs.trimmed_sequences)
        results += [read_seqs]

    results += table_viz

    return tuple(results)


def prep_sequence_reads_single(ctx, sequences=None, sequences_directory=None, metadata_template_dir=None, primer_f=None,
                               primer_r=None):
    results = []

    # importing external plugins to be used later
    cut_adapt = ctx.get_action('cutadapt', 'trim_single')
    create_table_viz = ctx.get_action('demux', 'summarize')

    # importing sequences

    if sequences and sequences_directory:
        raise ValueError("Please only provide a sequence qza or sequence_directory path not both.")
    elif sequences:
        read_seqs = sequences

    elif sequences_directory:

        sequences_directory = os.path.abspath(sequences_directory)
        temp_dir = tempfile.TemporaryDirectory()

        manifest_file_path, names = _generate_manifest_file_(sequences_directory, temp_dir.name)

        read_seqs = qiime2.Artifact.import_data(type='SampleData[PairedEndSequencesWithQuality]',
                                                view=manifest_file_path,
                                                view_type='PairedEndFastqManifestPhred33V2')

        # write metadata template
        _write_metadata_template_(sample_names=names, output_dir=metadata_template_dir)

    else:
        raise ValueError("Please provide either a sequence qza or sequence directory path"
                         )

    # using plugins to trim reads and create reads visualization

    if primer_f and primer_r and sequences:
        trimmed_reads = cut_adapt(demultiplexed_sequences=read_seqs,
                                  front=[primer_f, primer_r],
                                  cores=_get_cpus())

        table_viz = create_table_viz(data=trimmed_reads.trimmed_sequences)

        results += trimmed_reads
    else:
        table_viz = create_table_viz(data=read_seqs.trimmed_sequences)
        results += [read_seqs]

    results += table_viz

    return tuple(results)
