import biom
import qiime2
import os
import numpy as np
import pandas as pd
from biom.table import Table
import glob
import tempfile
from Bio import SeqIO
import re


def filter_fastq_by_size(fastq_file_path: str, filter_size: int, sample_name: str) -> list:
    out_records = list()
    with open(fastq_file_path) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if len(record) >= filter_size:
                record.description = record.description.replace("no_sample", sample_name)
                out_records.append(record[0:filter_size])

    return out_records


def _get_sample_names_(barcode_names: list, metadata_file: str) -> list:
    metadata = pd.read_csv(metadata_file, sep='\t')
    sample_names = metadata["sampleid"]
    barcode_id = metadata["Barcode_ID"]

    bar_name_tmp = ['BC' + x.replace('barcode', '') for x in barcode_names]

    sample_out = list()

    sample_out += [sample_names[barcode_id.str.contains(bar)].astype(str).values.tolist()[0] for bar in bar_name_tmp]

    return sample_out


def _concatenate_fastq_(barcodes_directory: str, metadata_file: str = None) -> str:
    output_sequence_directory: str = os.path.join(barcodes_directory, "output_sequences")
    bar_dirs = glob.glob(os.path.join(barcodes_directory, "barcode*"))
    bar_names = [os.path.basename(x) for x in bar_dirs]
    if metadata_file:
        bar_tmp = [os.path.basename(x) for x in bar_dirs]
        sample_names = _get_sample_names_(barcode_names=bar_tmp, metadata_file=metadata_file)
    else:
        sample_names = [os.path.basename(x) for x in bar_dirs]
    combined_file_path = os.path.join(output_sequence_directory, "combined_all.fq")
    with open(combined_file_path, 'w') as combined_file:
        combined_seqs = list()
        for bar_name, sample_name in zip(bar_names, sample_names):
            seq_file_name = os.path.join(output_sequence_directory, sample_name + ".fastq")

            print(seq_file_name)

            with open(seq_file_name, 'w') as outfile:

                out_seqs = list()

                fastq_files = glob.glob(os.path.join(barcodes_directory, bar_name, "*.fastq"))

                for fastq in fastq_files:
                    out_seqs += filter_fastq_by_size(fastq, 16000, sample_name)

                combined_seqs += out_seqs

                SeqIO.write(out_seqs, outfile, "fastq")
        SeqIO.write(combined_seqs, combined_file, "fastq")

    return output_sequence_directory


def _generate_feature_table_(sequences_directory: str, path_to_combined_fasta: str) -> Table:
    sample_files = glob.glob(os.path.join(sequences_directory, "*.fastq"))
    sample_names = [(os.path.basename(i)).replace(".fastq", "") for i in sample_files]
    sequence_data = dict()

    with open(path_to_combined_fasta, "r") as combined_file:

        for record in SeqIO.parse(combined_file, "fastq"):
            sequence_data.update({record.id: 0})

    sequence_names = [i for i in sequence_data.keys()]

    table_data = []

    for indx, sample_file in enumerate(sample_files):
        sd_tmp_dict = sequence_data.copy()
        with open(sample_file, 'r') as handle:
            for record in SeqIO.parse(handle, 'fastq'):
                sd_tmp_dict[record.id] = 1

            sd_tmp_values = [i for i in sd_tmp_dict.values()]
            table_data.append(sd_tmp_values)

    seq_names_array = np.array(sequence_names)
    samp_names_array = np.array(sample_names)
    td_array = np.transpose(np.array(table_data))

    table_out = Table(data=td_array, observation_ids=seq_names_array, sample_ids=samp_names_array,
                      table_id="minion_frequency_table")

    return table_out


def _minion_generate_manifest_file_(barcode_directory: str, manifest_file_dir: str, metadata_file: str = None) -> None:
    HEADER = ['sample-id', 'absolute-filepath']
    sample_directory = _concatenate_fastq_(barcodes_directory=barcode_directory, metadata_file=metadata_file)

    sample_paths = glob.glob(os.path.join(sample_directory, '*.fastq'))
    sample_names = [(os.path.basename(x)).replace(".fastq", "").strip() for x in sample_paths]

    manifest_file = os.path.join(manifest_file_dir, 'MANIFEST')

    with open(manifest_file, 'w') as man_file:
        man_file.write('\t'.join(HEADER) + '\n')
        for out in zip(sample_names, sample_paths):
            man_file.write('\t'.join(out) + '\n')


def minion_sequence_import(barcode_directory: str):
    pass


def isonclust_feature_table(isonclust_directory: str, sample_names: list) -> Table:
    clustered_final_file = os.path.join(isonclust_directory, 'final_clusters.tsv')
    cluster_origins_file = os.path.join(isonclust_directory, 'final_cluster_origins.tsv')

    cluster_sequence_names = list()

    with open(cluster_origins_file, 'r') as co_file:
        line = co_file.readline()
        while line:
            cluster_sequence_names.append(line.split('\t')[0])
            line = co_file.readline()

    clust_seq_dict = {i: 0 for i in cluster_sequence_names}
    sample_clust_dict = {i: clust_seq_dict.copy() for i in sample_names}

    with open(clustered_final_file, 'r') as cf:
        line = cf.readline()
        while line:
            line_split = line.split('\t')
            clust_seq_id = line_split[0]
            sample_id_all = line_split[1]
            match = re.search('sampleid=([0-9A-Z]{4})_', sample_id_all)
            sample_id = match.group(1)

            sample_clust_dict[sample_id][clust_seq_id] += 1

            line = cf.readline()

    data_list = list()
    for k in sample_clust_dict:
        for i in sample_clust_dict[k].values():
            data_list.append(i)

    seq_names_array = np.array(cluster_sequence_names)
    samp_names_array = np.array(sample_names)
    td_array = np.transpose(np.array(data_list).reshape(len(sample_names), len(cluster_sequence_names)))

    table_out = Table(data=td_array, observation_ids=seq_names_array, sample_ids=samp_names_array,
                      table_id="isonclust_frequency_table")

    return table_out


def isonclust_representative_sequences(isonclust_directory: str, output_directory: str) -> None:
    cluster_origins_file = os.path.join(isonclust_directory, 'final_cluster_origins.tsv')

    rep_seqs_file_path = os.path.join(output_directory, 'isonclust-rep-seqs.fasta')

    with open(cluster_origins_file, 'r') as clust_file, open(rep_seqs_file_path, 'w') as rep_file:
        line = clust_file.readline()

        while line:
            line_split = line.split('\t')
            clust_seq_id = line_split[0]
            sequence = line_split[2]

            rep_file.write(f'>{clust_seq_id}\n{sequence}\n')

            line = clust_file.readline()


if __name__ == '__main__':
    import biom

    barcode_directory = '/home/chase/DissertationProjects/Qiime2SummerCamp/PanViral/CulexMitoGenome/sequence_reads/'
    output_sequence_directory = '/home/chase/DissertationProjects/Qiime2SummerCamp/PanViral/CulexMitoGenome/sequence_reads/output_sequences/'
    metadata_file = '/home/chase/DissertationProjects/Qiime2SummerCamp/PanViral/CulexMitoGenome/Qiime_MetaData_MIP01_09162021_ZAB.tsv'

    # _concatenate_fastq_(barcode_directory)
    # _minion_generate_manifest_file_(barcode_directory=barcode_directory, manifest_file_dir=barcode_directory, metadata_file=metadata_file)

    # combined_fasta = os.path.join(output_sequence_directory, 'combined_all.fas')
    # t_out = _generate_feature_table_(output_sequence_directory, combined_fasta)
    #
    # with biom.util.biom_open(os.path.join(output_sequence_directory, "minion_table.biom"), 'w') as ff:
    #     t_out.to_hdf5(ff, "minion frequency table")

    isonclust_directory = '/home/chase/DissertationProjects/Qiime2SummerCamp/PanViral/CulexMitoGenome/sequence_reads/output_sequences/cluster_test'

    isonclust_representative_sequences(isonclust_directory=isonclust_directory,
                                       output_directory=output_sequence_directory)
    isonclut_table = isonclust_feature_table(isonclust_directory=isonclust_directory,
                                             sample_names=['0M32', '0M43', '12G8', '2M49', '54G7', '63D7', '7Y91',
                                                           '8Z26', '9Y89'])

    with biom.util.biom_open(os.path.join(output_sequence_directory, "isonclust_minion_table.biom"), 'w') as ff:
        isonclut_table.to_hdf5(ff, "minion frequency table")
