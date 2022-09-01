import os

import pandas as pd
import qiime2


#Prep Sequences Action Helpers
def _return_names_(file_path_names: list) -> list:
    names_raw = [os.path.basename(x) for x in file_path_names]

    names_out = list()

    for name in names_raw:
        if '-' in name:
            name_split = name.split('-')
            name_1_len = len(name_split[0])

            if name_1_len < 8:
                names_out.append(name_split[1])
            else:
                names_out.append(name_split[0])
        else:
            names_out = name.split('_jk')[0]
    return names_out

def _generate_manifest_file_(sample_dir_path: str, manifest_file_dir: str) -> tuple:
    """Generates a manifest file in a temporary directory to be used in sequence reads upload"""

    HEADER = ['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath']

    sequence_paths = os.listdir(sample_dir_path)
    forward_paths = [os.path.join(sample_dir_path, x) for x in sequence_paths if "_R1" in x]
    reverse_paths_tmp = [os.path.join(sample_dir_path, x) for x in sequence_paths if "_R2" in x]

    names = _return_names_(forward_paths)

    reverse_paths = []

    for name in names:
        for rev in reverse_paths_tmp:
            if name in rev:
                reverse_paths.append(rev)

    with open(os.path.join(manifest_file_dir, "manifest"), "w") as manifest_file:

        manifest_file.write('\t'.join(HEADER) + '\n')

        for out in zip(names, forward_paths, reverse_paths):
            manifest_file.write('\t'.join(out) + '\n')

    return tuple([os.path.join(manifest_file_dir, "manifest"), names])


def _write_metadata_template_(sample_names: list, output_dir: str) -> None:
    HEADER = "sample-id"

    if output_dir and not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    elif not output_dir:
        output_dir = os.getcwd()

    with open(os.path.join(output_dir, "metadata_template.tsv"), "w") as metadata_file:
        metadata_file.write(HEADER + '\n')
        metadata_file.write('\n'.join(sample_names))

# Classify Reads action helpers

def _merge_table_(transpose, dada2_table_out, dada2_rep_seqs_out, classified):
    # transposing table and getting metadata
    tt, = transpose(table=dada2_table_out)
    tt_dt = tt.view(pd.DataFrame)
    tt_dt["Total"] = tt_dt.sum(axis=1)
    tt_dt.index.name = "feature id"

    tt_m = qiime2.Metadata(tt_dt)
    dr_m = dada2_rep_seqs_out.view(view_type=qiime2.Metadata)
    c_m = classified.view(view_type=qiime2.Metadata)

    merged_table = c_m.merge(dr_m, tt_m)

    return merged_table

def _get_cpus() -> int:
    num_cpus = os.cpu_count()

    if num_cpus > 2:
        return num_cpus - 2
    else:
        return num_cpus - 1