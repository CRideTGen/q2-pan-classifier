import os
import re
import tempfile
from pathlib import Path
import qiime2
from q2_types.per_sample_sequences import SingleLanePerSamplePairedEndFastqDirFmt, PairedEndSequencesWithQuality
from q2_types.sample_data import SampleData

def _get_forward_reverse_(name: str, file_names: list) -> tuple:
    pass
def _return_manifest_sample_string_(file_path_names: str) -> str:
    #TODO
    # 1. Change name from _return_names_
    # 2. error handle regex mismatched name
    # 3. finish _get_forward_reverse_ function, this will also check if name is unique
    # 4. Have function return string for all samples in manifest
    #
    PAIRED_END_REGEX = r'(.+)_.+_L[0-9][0-9][0-9]_R[12]_.+\.fastq.*'

    names_raw = [os.path.basename(x) for x in file_path_names]

    names_out = list()

    for full_name in names_raw:
        try:
            name = re.search(PAIRED_END_REGEX, full_name).group(1)
            names_out.append(name)
        except AttributeError:
            #TODO Handle mismatch case
            print(f"{full_name} did not match fastq naming scheme")
    return names_out

def _generate_manifest_file_(sample_dir_path: str, manifest_file_dir: str) -> str:
    """Generates a manifest file in a temporary directory to be used in sequence reads upload"""

    HEADER = ['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath']

    #TODO

    sample_manifest_string = _return_manifest_sample_string_(sample_dir_path)

    with open(os.path.join(manifest_file_dir, "manifest"), "w") as manifest_file:

        manifest_file.write('\t'.join(HEADER) + '\n')
        #TODO Write output from  _return_manifest_sample_string_
        #manifest_file.write(sample_manifest_string)

    return os.path.join(manifest_file_dir)

def import_reads(sequence_directory: str) -> SingleLanePerSamplePairedEndFastqDirFmt:

    output_dir = SingleLanePerSamplePairedEndFastqDirFmt()

    sequence_directory = os.path.abspath(sequence_directory)
    temp_dir = tempfile.TemporaryDirectory()

    manifest_file_path, names = _generate_manifest_file_(sequence_directory, "/home/cridenour/PycharmProjects/q2-pan-classifier/test_data/")

    # manifest_file_path = "/home/cridenour/PycharmProjects/q2-pan-classifier/test_data/MANIFEST"
    #
    reads = qiime2.Artifact.import_data(type='SampleData[PairedEndSequencesWithQuality]',
                                        view=manifest_file_path,
                                        view_type='PairedEndFastqManifestPhred33V2')
    reads.export_data(output_dir.path)
    return output_dir
