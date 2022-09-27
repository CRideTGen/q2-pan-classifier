import os
import re
import tempfile
from pathlib import Path
from typing import Protocol

import qiime2
from q2_types.per_sample_sequences import SingleLanePerSamplePairedEndFastqDirFmt, PairedEndSequencesWithQuality
from q2_types.sample_data import SampleData
class SampleManifest(Protocol):
    _sample_name: str
    _sample_path: SamplePath

    def get_sample_manifest_str(self):
        ...

class PairedReads:
    def __init__(self, sample_name, sample_path=None, forward=None, reverse=None):
        self.sample_name = sample_name
        if sample_path and (forward or reverse):
            raise ValueError("Need sample_path or forward and reverse, NOT both")
        elif sample_path:
            self.sample_path = tuple(sample_path)
        elif forward and reverse:
            self.sample_path = tuple([forward, reverse])
        else:
            raise ValueError("Missing forward or reverse path")

    def get_paired_sample_name(self) -> str:
        PAIRED_END_REGEX = r'(.+)_.+_L[0-9][0-9][0-9]_R[12]_.+\.fastq.*'
        name = re.search(PAIRED_END_REGEX, self._sample_path.forward).group(1)

        return name

    def get_sample_manifest_str(self):
        return f"{self._sample_name}\t{self._sample_path.forward}\t{self._sample_path.reverse}\n"


class SingleRead:
    def __init__(self, sample_name, sample_path):
        if sample_path:
            self._sample_name = sample_name
            self._sample_path = sample_path
        else:
            raise ValueError("Missing forward or reverse path")

    def get_sample_manifest_str(self):
        pass


def _get_forward_reverse_(name: str, file_names: list) -> PairedReads:
    PAIRED_END_FORWARD_REGEX = f'{name}_.+_L[0-9][0-9][0-9]_R1_.+\.fastq.*'
    PAIRED_END_REVERSE_REGEX = f'{name}_.+_L[0-9][0-9][0-9]_R2_.+\.fastq.*'

    forward = [f_name for f_name in file_names if re.search(PAIRED_END_FORWARD_REGEX, f_name.name)]
    reverse = [r_name for r_name in file_names if re.search(PAIRED_END_REVERSE_REGEX, r_name.name)]

    if len(forward) == 1 and len(reverse) == 1:
        return PairedReads(sample_name=name, forward=forward, reverse=reverse)
    else:
        raise ValueError("Did not have matching pair")


def _return_manifest_sample_(sample_dir: Path) -> List[PairedReads]:
    # TODO
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

    return manifest_file_path

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


def seq_func(input_man: SampleManifest):
    return input_man.get_sample_manifest_str()


def main():
    x = PairedReads(sample_name="test", forward="./forward.fastq.gz", reverse="./reverse.fastq.gz")
    print(seq_func(input_man=x))


if __name__ == "__main__":
    main()
