from pathlib import Path
import gzip
from q2_types.per_sample_sequences import SingleLanePerSamplePairedEndFastqDirFmt


def import_reads(sequence_directory: str) -> SingleLanePerSamplePairedEndFastqDirFmt:
    output_directory = SingleLanePerSamplePairedEndFastqDirFmt()

    with gzip.open(Path(output_directory.path).joinpath("a_a_L111_R1_001.fastq.gz"), "wb") as ff:
        with open("/home/cridenour/PycharmProjects/q2-pan-classifier/test_data/ENTV_05Mar2022/trimmed_sequences/ENTV87Skyview20220228_S29_L001_R1_output.fastq", "rb") as fastq_file:
            ff.writelines(fastq_file)
    with gzip.open(Path(output_directory.path).joinpath("a_a_L111_R2_001.fastq.gz"), "wb") as ff:
        with open(
                "/home/cridenour/PycharmProjects/q2-pan-classifier/test_data/ENTV_05Mar2022/trimmed_sequences/ENTV87Skyview20220228_S29_L001_R2_output.fastq", "rb") as fastq_file:
            ff.writelines(fastq_file)
    with open(Path(output_directory.path).joinpath("MANIFEST"), "w") as ff:
        ff.write(','.join(['sample-id', 'filename', 'direction\n']))
        ff.write(','.join(['ENTV87Skyview20220228', '/home/cridenour/PycharmProjects/q2-pan-classifier/test_data/ENTV_05Mar2022/trimmed_sequences/ENTV87Skyview20220228_S29_L001_R1_output.fastq', 'forward\n']))
        ff.write(','.join(['ENTV87Skyview20220228', '/home/cridenour/PycharmProjects/q2-pan-classifier/test_data/ENTV_05Mar2022/trimmed_sequences/ENTV87Skyview20220228_S29_L001_R2_output.fastq', 'reverse\n']))

    with open(Path(output_directory.path).joinpath("metadata.yml"), "w") as ff:
        ff.write(','.join(['sample-id']))
    return output_directory