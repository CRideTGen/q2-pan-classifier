from pathlib import Path, PurePath
from Bio import SeqIO
import pandas as pd
import argparse

class MinionReads:

    def __init__(self, barcode_directory_path: str, output_directory: str, min_length: int,
                 metadata_filepath: str = None):

        self.bar_dir = Path(barcode_directory_path)
        self.min_length = min_length
        self.metadata_filepath = Path(metadata_filepath)
        self.output_directory = Path(output_directory)

        if not self.bar_dir.is_dir():
            raise NotADirectoryError(f'{self.bar_dir} is not a valid directory path')

        if not self.metadata_filepath.is_file():
            raise FileNotFoundError(f'{self.metadata_filepath} does not exist')

        Path.mkdir(self.output_directory, parents=False, exist_ok=True)

    def _filter_fastq_by_size_(self, fastq_file_path: PurePath, filter_size: int, sample_name: str) -> list:
        """

        :param filter_size:
        :param sample_name:
        :return: list of SeqIO sequence records
        """
        out_records = list()
        with open(fastq_file_path) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                if len(record) >= filter_size:
                    record.description = record.description.replace("no_sample", sample_name)
                    out_records.append(record[0:filter_size])

        return out_records

    def _get_sample_names_(self, barcode_names: list, metadata_file: Path) -> list:
        """

        :param metadata_file:
        :return:
        """
        metadata = pd.read_csv(metadata_file, sep='\t')
        sample_names = metadata["sampleid"]
        barcode_id = metadata["Barcode_ID"]

        bar_name_tmp = ['BC' + x.replace('barcode', '') for x in barcode_names]

        sample_out = list()

        sample_out += [sample_names[barcode_id.str.contains(bar)].astype(str).values.tolist()[0] for bar in
                       bar_name_tmp]

        return sample_out

    def concatenate_fastq(self) -> None:
        """

        :return:
        """
        print("starting concatenation")
        bar_dirs = [bar for bar in self.bar_dir.glob("barcode*")]
        bar_names = [PurePath(bar_name).name for bar_name in bar_dirs]

        if self.metadata_filepath:
            bar_tmp = [PurePath(x).name for x in bar_dirs]
            sample_names = self._get_sample_names_(barcode_names=bar_tmp, metadata_file=self.metadata_filepath)
        else:
            sample_names = [PurePath(x).name for x in bar_dirs]
        combined_file_path = PurePath(self.output_directory).joinpath("combined_all.fq")

        with open(combined_file_path, 'w') as combined_file:
            combined_seqs = list()
            for bar_name, sample_name in zip(bar_names, sample_names):
                seq_file_name = PurePath(self.output_directory).joinpath(sample_name + ".fastq")

                print(seq_file_name)

                with open(seq_file_name, 'w') as outfile:

                    out_seqs = list()
                    bar_dir_path = PurePath.joinpath(self.bar_dir, bar_name)
                    fastq_files = bar_dir_path.glob("*.fastq")

                    for fastq in fastq_files:
                        out_seqs += self._filter_fastq_by_size_(fastq, self.min_length, sample_name)

                    combined_seqs += out_seqs

                    SeqIO.write(out_seqs, outfile, "fastq")
            SeqIO.write(combined_seqs, combined_file, "fastq")





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenates and filters minion sequences")

    parser.add_argument("-d", "--barcode_directory_path",
                        dest="barcode_directory",
                        type=str,
                        help="path to directory which holds the barcode directories")
    parser.add_argument("-n", "--min_length",
                        type=int,
                        dest="min_length",
                        help="minimum length of read to keep")

    parser.add_argument('-m', '--metadata_filepath',
                        dest="metadata_filepath",
                        type=str,
                        help="Path to metadata file")

    parser.add_argument('-o', '--outpout_directory',
                        dest="output_directory",
                        type=str,
                        help="Path to directory to store concatenated sequences")

    args = parser.parse_args()

    minion_obj = MinionReads(barcode_directory_path=args.barcode_directory,
                             min_length=args.min_length,
                             metadata_filepath=args.metadata_filepath,
                             output_directory=args.output_directory)



    minion_obj.concatenate_fastq()

    print("Done")

