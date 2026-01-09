import os.path

import numpy as np


class GeneratingParsumFile:
    """"
        This object can be used for creating the required parsum.dat file, consisting of the partition sums of the
        molecules and isotopologues for the sqlite database.
        It requires the q{i}.txt files obtained from https://hitran.org/docs/iso-meta/, then it merely reforms the data.
    """
    def __init__(self, molecule_name: str, afgl_codes: list, file_names_partition_sum: list, parsum_path=''):
        """
        :param molecule_name:       e.g. O2
        :param afgl_codes:          see https://hitran.org/docs/iso-meta/, e.g. [626, 636]
        :param file_names:
        :param parsum_path:
        :return:
        """
        if len(afgl_codes) != len(file_names_partition_sum):
            raise ValueError("afgl_codes and file_names must have the same length")

        if all([os.path.exists(fn) for fn in file_names_partition_sum]):
            pass
        else:
            raise FileExistsError(f'{file_names_partition_sum=} at least ones of these files does not exist')

        self.molecule_name = molecule_name
        self.afgl_codes = afgl_codes
        self.file_names = file_names_partition_sum
        self._parsum_path = parsum_path
        self._data = [np.loadtxt(file_name) for file_name in file_names_partition_sum]
        self.write_to_file()

    @property
    def header(self) -> str:
        _header = ''
        for afgl_code in self.afgl_codes:
            _header += f'\t{self.molecule_name}_{afgl_code}'
        return 'temperature' + _header + '\n'

    def write_to_file(self) -> None:
        with open(fr'{self._parsum_path}\parsum.dat', 'w') as file_data:
            file_data.write(self.header)

            line_data = ''
            for i, x0 in enumerate(self._data[0][:, 0]):
                line_data += str(x0)
                for data_j in self._data:
                    if i >= len(data_j[:, 0]):
                        break                   # first file is longer than j'th series, so continuing to next line/file
                    elif x0 != data_j[:, 0][i]:
                        raise ValueError(f"{x0=} != {self._data[1][:, 0][i]=}!!!")  # check x_0 with x_j
                    line_data += f'\t{data_j[:, 1][i]}'
                line_data += '\n'
            line_data = line_data[:-1]
            file_data.write(line_data)
        return


if __name__ == "__main__":
    x = GeneratingParsumFile('NH3', [4111, 5111], [r"C:\Users\steij\Documents\Promotie\data\database\raw-nh3\q45.txt",
                                                   r"C:\Users\steij\Documents\Promotie\data\database\raw-nh3\q46.txt"],
                             r'C:\Users\steij\Desktop')


