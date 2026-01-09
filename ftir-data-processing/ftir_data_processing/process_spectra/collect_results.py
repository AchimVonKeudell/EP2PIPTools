import os.path
from dataclasses import dataclass
import numpy as np
from os import listdir
from configparser import ConfigParser
import matplotlib.pyplot as plt
import brukeropusreader
from ..file_extension import *


def create_timing_dat(directory_raw_data: str, directory_output_data: str, **kwargs) -> None:
    """Creates timing.dat containing the time stamps of the data points of a specific data set

    :param directory_raw_data:
    :param directory_output_data:
    :return:
    """

    output_full_path_name = f"{directory_output_data}\\timing.dat"

    def to_minutes(time_string: str):
        return 60 * int(time_string[:2]) + int(time_string[3:5]) + float(time_string[6:]) / 60

    background_time = None

    output = []
    for file in listdir(directory_raw_data):
        index = file.split('.')[-1]

        try:
            opus_data = brukeropusreader.read_file(directory_raw_data + '\\' + file)
        except PermissionError:
            continue

        if 'AB' not in opus_data:
            continue            # this is a background file, so lets skip it

        if ('background_file' in kwargs) & (background_time is None):
            background_file = kwargs['background_file']
            background_data = brukeropusreader.read_file(background_file)
            background_time = background_data['ScSm Data Parameter']['TIM'].split()[0]
        else:
            background_time = opus_data['ScRf Data Parameter']['TIM'].split()[0]

        spectrum_time = opus_data['AB Data Parameter']['TIM'].split()[0]
        corrected_time = to_minutes(spectrum_time) - to_minutes(background_time)
        output.append([index, corrected_time])

    np.savetxt(output_full_path_name, output, header='index\ttime [min]',
               delimiter='\t', fmt='%s')
    print(f'\tcreated {output_full_path_name}')
    return


@dataclass
class Output:
    directory_results: str
    data_point_indexes: any
    time: any
    n: any

    def omit(self, excluding_index: any) -> None:
        """Omits certain data points according to the given indexes.

        :param excluding_index: string values of the data_point_indexes to be omitted
        :return:
        """
        booleans = [index not in excluding_index for index in self.data_point_indexes]

        self.data_point_indexes = np.array(self.data_point_indexes)[booleans]
        self.time = np.array(self.time)[booleans]
        self.n = np.array(self.n)[booleans]
        return

    def reordering(self) -> None:
        """Reorders the data according to their time stamps.

        :return:
        """
        self.time = np.array(self.time)
        ordering_array = self.time.argsort()

        self.data_point_indexes = np.array(self.data_point_indexes)[ordering_array]
        self.time = self.time[ordering_array]
        self.n = np.array(self.n)[ordering_array]
        return

    def select(self, t_min=0, t_max=1e99):
        """Returns the same object containing the data enclosed by tmin and tmax.

        :param t_min: in minutes
        :param t_max: in minutes
        :return: Output
        """
        selecting_array = (self.time > t_min) & (self.time < t_max)

        return Output(
            directory_results=self.directory_results,
            data_point_indexes=self.data_point_indexes[selecting_array],
            time=self.time[selecting_array],
            n=self.n[selecting_array]
        )

    def __add__(self, other):
        def combine(attribute):
            if isinstance(getattr(self, attribute), list):
                return getattr(self, attribute) + getattr(other, attribute)
            elif isinstance(getattr(self, attribute), np.ndarray):
                return np.append(getattr(self, attribute), getattr(other, attribute))
            raise TypeError(f"{type(attribute)=} not known/implemented")

        return Output(
            directory_results=self.directory_results,
            data_point_indexes=combine('data_point_indexes'),
            time=combine('time'),
            n=combine('n')
        )

    def __getitem__(self, key) -> (float, float):
        return self.time[key], self.n[key]


def get_values_over_time(
        output_fitting_method_path: str,
        timing_data_filepath: str,
        parameter_name: str,
) -> Output:
    """
    Obtains the molar fraction over time for the various measurements.

    :param output_fitting_method_path:
    :param timing_data_filepath:
    :param parameter_name:
    :return:
    """

    print(f'{timing_data_filepath=}')
    timing_data = np.genfromtxt(timing_data_filepath, dtype=str)

    def get_ftir_timing() -> dict:
        """Retrieves the time stamp

        :return: dictionary with indexes and time stamps
        """
        return {index: float(timing) for index, timing in timing_data}

    output_values = Output(output_fitting_method_path, [], [], [])
    ftir_timing_dict = get_ftir_timing()
    for file in listdir(output_fitting_method_path):
        if '.metadata' in file:

            data_point_index = os.path.splitext(os.path.splitext(file)[0])[1][1:]
            if 'background' in data_point_index:
                continue
            output_values.data_point_indexes.append(data_point_index)

            # set time stamp
            output_values.time.append(ftir_timing_dict[data_point_index])

            # obtain FitSpectrum results
            results = ConfigParser()
            results.read(fr"{output_fitting_method_path}\{file}")

            string_value = results.get("Parameters", parameter_name)
            # elif value_type == 'uncertainties':
            #     string_value = results.get("Variables", 'uncertainty')
            output_values.n.append(float(string_value))

    output_values.reordering()
    return output_values


def collect_results(
        directory_output_data: str,
        output_directory: str,
        molecule_label: str,
        parameter_name: str,
):
    print(f'Collect results: {parameter_name}')
    data = get_values_over_time(directory_output_data, output_directory + '\\timing.dat', parameter_name)

    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111)
    ax.grid(True)
    ax.set(xlabel='time / min', ylabel=parameter_name)
    ax.plot(data.time, data.n, marker='o', ls=':')
    fig.tight_layout()

    output_fn = f'{output_directory}\\{molecule_label}_{parameter_name}.{DATA_EXTENSION}'
    np.savetxt(output_fn, np.stack([data.time, data.n], axis=1),
               delimiter='\t', header=f'time/min\t{parameter_name}')
    fig.savefig(f'{output_fn[:-4]}.{FIGURE_EXTENSION}')
    return


def collect_results_incl_uncertainty(directory_output_data: str, specie: str, parameter_name='molar_fraction'):
    plt.rc('axes', labelweight='bold')

    concentration_values = get_values_over_time(directory_output_data, specie, value_type='fits',
                                                parameter_name=parameter_name)
    uncertainty_values = get_values_over_time(directory_output_data, specie, value_type='uncertainties',
                                              parameter_name=parameter_name)

    fig, ax = plt.subplots(1, 1)
    ax.plot(concentration_values.time, concentration_values.n, marker='o', ls=':')
    ax.fill_between(
        concentration_values.time,
        concentration_values.n - uncertainty_values.n, concentration_values.n + uncertainty_values.n,
        color='tab:blue', alpha=0.5
    )

    ax.set(xlabel='time [min]', ylabel=f'{parameter_name} [-]')
    ax.grid(True)
    np.savetxt(directory_output_data + fr'\{specie.upper()}\{parameter_name}_incl_uncertainty.txt',
               np.transpose((concentration_values.time, concentration_values.n, uncertainty_values.n)),
               delimiter='\t', header=f'time [min]\t{parameter_name}\tuncertainty {parameter_name}')
    fig.savefig(directory_output_data + fr'\{specie.upper()}\{parameter_name}_incl_uncertainty.png')
    return

