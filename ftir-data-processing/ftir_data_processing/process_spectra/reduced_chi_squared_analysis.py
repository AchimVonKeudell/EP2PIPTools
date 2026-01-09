import numpy as np
import matplotlib.pyplot as plt
import configparser, os
from lmfit import Parameters
from dataclasses import dataclass

from vibrationalspectrum import GetSimulator, GetDistibution, GetSpectrum
from .paths import get_database
from .uncertainty_estimation import uncertainty_analysis


def get_specie(specie_name) -> str:
    """Converts the specie_name to an actual specie's name. For instance, if the v1+v2 transition is fitted for
    CO2, the specie name is 'co2-v1+v2'. Then this function converts is back to 'co2'.
    """
    if '-' in specie_name:
        return specie_name.split('-')[0]
    return specie_name


@dataclass
class FitData:
    x: np.ndarray
    transmittance: np.ndarray
    fit: np.ndarray

    spectrum_file_name: str
    number_of_variables = None
    thermal_equilibrium_spectrum: bool

    parameters: Parameters = None

    @staticmethod
    def from_txt_file(full_file_name: str, **kwargs):
        if os.path.exists(full_file_name):
            data = np.loadtxt(full_file_name, delimiter=',')
            return FitData(
                x=data[:, 0],
                transmittance=data[:, 1],
                fit=data[:, 2],
                spectrum_file_name=full_file_name,
                **kwargs
            )
        raise FileExistsError(f"{full_file_name=} not found")

    def read_metadata(self, specie: str, *, section_name='Parameters', **kwargs):

        def get_metadata(parameter_name):
            if parameter_name in kwargs:
                return kwargs[parameter_name]
            return float(metadata[section_name][parameter_name])

        metadata = configparser.ConfigParser()
        metadata.read(self.spectrum_file_name[:-4] + '.metadata')

        def to_bool(value: str) -> bool:
            if value == 'yes':
                return True
            elif value == 'no':
                return False
            raise ValueError

        print(self.spectrum_file_name[:-4] + '.metadata')
        self.number_of_variables = sum([to_bool(item) for item in metadata['Fixed'].values()])

        self.parameters = Parameters()
        self.parameters.add("pressure_atm", value=get_metadata('pressure_atm'), vary=False)
        self.parameters.add("path_length_cm", value=get_metadata('path_length_cm'), vary=False)
        self.parameters.add('wavenumber_offset', value=get_metadata('wavenumber_offset'), vary=False)
        self.parameters.add("instrumental_broadening_sigma", value=get_metadata('instrumental_broadening_sigma'),
                            min=0., max=1.0, vary=False)
        self.parameters.add(f"fraction_{specie}", value=get_metadata(f'fraction_{specie}'), min=0, max=1)
        self.parameters.add('temperature_rotation_kelvin', value=get_metadata('temperature_rotation_kelvin'),
                            min=200, max=1500, vary=True)
        if not self.thermal_equilibrium_spectrum:
            self.parameters.add('temperature_vib_12_kelvin', value=get_metadata('temperature_vib_12_kelvin'),
                                min=200, max=5000)
            self.parameters.add('temperature_vib_3_kelvin', value=get_metadata('temperature_vib_12_kelvin'),
                                min=200, max=5000)
        i = 0
        while f"p{i}" in metadata['Parameters']:
            self.parameters.add(f'p{i}', value=get_metadata(f'p{i}'), vary=False)
            i += 1
        return

    def adjust_fraction(self, specie: str, **kwargs):

        def get_value(parameter_name):
            if parameter_name in kwargs:
                return kwargs[parameter_name]
            return self.parameters[parameter_name].value

        parameters = Parameters()
        parameters.add("pressure_atm", value=get_value('pressure_atm'), vary=False)
        parameters.add("path_length_cm", value=get_value('path_length_cm'), vary=False)
        parameters.add('wavenumber_offset', value=get_value('wavenumber_offset'), vary=False)
        parameters.add("instrumental_broadening_sigma", value=get_value('instrumental_broadening_sigma'),
                       min=0., max=1.0, vary=False)
        parameters.add(f"fraction_{specie}", value=get_value(f'fraction_{specie}'), min=0, max=1)
        parameters.add('temperature_rotation_kelvin', value=get_value('temperature_rotation_kelvin'),
                       min=200, max=1500, vary=True)
        if not self.thermal_equilibrium_spectrum:
            parameters.add('temperature_vib_12_kelvin', value=get_value('temperature_vib_12_kelvin'),
                                min=200, max=5000)
            parameters.add('temperature_vib_3_kelvin', value=get_value('temperature_vib_12_kelvin'),
                                min=200, max=5000)
        i = 0
        while f"p{i}" in self.parameters:
            parameters.add(f'p{i}', value=get_value(f'p{i}'), vary=False)
            i += 1
        return parameters


def uncertainty_analysis_per_time_stamp(
        output_data_path: str, timing_index: str, specie_band_name: str, internal_error: float, resolution: float, *,
        parameter_name='fraction', equilibrium_distribution=True, verbose=False
):
    """Estimates the uncertainty using the reduced chi-squared method [1], where the best fit will result in a minimal
    chi-squared value. Then, the uncertainty is estimated for parameter X by solving [2]
        red_chi-squared(X_best_fit + X_error) = red_chi-squared(X_best_fit) + 1

    The criterion `+ 1` is chosen such that the change in the calculated spectrum by the uncertainty in the molar
    fraction equals the noise - i.e. internal error - of the measured spectrum.

    [1]: Drosg M (2007), Dealing with Uncertainties (Berlin: Springer)
    [2]: M A Damen et al. 2020 Plasma Sources Sci. Technol. 29 065016

    :param output_data_path:
    :param timing_index:
    :param specie_band_name:    can be the same as the specie name, but doesn't have to be. E.g. co2-v1+v3
    :param internal_error:      estimation of the noise on the transmittance spectrum obtained by the FTIR
    :param resolution:          of the uncertainty of the molar fraction
    :param parameter_name:
    :param equilibrium_distribution:
    :param verbose:
    :return:
    """

    def get_parameter_name():
        if parameter_name == 'fraction':
            return f'fraction_{specie_name}'
        elif parameter_name == 'rotational_temperature':
            return f'temperature_rotation_kelvin'
        else:
            raise ValueError(f"{parameter_name=} not implemented")

    # 0. define relevant full file paths
    file_name_best_fit_spectrum = output_data_path + fr'\{specie_band_name}\fits\{timing_index}.txt'
    uncertainty_output_directory = output_data_path + fr'\{specie_band_name}\uncertainties'
    # file_name_best_fit_spectrum = output_data_path + fr'\fits\{timing_index}.txt'
    # uncertainty_output_directory = output_data_path + fr'\uncertainties'
    output_files = {
        file_type: uncertainty_output_directory + fr'\{timing_index}.{file_type}' for file_type in ['metadata', 'png']
    }

    if os.path.exists(uncertainty_output_directory + fr'\{timing_index}.metadata'):
        return            # skip if already done
    elif os.path.exists(uncertainty_output_directory):
        pass
    else:
        os.makedirs(uncertainty_output_directory)

    # 1. define variables used for the calculations
    specie_name = get_specie(specie_band_name)

    best_fit_data = FitData.from_txt_file(file_name_best_fit_spectrum,
                                          thermal_equilibrium_spectrum=equilibrium_distribution)
    best_fit_data.read_metadata(specie=specie_name)
    best_fit_value = best_fit_data.parameters[get_parameter_name()].value      # [ppm]

    molecule_simulator = getattr(GetSimulator, f'get_{specie_name.upper()}_simulator')(
        get_database(specie_name), best_fit_data.x[0], best_fit_data.x[-1]
    )

    if equilibrium_distribution:
        molecule_distribution = GetDistibution.get_equilibrium_distributions(
            len(molecule_simulator._isotopologue_simulators)
        )
        spectrum_simulator = GetSpectrum.get_parameterised_equilibrium_fit_function(
            molecule_simulator, molecule_distribution, specie=specie_name
        )
    else:
        if specie_name == 'co2':
            molecule_distribution = GetDistibution.get_CO2_treanor_distributions()
        elif specie_name == 'co':
            molecule_distribution = GetDistibution.get_CO_treanor_distributions()
        else:
            raise ValueError(f"{specie_name=} has not (yet) a non-equilibrium spectrum")
        spectrum_simulator = GetSpectrum.get_parameterised_non_equilibrium_fit_function(
            molecule_simulator, molecule_distribution, specie=specie_name
        )

    def calculate_spectrum_function(x_values, updated_value):
        return spectrum_simulator(
            x_values,  best_fit_data.adjust_fraction(specie_name, **{get_parameter_name(): updated_value})
        )

    print(fr'uncertainty\{timing_index}')
    reduced_chi_squared_min, error_limits = uncertainty_analysis(
        best_fit_data.x, best_fit_data.transmittance, calculate_spectrum_function, best_fit_value,
        internal_error=internal_error, number_of_parameters=best_fit_data.number_of_variables,
        increment=resolution, verbose=True
    )

    def determine_uncertainty():
        if error_limits[0] == 0.0:
            return error_limits[1]
        return (error_limits[1] - error_limits[0]) / 2

    uncertainty_concentration = determine_uncertainty()
    if verbose:
        print(f'\t{uncertainty_concentration=}, from {error_limits=}')

    # 3. saving the results to a .metadata file
    metadata_results = configparser.ConfigParser()
    metadata_results.add_section('Variables')
    metadata_results.add_section('Files')

    metadata_results.set('Variables', 'number_of_variables', str(best_fit_data.number_of_variables))
    metadata_results.set('Variables', 'internal_error', str(internal_error))
    metadata_results.set('Variables', 'reduced_chi_squared_min', str(reduced_chi_squared_min))
    metadata_results.set('Variables', f'resolution_error_analysis', str(resolution))
    metadata_results.set('Variables', f'error_limit', str(error_limits))
    metadata_results.set('Variables', f'uncertainty', str(uncertainty_concentration))

    metadata_results.set('Files', 'best fit', file_name_best_fit_spectrum)
    metadata_results.set('Files', 'database', get_database(specie_name))

    with open(output_files['metadata'], 'w') as metadata_file:
        metadata_results.write(metadata_file)
    return
