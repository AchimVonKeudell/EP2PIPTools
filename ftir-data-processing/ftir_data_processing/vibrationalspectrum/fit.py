import os.path

import numpy as np
from lmfit import Parameters, minimize, report_fit
from .spectrumsimulator.spectrum import Spectrum
from . import GetParameterisedSpectrum, GetSimulator, GetDistibution, plot, metadata


def get_spectrum(full_file_path: str, **kwargs) -> Spectrum:
    """Acquire the spectrum, from the given path

    :param full_file_path:
    :param kwargs: optionally giving
        `resolution`: to change the wavenumber spacing
        `other_contributions`: removing these other contributions to isolate specific effects.
    :return:
    """

    loaded_data = np.loadtxt(full_file_path)

    if 'other_contributions' in kwargs:
        current_processing_step = os.path.dirname(full_file_path)
        current_spectrum_fn = os.path.basename(full_file_path)
        for other_contribution in kwargs['other_contributions']:
            other_contribution_fn = f'{os.path.dirname(current_processing_step)}\\{other_contribution}\\{current_spectrum_fn}'
            print(f'{other_contribution_fn=}')
            other_contribution_data = np.loadtxt(other_contribution_fn)
            loaded_data[:, 1] /= other_contribution_data[:, 2]


    if 'wavenumber_resolution' in kwargs:
        resolution = kwargs['wavenumber_resolution']
        x_values = np.arange(loaded_data[:, 0].min(), loaded_data[:, 0].max() + resolution, resolution)
        return Spectrum(x_values, np.interp(x_values, loaded_data[:, 0], loaded_data[:, 1]))
    return Spectrum(loaded_data[:, 0], loaded_data[:, 1])


def add_parameter(parameters: Parameters, dictionary: dict, label: str, key: str = None):
    """Add parameters to the Parameters object from the lmft routine.
    It looks in `dictionary` for `label`. If it's not there, then it tries to search for `label` + `_range` which
    indicates the range over which the parameter may be varied.

    :param parameters:
    :param dictionary:
    :param label:
    :param key:
    :return:
    """
    if key is None:
        key = label

    if key in dictionary:
        value = dictionary[key]
        parameters.add(label, value=float(value), vary=False)
        return

    key_range = key + '_range'
    if key_range in dictionary:
        values = [float(value) for value in dictionary[key_range]]
        parameters.add(label, value=np.mean(values), min=values[0], max=values[1])
        return

    raise ValueError(f"ISSUE with CONFIG file\n"
                     f"{label}: {dictionary=}\n should have {key=} or {key_range=}")


def obtain_best_fit(fit_function, measured_spectrum: Spectrum, parameters, **kwargs):
    def minimise_fit_function(params, x, data):
        if 'fit_derivative' in kwargs:
            if kwargs['fit_derivative']:
                return np.diff(fit_function(x, params)) - np.diff(data)
        return fit_function(x, params) - data

    output = minimize(minimise_fit_function, parameters,
                       args=(measured_spectrum.wavenumber_array, measured_spectrum.data_array))
    report_fit(output)
    return output


def save_results(output_files: iter,
                 files_used: dict,
                 best_fit_parameters: Parameters,
                 measured_spectrum: Spectrum, best_fit_spectrum: np.ndarray,
                 **kwargs):
    """

    :param output_files: a list the file name ending the extensions following file_extension.EXTENSIONS
    :param files_used:
    :param best_fit_parameters:
    :param measured_spectrum:
    :param best_fit_spectrum:
    :return;
    """

    # saving the measured and fitted spectra in .csv and .png formats
    np.savetxt(output_files[0],
               np.transpose((measured_spectrum.wavenumber_array, measured_spectrum.data_array, best_fit_spectrum)),
               delimiter='\t', header="wavenumber/cm^-1\ttransmittance,fit")

    # saving the metadata
    config = metadata.Metadata(file_path=output_files[1])
    with config as mtd:
        mtd.set_best_fit_parameters(best_fit_parameters)
        mtd.set_files(**files_used)

    # saving the figure
    plot.plot_spectra(output_files[2], measured_spectrum.wavenumber_array,
                      measured_spectrum.data_array, best_fit_spectrum, **kwargs)
    return


def fit_equilibrium_spectrum(spectrum_file, output_files: list, molecule_name: str, hitran_database_file, *,
                             pressure_atm, path_length_cm, **kwargs):
    """Fits the IR spectrum of a single specie (defined with 'molecule_name').

    It saves the resulting parameters and spectra in 3 different files.
        output_files:
            #1 will store the resulting fit parameters, e.g. fraction of NH3;
            #2 will store the measured and fit spectrum values;
            #3 will store the measured and fit spectrum in a figure.

    :param spectrum_file:
    :param output_files:            a list the file name ending the extensions following file_extension.EXTENSIONS
    :param molecule_name:           co2, co, nh3, n2o, ...
    :param hitran_database_file:    HITRAN SQlite database
    :param pressure_atm:
    :param path_length_cm:
    :return:
    """

    molecule_name = molecule_name.lower()       # making sure that no upper letters are present

    # 1) obtaining the transmittance spectrum
    experimental_data = get_spectrum(spectrum_file, **kwargs)

    # 2) acquiring the simulation function
    molecule_simulator = getattr(GetSimulator, f'get_{molecule_name.upper()}_simulator')(
        hitran_database_file, experimental_data.wavenumber_array.min(), experimental_data.wavenumber_array.max()
    )

    fit_function = GetParameterisedSpectrum.get_parameterised_equilibrium_fit_function(
        molecule_simulator,
        GetDistibution.get_equilibrium_distributions(),
        specie=molecule_name
    )

    # 3) define the parameters
    parameters = Parameters()
    parameters.add("pressure_atm", value=pressure_atm, vary=False)
    parameters.add("path_length_cm", value=path_length_cm, vary=False)

    add_parameter(parameters, kwargs, f"fraction_{molecule_name}")
    add_parameter(parameters, kwargs, 'instrumental_broadening_sigma')

    try:
        add_parameter(parameters, kwargs, 'wavenumber_offset')
    except ValueError:
        add_parameter(parameters, kwargs, 'wavenumber_stretch')


    if 'p0' in kwargs:
        parameters.add('p0', value=kwargs['p0'], vary=False)
        order = 1
        while f'p{order}' in kwargs:
            parameters.add(f'p{order}', value=kwargs[f'p{order}'], vary=False)
            order += 1
    else:
        parameters.add('p0', value=0, vary=False)

    try:
        add_parameter(parameters, kwargs, 'temperature_rotation_kelvin', 'gas_temperature_kelvin')
    except ValueError:
        add_parameter(parameters, kwargs, 'temperature_rotation_kelvin')

    # 4) fitting the spectrum
    results = obtain_best_fit(
        fit_function=fit_function,
        measured_spectrum=experimental_data,
        parameters=parameters,
        **kwargs
    )

    # 5) saving the fit results
    save_results(
        output_files=output_files,
        files_used={'input_file_name': spectrum_file, 'sqlite_data_file_name': hitran_database_file},
        best_fit_parameters=results.params,
        measured_spectrum=experimental_data,
        best_fit_spectrum=fit_function(experimental_data.wavenumber_array, results.params),
        **kwargs
    )
    return


def fit_non_equilibrium_spectrum(spectrum_file, output_files: iter, molecule_name: str, hitran_database_file, *,
                                 pressure_atm: float, path_length_cm: float, **kwargs):
    """

    :param spectrum_file:
    :param output_files:    a list the file name ending the extensions following file_extension.EXTENSIONS
    :param molecule_name:   co2, co, ...
    :param hitran_database_file:   HITRAN SQlite database
    :param pressure_atm:
    :param path_length_cm:
    :param kwargs:
    :return:
    """

    molecule_name = molecule_name.lower()       # making sure no uppercase letters are present

    # 1) obtaining the transmittance spectrum
    experimental_data = get_spectrum(spectrum_file, **kwargs)

    # 2a) acquiring the simulation function
    molecule_simulator = getattr(GetSimulator, f'get_{molecule_name.upper()}_simulator')(
        hitran_database_file, experimental_data.wavenumber_array[0], max(experimental_data.wavenumber_array))

    if molecule_name == 'co2':
        molecule_distribution = GetDistibution.get_CO2_treanor_distributions()
    elif molecule_name == 'co':
        molecule_distribution = GetDistibution.get_CO_treanor_distributions()
    else:
        raise ValueError(f"no distribution object found for {molecule_name=}")

    fit_function = GetParameterisedSpectrum.get_parameterised_non_equilibrium_fit_function(
        molecule_simulator,
        molecule_distribution,
        specie=molecule_name
    )

    # 2b) define the parameters
    parameters = Parameters()
    parameters.add("pressure_atm", value=pressure_atm, vary=False)
    parameters.add("path_length_cm", value=path_length_cm, vary=False)

    add_parameter(parameters, kwargs, f"fraction_{molecule_name}")
    add_parameter(parameters, kwargs,'wavenumber_offset')
    add_parameter(parameters, kwargs, 'instrumental_broadening_sigma')

    if 'p0' in kwargs:
        parameters.add('p0', value=kwargs['p0'], vary=False)
        order = 1
        while f'p{order}' in kwargs:
            parameters.add(f'p{order}', value=kwargs[f'p{order}'], vary=False)
            order += 1
    else:
        parameters.add('p0', value=0, vary=False)

    add_parameter(parameters, kwargs, 'temperature_rotation_kelvin')

    if molecule_name in ['co', 'no', 'n2', 'o2']:     # for molecules with 1 vibrational temperature
        add_parameter(parameters, kwargs, 'temperature_vibration_kelvin')
    elif molecule_name == 'co2':
        add_parameter(parameters, kwargs, 'temperature_vibration_12_kelvin')
        add_parameter(parameters, kwargs, 'temperature_vibration_3_kelvin')
    else:
        raise ValueError(f'{molecule_name=} not implemented when configuring the vibrational temperatures')

    # 2c) fitting the spectrum
    results = obtain_best_fit(
        fit_function=fit_function,
        measured_spectrum=experimental_data,
        parameters=parameters,
        **kwargs
    )

    # 3) saving the fit results
    save_results(
        output_files=output_files,
        files_used={'input_file_name': spectrum_file, 'sqlite_data_file_name': hitran_database_file},
        best_fit_parameters=results.params,
        measured_spectrum=experimental_data,
        best_fit_spectrum=fit_function(experimental_data.wavenumber_array, results.params)
    )
    return

# def fit_2compartment_equilibrium_spectrum(spectrum_file, output_files: iter, database_file, *,
#                                           pressure_atm, path_length_non_plasma_cm: float,
#                                           path_length_plasma_cm: float, molecule_name: str, **kwargs):
#     """Fits the IR spectrum of a single specie (defined with 'molecule_name').
#
#     It saves the resulting parameters and spectra in 3 different files.
#         output_files:
#             #1 will store the resulting fit parameters, e.g. fraction of NH3;
#             #2 will store the measured and fit spectrum values;
#             #3 will store the measured and fit spectrum in a figure.
#
#     :param spectrum_file:
#     :param output_files:                        3 files: .METADATA, .txt, .png
#     :param database_file:                       HITRAN SQlite database
#     :param pressure_atm:
#     :param path_length_non_plasma_cm:
#     :param path_length_plasma_cm:
#     :param molecule_name:                       co2, co, nh3, n2o
#     :return:
#     """
#
#     molecule_name = molecule_name.lower()  # making sure that no upper letters are present
#
#     # obtaining the transmittance spectrum
#     resolution = None
#     if 'resolution' in kwargs:
#         resolution = kwargs['resolution']
#     experimental_data = get_spectrum(spectrum_file, resolution)
#     # acquiring the simulation function
#     molecule_simulator = getattr(GetSimulator, f'get_{molecule_name.upper()}_simulator')(
#         database_file, experimental_data.wavenumber_array[0], max(experimental_data.wavenumber_array))
#
#     non_plasma_distribution = GetDistibution.get_equilibrium_distributions()
#
#     if molecule_name == 'co':
#         plasma_distribution = GetDistibution.get_CO_treanor_distributions()
#     else:
#         NotImplementedError(f'{molecule_name=} not implemented')
#
#     fit_function = GetSpectrum.get_parameterised_2_compartment_fit_function(
#         molecule_simulator, non_plasma_distribution, plasma_distribution, specie=molecule_name
#     )
#
#     # define the parameters for compartment 1
#     parameters = Parameters()
#     parameters.add("pressure_atm", value=pressure_atm, vary=False)
#     parameters.add("path_length_non_plasma_cm", value=path_length_non_plasma_cm, vary=False)
#     parameters.add("path_length_plasma_cm", value=path_length_plasma_cm, vary=False)
#
#     if 'wavenumber_offset' in kwargs:
#         parameters.add('wavenumber_offset', value=kwargs['wavenumber_offset'], vary=False)
#     else:
#         parameters.add('wavenumber_offset', value=0, min=-1, max=1)
#
#     if 'p0' in kwargs:
#         parameters.add('p0', value=kwargs['p0'], vary=False)
#         order = 1
#         while f'p{order}' in kwargs:
#             parameters.add(f'p{order}', value=kwargs[f'p{order}'], vary=False)
#             order += 1
#     else:
#         parameters.add('p0', value=0, vary=False)
#         # parameters.add('p1', value=0, min=-1e-5, max=1e-5)
#         # parameters.add('p2', value=0, min=-1e-5, max=1e-5)
#
#     if 'instrumental_broadening_sigma' in kwargs:
#         parameters.add("instrumental_broadening_sigma", value=kwargs['instrumental_broadening_sigma'], vary=False)
#     else:
#         parameters.add("instrumental_broadening_sigma", value=0.5, min=0., max=5.0, vary=True)
#
#     parameters.add(f"fraction_{molecule_name}", value=5e-6, min=0, max=0.1)
#
#     if 'temperature_rotation_kelvin' in kwargs:
#         parameters.add('temperature_rotation_kelvin', value=kwargs['temperature_rotation_kelvin'], vary=False)
#     else:
#         parameters.add('temperature_rotation_kelvin', value=300, min=200, max=3000)
#
#     parameters.add('delta_temperature_vibration_rotation_kelvin', value=500, min=0, max=3000)
#
#     # fitting the spectra
#     results = minimize(lambda params, x, data: fit_function(x, params) - data, parameters,
#                        args=(experimental_data.wavenumber_array, experimental_data.data_array))
#     report_fit(results)
#
#     # saving the fit results
#     metadata = ConfigParser()
#     metadata.add_section("Parameters")
#     metadata.add_section("Errors")
#     metadata.add_section("Fixed")
#     for key, parameter in results.params.items():
#         metadata.set("Parameters", key, str(parameter.value))
#         metadata.set("Errors", key, str(parameter.stderr))
#         metadata.set("Fixed", key, 'no' if parameter.vary else 'yes')
#     metadata.add_section("Files")
#     metadata.set("Files", 'spectrum', spectrum_file)
#     metadata.set("Files", 'database', database_file)
#     metadata.write(open(output_files[0], 'w'))
#
#     # saving the measured and fitted spectra in .csv and .png formats
#     fitted_spectrum = fit_function(experimental_data.wavenumber_array, results.params)
#     np.savetxt(output_files[1],
#                np.transpose((experimental_data.wavenumber_array, experimental_data.data_array, fitted_spectrum)),
#                delimiter=',', header="wavenumber [cm^-1],transmittance,fit")
#
#     plot.plot_spectra(output_files[2], experimental_data.wavenumber_array,
#                       experimental_data.data_array, fitted_spectrum)
#     return


def fit_reflection_spectrum(spectrum_file, output_files: iter, reference_molecular_spectra: dict, *,
                            pressure_atm: float, path_length_cm: float, gas_temperature_kelvin: float,
                            contributions_gas_phase: dict, contributions_surface: dict,
                            baseline_wavenumber_range, baseline_polynomial_order, **kwargs):
    """Fitting a gas phase specie + surface bands (approximated by a Gaussian or Lorentzian)

    :param spectrum_file:
    :param output_files: ending the extensions following file_extension.EXTENSIONS
    :param reference_molecular_spectra:
    :param pressure_atm:
    :param path_length_cm:
    :param gas_temperature_kelvin:
    :param contributions_gas_phase:
    :param contributions_surface:
    :param baseline_wavenumber_range:
    :param baseline_polynomial_order:
    :param kwargs:
    :return:
    """

    # 1) obtaining the transmittance spectrum
    experimental_data = get_spectrum(spectrum_file, **kwargs)

    if baseline_wavenumber_range is None:
        baseline_coefficients = [1.0]
    else:
        baseline_experimental_data = experimental_data.get_snip(baseline_wavenumber_range)
        baseline_coefficients = np.polyfit(baseline_experimental_data.wavenumber_array,
                                           baseline_experimental_data.data_array, baseline_polynomial_order)
    baseline_spectrum = np.polyval(baseline_coefficients, experimental_data.wavenumber_array)

    if 'fit_wavenumber_range' in kwargs:
        fit_experimental_data = experimental_data.get_snip(kwargs['fit_wavenumber_range'])
    else:
        fit_experimental_data = experimental_data

    fit_experimental_data.data_array /= np.polyval(baseline_coefficients, fit_experimental_data.wavenumber_array)

    # 2) acquiring the simulation function
    parameters = Parameters()
    parameters.add("pressure_atm", value=pressure_atm, vary=False)
    parameters.add("path_length_cm", value=path_length_cm, vary=False)
    parameters.add('temperature_rotation_kelvin', value=gas_temperature_kelvin, vary=False)

    if 'wavenumber_offset' in kwargs:
        parameters.add('wavenumber_offset', value=kwargs['wavenumber_offset'], vary=False)
    else:
        parameters.add('wavenumber_offset', value=0, min=-1, max=1)

    if 'wavenumber_stretch_range' in kwargs:
        add_parameter(parameters, kwargs, 'wavenumber_stretch')

    ## 2a. get the information for the gaseous molecules
    files_used_dict = {'input_file_name': spectrum_file}
    for molecule_label, contribution_gas_phase in contributions_gas_phase.items():
        files_used_dict[f'reference_spectrum_{molecule_label}'] = reference_molecular_spectra[molecule_label].file_name

        add_parameter(parameters, contribution_gas_phase,f'fraction_{molecule_label.lower()}', 'molar_fraction')

    ## 2b. acquire the surface absorption bands
    band_line_shapes = {}
    for band_label, band_parameters in contributions_surface.items():
        band_line_shapes[band_label] = band_parameters['line_shape']

        add_parameter(parameters, band_parameters,f'central_position_{band_label}_cm','central_wavenumber_cm')
        add_parameter(parameters, band_parameters,f'half_width_half_maximum_{band_label}_cm','half_width_half_maximum_cm')
        add_parameter(parameters, band_parameters,f'amplitude_{band_label}','amplitude')


    fit_function = GetParameterisedSpectrum.get_parameterised_reflection_spectrum(
        reference_molecular_spectra, band_line_shapes
    )

    # 4) fitting the spectrum
    results = obtain_best_fit(
        fit_function=fit_function,
        measured_spectrum=fit_experimental_data,
        parameters=parameters,
        **kwargs
    )

    # 5) saving the fit results
    save_results(
        output_files=output_files,
        files_used=files_used_dict,
        best_fit_parameters=results.params,
        measured_spectrum=experimental_data,
        best_fit_spectrum=fit_function(experimental_data.wavenumber_array, results.params),
        zoom=1e3, baseline_spectrum=baseline_spectrum,
        **kwargs
    )
    return

