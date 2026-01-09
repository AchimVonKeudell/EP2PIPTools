import numpy as np
from .spectrumsimulator.spectrum import Spectrum
from .spectrumsimulator import constants
from . import calculate_spectrum


def get_parameter_value(params, key, default=None):
    if key in params:
        return params[key].value
    elif default is None:
        raise ValueError(f'{key} not found in params')
    return default



def get_parameterized_co2_and_co_fit_function(molecule_simulator_co2, molecular_simulator_co,
                                              distributions_co2, distributions_co,
                                              distributions_co2_non_plasma=None, distributions_co_non_plasma=None):
    """
    Retrieve a fit function for a given molecule simulator that takes wavenumbers and a Parameters object as input to
    calculate absorption spectra.
    :param molecule_simulator_co2:
    :param molecular_simulator_co:
    :param distributions_co2:
    :param distributions_co:
    :param distributions_co2_non_plasma:
    :param distributions_co_non_plasma:
    :return:
    """
    if distributions_co2_non_plasma is None:
        distributions_co2_non_plasma = distributions_co2

    if distributions_co_non_plasma is None:
        distributions_co_non_plasma = distributions_co

    def get_fit_outcome(wavenumbers, params):
        # check if co2 and co should be fitted
        fit_co2 = False
        fit_co = False
        parameters = params.valuesdict()

        if "n_CO2" in parameters:
            fit_co2 = True
        if "n_CO" in parameters:
            fit_co = True

        # extract general properties
        pressure_atm = get_parameter_value(params, "pressure_atm")
        path_length_plasma_cm = get_parameter_value(params, "path_length_plasma_cm")
        path_length_non_plasma_cm = get_parameter_value(params, "path_length_non_plasma_cm", default=0)
        temperature_plasma_gas_K = get_parameter_value(params, "temperature_plasma_gas_K")
        temperature_non_plasma_gas_K = get_parameter_value(params, "temperature_non_plasma_gas_K", default=temperature_plasma_gas_K)
        wavenumber_shift = get_parameter_value(params, "wavenumber_shift", default=0)
        wavenumber_stretch_factor = 1 + get_parameter_value(params, "wavenumber_stretch_factor", default=0)

        # calculate total number density
        n_total = (pressure_atm * 101325) / (constants.k_B * temperature_plasma_gas_K)

        # stretch wavneumbers
        wavenumbers_center = np.average(wavenumbers)
        wavenumbers_stretched = wavenumber_stretch_factor * (wavenumbers - wavenumbers_center) + wavenumbers_center + wavenumber_shift

        # check if polynomial fitting parameters for the baseline are present
        if "p0" in params.valuesdict():
            poly_coefficients = [params["p0"].value + 1]
            poly_order = 0

            while "p%i" % (poly_order+1) in params.valuesdict():
                poly_coefficients.insert(0, params["p%i" % (poly_order+1)].value)
                poly_order += 1

            baseline = np.polyval(poly_coefficients, wavenumbers_stretched - wavenumbers_center)
        else:
            baseline = np.ones_like(wavenumbers_stretched)

        # create spectrum object to hold output
        spectrum_total = Spectrum(wavenumbers_stretched, baseline)

        # if co2 is enabled
        if fit_co2:
            n_CO2 = get_parameter_value(params, "n_CO2") * 1e20
            fraction_CO2 = n_CO2 / n_total

            temperature_plasma_rot_K = get_parameter_value(params, "temperature_plasma_rot_K", default=temperature_plasma_gas_K)
            temperature_plasma_vib_12_K = get_parameter_value(params, "temperature_plasma_vib_12_K", default=temperature_plasma_rot_K)
            temperature_plasma_vib_3_K = get_parameter_value(params, "temperature_plasma_vib_3_K", default=temperature_plasma_rot_K)
            fit_outcome_co2 = calculate_spectrum.calculate_non_equilibrium_spectrum(
                wavenumbers_stretched, molecule_simulator_co2, distributions_co2,
                pressure_atm, path_length_plasma_cm, fraction_CO2, temperature_plasma_gas_K,
                temperature_plasma_rot_K, (temperature_plasma_vib_12_K, temperature_plasma_vib_3_K))

            if path_length_non_plasma_cm > 0:
                temperature_non_plasma_rot_K = get_parameter_value(params, "temperature_non_plasma_rot_K", default=temperature_non_plasma_gas_K)
                temperature_non_plasma_vib_12_K = get_parameter_value(params, "temperature_non_plasma_vib_12_K", default=temperature_non_plasma_rot_K)
                temperature_non_plasma_vib_3_K = get_parameter_value(params, "temperature_non_plasma_vib_3_K", default=temperature_non_plasma_rot_K)
                fit_outcome_co2 *= calculate_spectrum.calculate_non_equilibrium_spectrum(
                    wavenumbers_stretched, molecule_simulator_co2, distributions_co2_non_plasma,
                    pressure_atm, path_length_non_plasma_cm, fraction_CO2, temperature_non_plasma_gas_K,
                    temperature_non_plasma_rot_K, (temperature_non_plasma_vib_12_K, temperature_non_plasma_vib_3_K))

            spectrum_total.data_array *= fit_outcome_co2

        # if co is enabled
        if fit_co:
            n_CO = get_parameter_value(params, "n_CO") * 1e20
            fraction_CO = n_CO / n_total

            temperature_plasma_rot_K = get_parameter_value(params, "temperature_plasma_rot_K", default=temperature_plasma_gas_K)
            temperature_plasma_vib_CO_K = get_parameter_value(params, "temperature_plasma_vib_CO_K", default=temperature_plasma_rot_K)
            fit_outcome_co = calculate_spectrum.calculate_non_equilibrium_spectrum(
                wavenumbers_stretched, molecular_simulator_co, distributions_co, pressure_atm, path_length_plasma_cm,
                fraction_CO, temperature_plasma_gas_K, temperature_plasma_rot_K, (temperature_plasma_vib_CO_K,))

            if path_length_non_plasma_cm > 0:
                temperature_non_plasma_rot_K = get_parameter_value(params, "temperature_non_plasma_rot_K", default=temperature_non_plasma_gas_K)
                temperature_non_plasma_vib_CO_K = get_parameter_value(params, "temperature_non_plasma_vib_CO_K", default=temperature_non_plasma_rot_K)
                fit_outcome_co *= calculate_spectrum.calculate_non_equilibrium_spectrum(
                    wavenumbers_stretched, molecular_simulator_co, distributions_co_non_plasma, pressure_atm,
                    path_length_non_plasma_cm, fraction_CO, temperature_non_plasma_gas_K, temperature_non_plasma_rot_K,
                    (temperature_non_plasma_vib_CO_K,))

            spectrum_total.data_array *= fit_outcome_co
        return spectrum_total.data_array
    return get_fit_outcome


def get_parameterised_equilibrium_fit_function(molecule_simulator, molecule_distribution, *, specie: str):
    """Retrieve a function to fit an IR spectrum in thermal equilibrium.

    :param molecule_simulator:
    :param molecule_distribution:
    :param specie:
    :return:
    """

    if specie.lower() not in ['co2', 'co', 'nh3', 'n2o', 'no', 'no2', 'o3', 'ch4']:
        raise NotImplementedError(f"{specie.lower()=} not implemented")

    def get_fit_outcome(wavenumbers, params):
        """Function that calculates the spectrum using the given wavenumbers and parameters (in params).

        The wavenumber_offset is an unrealistic parameter, but very useful. An apparent shift might occur in FTIR spectra
        when the optical path difference is miscalculated. This happens when the IR beam enters under and angle. Or
        the beam is diverging too much, which introduces a convolution of various spectra that each have a different
        wavenumber scale.

        :param wavenumbers:
        :param params:
        :return:
        """

        # Correcting the wavenumber-axis
        wavenumber_offset = get_parameter_value(params, 'wavenumber_offset', default=0)
        wavenumber_stretch = get_parameter_value(params, 'wavenumber_stretch', default=1)
        wavenumbers_corrected = wavenumber_stretch * (wavenumbers - wavenumber_offset)

        # Correcting the baseline
        if "p0" in params.valuesdict():
            poly_coefficients = [params["p0"].value + 1]
            poly_order = 0

            while "p%i" % (poly_order+1) in params.valuesdict():
                poly_coefficients.insert(0, params["p%i" % (poly_order+1)].value)
                poly_order += 1

            baseline = np.polyval(poly_coefficients, wavenumbers_corrected)
        else:
            baseline = np.ones_like(wavenumbers_corrected)

        spectrum_total = baseline

        # Acquiring the other parameters
        pressure_atm = get_parameter_value(params, 'pressure_atm')
        path_length_cm = get_parameter_value(params, 'path_length_cm')
        rotational_temperature_kelvin = get_parameter_value(params, 'temperature_rotation_kelvin')
        specie_fraction = get_parameter_value(params, f'fraction_{specie.lower()}')

        # Define the instrumental broadening
        molecule_simulator.set_instrumental_broadening_coefficient(
            sigma=get_parameter_value(params, 'instrumental_broadening_sigma', default=0)
        )

        spectrum_total *= calculate_spectrum.calculate_equilibrium_spectrum(
            wavenumbers_corrected, molecule_simulator, molecule_distribution,
            pressure_atm=pressure_atm, path_length_cm=path_length_cm, molecular_fraction=specie_fraction,
            temperature_gas_kelvin=rotational_temperature_kelvin
        )
        return spectrum_total
    return get_fit_outcome


def get_parameterised_non_equilibrium_fit_function(molecule_simulator, molecule_distributions, *, specie: str):
    """Retrieve a function to fit an IR spectrum in thermal non-equilibrium.

    :param molecule_simulator:
    :param molecule_distributions:
    :param specie:
    :return:
    """

    if specie.lower() not in ['co2', 'co', 'nh3', 'n2o']:
        raise ValueError(f"{specie.lower()=} not known")

    def get_fit_outcome(wavenumbers, params):
        """Function that calculates the spectrum using the given wavenumbers and parameters (in params)

        :param wavenumbers:
        :param params:
        :return:
        """

        # correcting the wavenumber-axis
        wavenumber_offset = get_parameter_value(params, 'wavenumber_offset')
        wavenumbers_corrected = wavenumbers - wavenumber_offset

        # Correcting the baseline
        if "p0" in params.valuesdict():
            poly_coefficients = [params["p0"].value + 1]
            poly_order = 0

            while "p%i" % (poly_order+1) in params.valuesdict():
                poly_coefficients.insert(0, params["p%i" % (poly_order+1)].value)
                poly_order += 1

            baseline = np.polyval(poly_coefficients, wavenumbers_corrected)
        else:
            baseline = np.ones_like(wavenumbers_corrected)

        spectrum_total = baseline

        # acquiring the other parameters
        pressure_atm = get_parameter_value(params, 'pressure_atm')
        path_length_cm = get_parameter_value(params, 'path_length_cm')
        specie_fraction = get_parameter_value(params, f'fraction_{specie.lower()}')
        rotational_temperature = get_parameter_value(params, 'temperature_rotation_kelvin')
        if specie.lower() == 'co2':
            vibrational_temperature_12 = get_parameter_value(params, "temperature_vibration_12_kelvin")
            vibrational_temperature_3 = get_parameter_value(params, "temperature_vibration_3_kelvin")
            vibrational_temperatures = vibrational_temperature_12, vibrational_temperature_3
        else:
            vibrational_temperature = get_parameter_value(params, 'temperature_vibration_kelvin')
            vibrational_temperatures = vibrational_temperature,

        # define the instrumental broadening
        molecule_simulator.set_instrumental_broadening_coefficient(
            sigma=get_parameter_value(params, 'instrumental_broadening_sigma', default=0)
        )
        spectrum_total *= calculate_spectrum.calculate_non_equilibrium_spectrum(
            wavenumbers_corrected, molecule_simulator, molecule_distributions,
            pressure_atm=pressure_atm, path_length_cm=path_length_cm, molecular_fraction=specie_fraction,
            temperature_gas_kelvin=rotational_temperature, temperature_rot_kelvin=rotational_temperature,
            temperatures_vib_kelvin=vibrational_temperatures)
        return spectrum_total

    return get_fit_outcome
#
#
# def get_parameterised_2_compartment_fit_function(molecule_simulator, distribution_non_plasma, distribution_plasma,
#                                                  *, specie: str):
#     if specie not in ['co']:
#         raise ValueError(f'{specie=} not implemented')
#
#     def get_fit_outcome(wavenumbers, params):
#         """Function that calculates the spectrum using the given wavenumbers and parameters (in params)
#
#         :param wavenumbers:
#         :param params:
#         :return:
#         """
#
#         # Correcting the wavenumber-axis
#         wavenumber_offset = get_parameter_value(params, 'wavenumber_offset')
#         wavenumbers_corrected = wavenumbers - wavenumber_offset
#
#         # Correcting the baseline
#         if "p0" in params.valuesdict():
#             poly_coefficients = [params["p0"].value + 1]
#             poly_order = 0
#
#             while "p%i" % (poly_order+1) in params.valuesdict():
#                 poly_coefficients.insert(0, params["p%i" % (poly_order+1)].value)
#                 poly_order += 1
#
#             baseline = np.polyval(poly_coefficients, wavenumbers_corrected)
#         else:
#             baseline = np.ones_like(wavenumbers_corrected)
#
#         spectrum_total = baseline
#
#         # Acquiring the other parameters
#         pressure_atm = get_parameter_value(params, 'pressure_atm')
#         path_length_non_plasma_cm = get_parameter_value(params, 'path_length_non_plasma_cm')
#         path_length_plasma_cm = get_parameter_value(params, 'path_length_plasma_cm')
#
#         rotational_temperature_kelvin = get_parameter_value(params, 'temperature_rotation_kelvin')
#
#         delta_vibration_rotation = get_parameter_value(params, 'delta_temperature_vibration_rotation_kelvin')
#         vibrational_temperature_kelvin = rotational_temperature_kelvin + delta_vibration_rotation
#
#         specie_fraction = get_parameter_value(params, f'fraction_{specie.lower()}')
#
#         # Define the instrumental broadening
#         molecule_simulator.set_instrumental_broadening_coefficient(
#             sigma=get_parameter_value(params, 'instrumental_broadening_sigma', default=0)
#         )
#
#         spectrum_total *= calculate_2_compartment_spectrum(
#             molecule_simulator, distribution_non_plasma, distribution_plasma, wavenumbers_corrected,
#             pressure_atm, path_length_non_plasma_cm, path_length_plasma_cm, specie_fraction,
#             temperature_gas_kelvin=rotational_temperature_kelvin,
#             temperature_vib_kelvin=vibrational_temperature_kelvin,
#         )
#         return spectrum_total
#     return get_fit_outcome


def get_parameterised_reflection_spectrum(reference_molecular_spectra: dict, band_line_shape: dict):
    def get_fit_outcome(wavenumbers, params):
        """Function that calculates the spectrum using the given wavenumbers and parameters (in params)

        :param wavenumbers:
        :param params:
        :return:
        """

        # Correcting the wavenumber-axis
        wavenumber_offset = get_parameter_value(params, 'wavenumber_offset', default=0)
        wavenumber_stretch = get_parameter_value(params, 'wavenumber_stretch', default=1)
        wavenumbers_corrected = wavenumber_stretch * (wavenumbers - wavenumber_offset)

        # Correcting the baseline
        if "p0" in params.valuesdict():
            poly_coefficients = [params["p0"].value + 1]
            poly_order = 0

            while "p%i" % (poly_order+1) in params.valuesdict():
                poly_coefficients.insert(0, params["p%i" % (poly_order+1)].value)
                poly_order += 1

            baseline = np.polyval(poly_coefficients, wavenumbers_corrected)
        else:
            baseline = np.ones_like(wavenumbers_corrected)

        spectrum_total = baseline

        # Acquiring the other parameters
        pressure_atm = get_parameter_value(params, 'pressure_atm')
        path_length_cm = get_parameter_value(params, 'path_length_cm')
        rotational_temperature_kelvin = get_parameter_value(params, 'temperature_rotation_kelvin')

        molar_fractions = {}
        for molecule_label in reference_molecular_spectra.keys():
            molar_fractions[molecule_label] = get_parameter_value(params, f'fraction_{molecule_label.lower()}')

        surface_group_bands = {}
        for band_label, line_shape in band_line_shape.items():
            surface_group_bands[band_label] = {
                'line_shape': line_shape,
                'central_position_cm': get_parameter_value(params, f'central_position_{band_label}_cm'),
                'half_width_half_maximum_cm': get_parameter_value(
                    params, f'half_width_half_maximum_{band_label}_cm'),
                'amplitude': get_parameter_value(params, f'amplitude_{band_label}'),
            }

        spectrum_total *= calculate_spectrum.calculate_reflection_spectrum(
            wavenumbers_corrected, reference_molecular_spectra, surface_group_bands,
            pressure_atm=pressure_atm, path_length_cm=path_length_cm,
            temperature_gas_kelvin=rotational_temperature_kelvin, molar_fractions=molar_fractions,
        )
        return spectrum_total
    return get_fit_outcome