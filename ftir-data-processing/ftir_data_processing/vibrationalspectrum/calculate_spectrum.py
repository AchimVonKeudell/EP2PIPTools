import numpy as np
from .spectrumsimulator.spectrum import Spectrum
from .spectrumsimulator import conditions, simulators, distributions, lineshapes, constants


def calculate_equilibrium_spectrum(
        wavenumbers,
        molecule_simulator: simulators.MoleculeSimulator,
        molecule_distribution: distributions.Distribution,
        *, pressure_atm, path_length_cm, molecular_fraction, temperature_gas_kelvin,
        simulation_transmission_threshold=0.9999
):
    """Calculates a transmittance spectrum assuming rotational-vibrational equilibrium that follows the given
    temperature_gas_kelvin.

    :param wavenumbers:
    :param molecule_simulator:
    :param molecule_distribution:
    :param pressure_atm:
    :param path_length_cm:
    :param molecular_fraction:
    :param temperature_gas_kelvin:
    :param simulation_transmission_threshold:
    :return:
    """
    spectrum_eq = Spectrum(wavenumbers, np.ones_like(wavenumbers))

    partial_condition = conditions.TransmissionMeasurementGasConditions(
        pressure_atm=pressure_atm,
        temperature_K=temperature_gas_kelvin,
        path_length_cm=path_length_cm
    ).get_species_conditions_from_mixing_ratio(molecular_fraction)

    # set rotational temperature:
    molecule_distribution.set_temperature_K(temperature_gas_kelvin)
    molecule_simulator.set_distribution(molecule_distribution)

    # calculating the spectrum
    spectrum_eq.data_array *= molecule_simulator.add_transmission_to_spectrum(
        spectrum=spectrum_eq, conditions=partial_condition,
        inclusion_strength_threshold=1e-25,
        simulation_transmission_threshold=simulation_transmission_threshold,
    )
    molecule_simulator.apply_instrumental_broadening(spectrum_eq)
    return spectrum_eq.data_array


def calculate_non_equilibrium_spectrum(wavenumbers, molecule_simulator, molecule_distributions,
                                       pressure_atm, path_length_cm, molecular_fraction, temperature_gas_kelvin,
                                       temperature_rot_kelvin, temperatures_vib_kelvin: tuple):
    """Calculate the transmittance of a non-thermal plasma with the provided temperatures and CO2 fraction with the
    addition of a slab of thermal CO2 gas.

    :param wavenumbers:
    :param molecule_simulator:
    :param molecule_distributions:
    :param pressure_atm:
    :param path_length_cm:
    :param molecular_fraction:
    :param temperature_gas_kelvin:
    :param temperature_rot_kelvin:
    :param temperatures_vib_kelvin:
    :return:
    """
    # Create a spectrum object to be able to simulate a spectrum
    spectrum_non_eq = Spectrum(wavenumbers, np.ones_like(wavenumbers))

    # Define the conditions of the measurement, for the partial gas of the given molecule
    partial_conditions = conditions.TransmissionMeasurementGasConditions(
        pressure_atm=pressure_atm,
        temperature_K=temperature_gas_kelvin,
        path_length_cm=path_length_cm
    ).get_species_conditions_from_mixing_ratio(molecular_fraction)

    # Set a Treanor distribution for the isotopologues
    for index in range(len(molecule_simulator._isotopologue_simulators)):
        molecule_distribution = molecule_distributions[index]
        molecule_distribution.set_temperature_K(temperature_rot_kelvin, *temperatures_vib_kelvin)
        molecule_simulator._isotopologue_simulators[index].set_distribution(molecule_distribution)

    # Add the transmission data to the spectrum object, and define the thresholds for line strength inclusion and
    # the threshold at which the simulation of the line will stop
    spectrum_non_eq.data_array *= molecule_simulator.add_transmission_to_spectrum(
        spectrum_non_eq, partial_conditions,
        inclusion_strength_threshold=1e-25, simulation_transmission_threshold=0.9999)
    molecule_simulator.apply_instrumental_broadening(spectrum_non_eq)
    return spectrum_non_eq.data_array
#
#
# def calculate_2_compartment_spectrum(molecule_simulator, non_plasma_distribution, plasma_distributions, wavenumbers,
#                                      pressure_atm, path_length_non_plasma_cm, path_length_plasma_cm, molar_fraction,
#                                      temperature_gas_kelvin, temperature_vib_kelvin):
#     spectrum_eq = Spectrum.from_wavenumber_array(wavenumbers)
#     spectrum_non_eq = Spectrum.from_wavenumber_array(wavenumbers)
#
#     partial_condition_1 = TransmissionMeasurementGasConditions(
#         pressure_atm=pressure_atm,
#         temperature_K=temperature_gas_kelvin,
#         path_length_cm=path_length_non_plasma_cm
#     ).get_species_conditions_from_mixing_ratio(molar_fraction)
#
#     partial_condition_2 = TransmissionMeasurementGasConditions(
#         pressure_atm=pressure_atm,
#         temperature_K=temperature_gas_kelvin,
#         path_length_cm=path_length_plasma_cm
#     ).get_species_conditions_from_mixing_ratio(molar_fraction)
#
#     # calculate the non-plasma compartment:
#     non_plasma_distribution.set_temperature_K(temperature_gas_kelvin)
#     molecule_simulator.set_distribution(non_plasma_distribution)
#
#     molecule_simulator.set_distribution(non_plasma_distribution)
#     molecule_simulator.add_transmission_to_spectrum(spectrum_eq, partial_condition_1,
#                                                     inclusion_strength_threshold=1e-25,
#                                                     simulation_transmission_threshold=0.9999)
#
#     # calculate the plasma (hot) compartment
#     for _plasma_distribution in plasma_distributions:
#         _plasma_distribution.set_temperature_K(temperature_gas_kelvin, temperature_vib_kelvin)
#
#     molecule_simulator.set_distribution(plasma_distributions)
#     molecule_simulator.add_transmission_to_spectrum(spectrum_non_eq, partial_condition_2,
#                                                     inclusion_strength_threshold=1e-25,
#                                                     simulation_transmission_threshold=0.9999)
#     total_spectrum = Spectrum.from_wavenumber_array(wavenumbers)
#     total_spectrum.data_array = spectrum_non_eq.data_array * spectrum_eq.data_array
#
#     molecule_simulator.apply_instrumental_broadening(total_spectrum)
#     return total_spectrum.data_array


def calculate_reflection_spectrum(
        wavenumbers,
        reference_molecular_spectra: dict, surface_group_bands: dict, *,
        pressure_atm: float, path_length_cm: float, temperature_gas_kelvin: float, molar_fractions: dict
):
    total_spectrum = Spectrum(wavenumbers, np.ones_like(wavenumbers))

    for molecule_label, reference_spectrum in reference_molecular_spectra.items():
        density_cm3 = (101325 * pressure_atm) / (constants.k_B * temperature_gas_kelvin) * 1e-6
        density_cm3 *= molar_fractions[molecule_label]
        concentration_ratio = (density_cm3 * path_length_cm) / (
                reference_spectrum.reference_number_density_cm3 * reference_spectrum.path_length_cm)

        total_spectrum.data_array *= reference_spectrum.interpolate(wavenumbers) ** concentration_ratio

    for surface_group_label, band_values in surface_group_bands.items():
        if band_values['line_shape'] == 'gaussian':
            band_line_shape = lineshapes.gaussian(
                wavenumbers,
                center=band_values['central_position_cm'],
                sigma=band_values['half_width_half_maximum_cm'] / np.sqrt(2 * np.log(2)),
                amplitude=band_values['amplitude']
            )
        elif band_values['line_shape'] == 'lorentzian':
            band_line_shape = lineshapes.lorentzian(
                wavenumbers,
                center=band_values['central_position_cm'],
                gamma=band_values['half_width_half_maximum_cm'],
                amplitude = band_values['amplitude']
            )
        else:
            raise ValueError(f"{band_values['line_shape']=} not known/implemented")

        total_spectrum.data_array *= np.exp(-band_line_shape)
    return total_spectrum.data_array