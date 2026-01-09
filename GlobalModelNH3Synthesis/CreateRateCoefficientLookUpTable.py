import os
import numpy as np

from kinetic_model.constants import (cross_sections_path, electron_energy_distribution_directory,
                                     electron_volt_to_speed_squared)
from kinetic_model.kinetic_reactions import *
from kinetic_model.kinetic_reactions.basic_reaction import Constant
from kinetic_model.kinetic_reactions.ElectronImpact.CrossSectionBasedMaxwellian import CrossSectionSet


"""
    create a file per reaction, with the following columns: E/N [Td], T_e [eV], rate_coefficient [cm3/s]
    In the header, note the 
        cross-section file
        location calculated EEDFs
"""


def calculate_reaction_rate_cm3_s(
        cross_section_data_m2: CrossSectionSet,
        electron_energy_density_function: tuple
):
    electron_speed_m_s = np.sqrt(electron_volt_to_speed_squared * cross_section_data_m2.energies_ev)
    eedf_values = np.interp(cross_section_data_m2.energies_ev, *electron_energy_density_function)
    eepf_values = np.sqrt(cross_section_data_m2.energies_ev) * eedf_values
    return 1e6 * np.sum(
        cross_section_data_m2.cross_sections_m2 * electron_speed_m_s * eepf_values * cross_section_data_m2.energy_steps_ev
    )


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    reaction_name = 'He_e_to_He23S_e'
    cross_section_file = f'{cross_sections_path}\\Helium\\{reaction_name}.dat'

    loki_subdirectory = 'He-rf-plasma'
    directory_eedfs = f"{electron_energy_distribution_directory}\\{loki_subdirectory}"
    swarm_parameters_file_format = f"{directory_eedfs}\\%s\\swarmParameters.txt"
    eedf_file_format = f"{directory_eedfs}\\%s\\eedf.txt"

    output_subdirectory = 'HeliumPlasma'
    # TODO how to save it, can we put more reactions together?

    reduced_electric_field_values = []
    electron_temperatures = []
    rate_coefficient_cm3s = []
    for reduced_e_field_dir in os.listdir(directory_eedfs):
        if '.txt' in reduced_e_field_dir:
            continue

        reduced_electric_field_values.append(
            float(reduced_e_field_dir.split("_")[1])
        )
        # print(f'E/N = {reduced_electric_field_values[-1]} Td')

        # get electron temperature
        parameters = {}
        with open(swarm_parameters_file_format % reduced_e_field_dir, 'r') as file_data:
            for line_data in file_data:
                key, item = line_data.split(' = ')
                try:
                    parameters[key.strip()] = Constant(*item.split())
                except ValueError as e:
                    pass
        electron_temperatures.append(parameters['Electron temperature'].value)

        # get EEDF and calculate the rate_coefficient k [cm3/s]
        eedf_data = np.loadtxt(eedf_file_format % reduced_e_field_dir, skiprows=1)
        rate_coefficient_cm3s.append(calculate_reaction_rate_cm3_s(
            cross_section_data_m2=CrossSectionSet.from_file(cross_section_file),
            electron_energy_density_function=(eedf_data[:, 0], eedf_data[:, 1])
        ))
        # TODO get cross section
        # TODO calculate rate_coefficient
    output = np.array([reduced_electric_field_values, electron_temperatures, rate_coefficient_cm3s]).T

    np.savetxt()


