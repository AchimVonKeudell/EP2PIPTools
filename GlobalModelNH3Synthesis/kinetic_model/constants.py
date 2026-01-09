from math import pi

# physical constants
avogadros_constant = 6.02214076e23                          # molecules / mol
boltzmann_constant = 1.380649e-23                           # [J / K]
electron_charge_coulomb = 1.60217663e-19                    # [C]
plancks_constant = 6.62607015e-34                           # [J * s]
mass_electron_kg = 9.109e-31

# conversion of units
constant_kJmol_to_kelvin = 1 / (1e-3 * avogadros_constant * boltzmann_constant)
constant_eV_to_kelvin = electron_charge_coulomb / boltzmann_constant
electron_volt_to_speed_squared = 2 * electron_charge_coulomb / mass_electron_kg

# other constants
NORM_FACTOR_MAXWELLIAN_EEDF = 2 / pi ** 0.5

# paths of the project
project_path = __file__[:__file__.find('kinetic_model')]
cross_sections_path = f"{project_path}\\Input\\LxCat"
electron_energy_distribution_directory = 'C:\\Users\\steij\\Documents\\Promotie\\scripts\\LoKI-master\\Code\\Output'
