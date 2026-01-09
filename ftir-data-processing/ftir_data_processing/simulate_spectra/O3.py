import matplotlib.pyplot as plt
from context import *

hitran_database_file = r"D:\Promotie\data\database\o3_database.sqlite"
spectrum_data = spectrum.Spectrum.from_wavenumber_array_params(900, 1900, 0.02)

molecule_simulator = GetSimulator.get_O3_simulator(
    hitran_database_file, spectrum_data.wavenumber_array.min(), spectrum_data.wavenumber_array.max()
)
molecule_simulator.set_instrumental_broadening_coefficient(0.15)

spectrum_data.data_array = calculate_spectrum.calculate_equilibrium_spectrum(
    spectrum_data.wavenumber_array,
    molecule_simulator, distributions.EquilibriumDistribution(),
    pressure_atm=1,
    molecular_fraction=2e-2,
    temperature_gas_kelvin=300,
    path_length_cm=1
)

plt.plot(spectrum_data.wavenumber_array, spectrum_data.data_array)
plt.title('O3')
plt.xlabel('wavenumber / cm$^{-1}$')
plt.ylabel('transmittance')
plt.grid(True)
plt.show()