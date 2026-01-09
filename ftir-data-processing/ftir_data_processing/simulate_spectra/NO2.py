import matplotlib.pyplot as plt
from context import *


hitran_database_file = r"D:\Promotie\data\database\no2_database.sqlite"
spectrum_data = spectrum.Spectrum.from_wavenumber_array_params(500, 3500, 0.04)

molecule_simulator = GetSimulator.get_NO2_simulator(
    hitran_database_file, spectrum_data.wavenumber_array.min(), spectrum_data.wavenumber_array.max()
)
molecule_simulator.set_instrumental_broadening_coefficient(0.2)

spectrum_data.data_array = calculate_spectrum.calculate_equilibrium_spectrum(
    spectrum_data.wavenumber_array,
    molecule_simulator, distributions.EquilibriumDistribution(),
    pressure_atm=1.,
    molecular_fraction=1e-3,
    temperature_gas_kelvin=300,
    path_length_cm=30
)

plt.plot(spectrum_data.wavenumber_array, spectrum_data.data_array)
plt.title('NO2')
plt.xlabel('wavenumber / cm$^{-1}$')
plt.ylabel('transmittance')
plt.grid(True)
plt.show()