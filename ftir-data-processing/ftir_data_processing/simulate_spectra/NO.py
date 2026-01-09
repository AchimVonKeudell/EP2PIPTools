import matplotlib.pyplot as plt
from context import *

hitran_database_file = r"D:\Promotie\data\database\no_database.sqlite"
spectrum_data = spectrum.Spectrum.from_wavenumber_array_params(1700, 2100, 0.001)

molecule_simulator = GetSimulator.get_NO_simulator(
    hitran_database_file, spectrum_data.wavenumber_array.min(), spectrum_data.wavenumber_array.max()
)
molecule_simulator.set_instrumental_broadening_coefficient(1.05)

spectrum_data.data_array = calculate_spectrum.calculate_equilibrium_spectrum(
    spectrum_data.wavenumber_array,
    molecule_simulator, distributions.EquilibriumDistribution(),
    pressure_atm=4.93e-3,
    molecular_fraction=1e-6,
    temperature_gas_kelvin=300,
    path_length_cm=42.6
)

plt.plot(spectrum_data.wavenumber_array, spectrum_data.data_array)
plt.title('NO')
plt.xlabel('wavenumber / cm$^{-1}$')
plt.ylabel('transmittance')
plt.grid(True)
plt.show()