import matplotlib.pyplot as plt
from context import *

hitran_database_file = r"D:\Promotie\data\database\n2o_database.sqlite"
spectrum_data = spectrum.Spectrum.from_wavenumber_array_params(2000, 2400, 0.001)

molecule_simulator = GetSimulator.get_N2O_simulator(
    hitran_database_file, spectrum_data.wavenumber_array.min(), spectrum_data.wavenumber_array.max()
)
molecule_simulator.set_instrumental_broadening_coefficient(1.05)

spectrum_data.data_array = calculate_spectrum.calculate_equilibrium_spectrum(
    spectrum_data.wavenumber_array,
    molecule_simulator, distributions.EquilibriumDistribution(),
    pressure_atm=10/1013.25,
    molecular_fraction=1e-4,
    temperature_gas_kelvin=300,
    path_length_cm=30
)

plt.plot(spectrum_data.wavenumber_array, spectrum_data.data_array)
plt.title('N2O')
plt.xlabel('wavenumber / cm$^{-1}$')
plt.ylabel('transmittance')
plt.grid(True)
plt.show()