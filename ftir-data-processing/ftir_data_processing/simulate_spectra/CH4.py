import matplotlib.pyplot as plt
from context import *
import time

hitran_database_file = r"D:\Promotie\data\database\ch4_database.sqlite"
spectrum_data = spectrum.Spectrum.from_wavenumber_array_params(2700, 3300, 0.005)

molecule_simulator = GetSimulator.get_CH4_simulator(
    hitran_database_file, spectrum_data.wavenumber_array.min(), spectrum_data.wavenumber_array.max()
)
molecule_simulator.set_instrumental_broadening_coefficient(1.05)

for gas_temperature_kelvin in [293]:
    start_time = time.time()
    spectrum_data.data_array = calculate_spectrum.calculate_equilibrium_spectrum(
        spectrum_data.wavenumber_array,
        molecule_simulator, distributions.EquilibriumDistribution(),
        pressure_atm=4e-3 / 1.01325,
        molecular_fraction=.5 / 25.5,
        temperature_gas_kelvin=300,
        path_length_cm=42.6
    )
    end_time = time.time()
    print(f'{end_time - start_time:.1f} seconds')
    plt.plot(spectrum_data.wavenumber_array, spectrum_data.data_array)
plt.title('CH4')
plt.xlabel('wavenumber / cm$^{-1}$')
plt.ylabel('transmittance')
plt.grid(True)
plt.show()