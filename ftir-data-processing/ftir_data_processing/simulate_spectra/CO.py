import matplotlib.pyplot as plt
import numpy as np
from context import *

hitran_database_file = r"D:\Promotie\data\database\co_database-He.sqlite"


spectrum_data = spectrum.Spectrum.from_wavenumber_array_params(2000, 2250, 0.04)

molecule_simulator = GetSimulator.get_CO_simulator(
    hitran_database_file, spectrum_data.wavenumber_array.min(), spectrum_data.wavenumber_array.max()
)
molecule_simulator.set_instrumental_broadening_coefficient(0.15)

spectrum_data.data_array = calculate_spectrum.calculate_equilibrium_spectrum(
    spectrum_data.wavenumber_array,
    molecule_simulator, distributions.EquilibriumDistribution(),
    pressure_atm=1,
    molecular_fraction=1e-4,
    temperature_gas_kelvin=300,
    path_length_cm=1
)

print(np.concatenate([[spectrum_data.wavenumber_array], [spectrum_data.data_array]]).T.shape)
plt.rc('axes', labelweight='bold')

plt.plot(spectrum_data.wavenumber_array, spectrum_data.data_array)
# plt.plot([1800, 2400], [.9999]*2, ls='--', c='black')
plt.xlim(2000, 2250)
# plt.ylim(1-1.1e-3, 1+1e-4)
plt.xlabel('wavenumber / cm$^{-1}$')
plt.ylabel('transmittance / -')

np.savetxt(
    fname=r"D:\Promotie\research output\01_presentations+posters\20260112 - talk-Nijmegen\figures\example_CO-gas.dat",
    X=np.concatenate([[spectrum_data.wavenumber_array], [spectrum_data.data_array]]).T,
    header='1 atm, 1e-4, 300K, 1cm\nwavenumber\ttransmittance'
)

plt.show()