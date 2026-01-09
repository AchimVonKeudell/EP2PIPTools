import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from kinetic_model import *

reactions_raw_string = \
    r"""1,He -> He*,ElectronImpactMaxwellianEEDF(Helium\He_e_to_He23S_e.dat)
    2,He -> He+,ElectronImpactMaxwellianEEDF(Helium\He_e_to_He-ion_e.dat)
    3,He* -> He+,ElectronImpactMaxwellianEEDF(Helium\He23S_e_to_He-ion_e.dat)
    4,He2+ -> He* + He,ElectronImpactTemperatureReaction(8.9e-9 cm3/s;1.0 eV;-1.5 -;1. K;1.5 -)
    5,N2+ -> N2,ElectronImpactTemperatureReaction(4.8e-7 cm3/s;1. eV;-0.5 -;1. K;0.5 -)
    6,N4+ -> 2 N2,ElectronImpactTemperatureReaction(2e-6 cm3/s;1 eV;-0.5 -;1 K;0.5 -)
    7,He+ + 2 He -> He2+ + He,Reaction(1.1e-31 cm6/s)
    8,He2+ + N2 -> He2* + N2+,Reaction(1.4e-9 cm3/s)
    9,N2+ + 2 N2 -> N4+ + N2,Reaction(1.9e-29 cm6/s)
    10,N2+ + N2 + He -> N4+ + He,Reaction(1.9e-29 cm6/s)
    11,N4+ + N2 -> N2+ + 2 N2,Reaction(2.5e-15 cm3/s)
    12,N4+ + He -> N2+ + N2 + He,Reaction(2.5e-15 cm3/s)
    13,He* + 2 He -> He2* + He,Reaction(2e-34 cm6/s)
    14,2 He* -> He2+,Reaction(1.5e-9 cm3/s)
    15,He* + N2 -> He + N2+,Reaction(5e-11 cm3/s)
    16,He2* -> 2 He,Reaction(1e-4 1/s)
    17,2 He2* -> He2+ + 2 He,Reaction(1.5e-9 cm3/s)
    18,He2* + N2 -> 2 He + N2+,Reaction(3e-11 cm3/s)""".split('\n')

conditions = {
    'pressure_torr': 760,
    'surface_temperature_kelvin': 300,
    'gas_temperature_kelvin': 300,
    'electron_temperature_ev': 2.0,
    'electron_density_cm3': 1e10
}
conditions['density_cm3'] = ((conditions['pressure_torr'] / 760 * 101325) /
                             (1.38e-23 * conditions['gas_temperature_kelvin']) * 1e-6)

model = PlasmaKineticModel(reactions_raw_string, conditions)

initial_concentrations = {
    'He': 1.0, 'He*': 0.0, 'He+': 0., 'He2*': 0., 'He2+': 0., 'N2': 0.0, 'N2+': 0., 'N4+': 0.
}
initial_densities = [
    conditions['density_cm3'] * initial_concentrations[specie] for specie in model.species
]   # [1/cm3] re-configure to list following the order of model.species

t_values = np.logspace(-11, 1)
sol = solve_ivp(
    model.ode, t_span=(t_values[0], t_values[-1]), t_eval=t_values,
    y0=initial_densities, method='Radau', atol=1e-21
)
print(sol)

for i, specie in enumerate(model.species):
    plt.loglog(t_values, sol.y[i], label=specie,
               ls='-' if i < 2 else '--')
plt.legend(loc=3)
plt.ylim([1, conditions['density_cm3']*2.])
plt.show()
