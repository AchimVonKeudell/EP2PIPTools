import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from kinetic_model import *

file = f"{os.path.dirname(__file__)}\\Input\\reaction_sets\\helium-model_Martens-et-al-2008.csv"
print(f'{file=}')
with open(file, 'r') as file_data:
    reactions_raw_string = file_data.readlines()

conditions = {
    'pressure_mbar': 8,
    'surface_temperature_kelvin': 300,
    'gas_temperature_kelvin': 300,
    'electron_temperature_ev': 0.1,
    'electron_density_cm3': 1e10
}
conditions['density_cm3'] = ((conditions['pressure_mbar'] * 101.325) /
                             (1.38e-23 * conditions['gas_temperature_kelvin']) * 1e-6)

model = PlasmaKineticModel(reactions_raw_string)

initial_concentrations = {
    'He': .98, 'He*': 0.01, 'He+': 0., 'He2*': 0., 'He2+': 0.01, 'N2': 0.0, 'N2+': 0., 'N4+': 0.
}
initial_densities = [
    conditions['density_cm3'] * initial_concentrations[specie] for specie in model.species
]   # [1/cm3] re-configure to list following the order of model.species

t_values = np.logspace(-11, 1)
sol = solve_ivp(
    model.ode, t_span=(t_values[0], t_values[-1]), t_eval=t_values,
    y0=initial_densities, args=(conditions,), method='Radau', atol=1e-21
)
print(sol)

for i, specie in enumerate(model.species):
    plt.loglog(t_values, sol.y[i], label=specie,
               ls='-' if i < 2 else '--')
plt.legend(loc=3)
plt.ylim([1, conditions['density_cm3']*2.])
plt.show()
