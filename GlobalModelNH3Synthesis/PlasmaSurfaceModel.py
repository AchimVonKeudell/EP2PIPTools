import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from kinetic_model import *

# TODO include LookUpTables for e-impact reactions obtained from EEDF calculations
file = f"{os.path.dirname(__file__)}\\Input\\reaction_sets\\NH3-plasma-surface-reactions-metal.csv"
print(f'{file=}')
with open(file, 'r') as file_data:
    reactions_raw_string = file_data.readlines()

conditions = {
    'pressure_mbar': 1000,
    'surface_temperature_kelvin': 300,
    'gas_temperature_kelvin': 300,
    'volume_to_surface_area_ratio_cm': 1.0,
    'surface_site_density_cm2': 1e15,
    'characteristic_diffusion_length_cm': 1.0 / 2,
    'electron_temperature_ev': 2.0,
    'electron_density_cm3': 1e10
}
conditions['density_cm3'] = 1e-6 * (
    1e2 * (conditions['pressure_mbar']) /
    (1.38e-23 * conditions['gas_temperature_kelvin'])
)

model = PlasmaKineticModel(reactions_raw_string, conditions)

initial_concentrations = {
    's': 0.05, 'Hs': 0., 'H2s': 0., 'Ns': 0.95, 'N2s': 0.0, 'NHs': 0., 'NH2s': 0., 'NH3s': 0.,
    'N2': 0.25 * conditions['density_cm3'], 'N': 0., 'H2': 0.75 * conditions['density_cm3'], 'H': 0,
    'NH': 0., 'NH2': 0., 'NH3': 0.}
initial_densities = [
    initial_concentrations[specie] for specie in model.species
]   # [1/cm3] re-configure to list following the order of model.species

t_values = np.logspace(-11, 1)
sol = solve_ivp(
    model.ode, t_span=(t_values[0], t_values[-1]), t_eval=t_values,
    y0=initial_densities, method='Radau', atol=1e-21
)
print(sol)

fig, ax = plt.subplots(2, 1, sharex='all')
for i, specie in enumerate(model.species):
    if 's' in specie:
        ax[1].loglog(t_values, sol.y[i], label=specie)
    else:
        ax[0].loglog(t_values, sol.y[i], label=specie)
ax[0].set_title('T$_{surface}$ = %.1f C, p = %.1f mbar' % (
        conditions['surface_temperature_kelvin']-273.15, conditions['pressure_mbar']))
ax[0].set(xlim=(t_values[0], t_values[-1]), ylim=(1, 2*conditions['density_cm3']), ylabel='density [cm$^{-3}$]')
ax[1].set(ylim=(1e-10, 2), xlabel='time [s]', ylabel='surface coverage [-]')
ax[0].legend(loc=(1.01, 0.0)), ax[1].legend(loc=(1.01, 0.0))
[axi.grid(True) for axi in ax]
fig.tight_layout()
plt.show()
