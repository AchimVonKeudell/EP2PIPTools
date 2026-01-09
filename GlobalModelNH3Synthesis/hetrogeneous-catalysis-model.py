import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from kinetic_model import SurfaceKineticModel

# reactor configurations
surface_site_density_cm2 = 1e15
volume_to_surface_ratio_cm = 1.0

# condition
pressure_mbar = 1000.
conditions = {
    'surface_temperature_kelvin': 273.15 + 500
}
conditions['gas_temperature_kelvin'] = conditions['surface_temperature_kelvin']
density_cm3 = (pressure_mbar * 1e2) / (1.38e-23 * conditions['surface_temperature_kelvin']) * 1e-6

# reaction set
file = f"{os.path.dirname(__file__)}\\Input\\reaction_sets\\NH3-Dumesic-and-Trevino-1989-ModelTable2.csv"
print(f'{file=}')
with open(file, 'r') as file_data:
    reactions_raw_string = file_data.readlines()

model = SurfaceKineticModel(reactions_raw_string,
                            surface_to_volume_ratio=surface_site_density_cm2 / volume_to_surface_ratio_cm)

initial_surface_coverages = {
    's': 0., 'Hs': 0., 'H2s': 0., 'Ns': 1.0, 'N2s': 0.0, 'NHs': 0., 'NH2s': 0., 'NH3s': 0.,
    'N2': 0.0, 'H2': 1.0 * density_cm3, 'NH3': 0.
}

t_values = np.logspace(-5, 4)
sol = solve_ivp(
    model.ode, t_span=(t_values[0], t_values[-1]), t_eval=t_values,
    y0=[initial_surface_coverages[key] for key in model.species], args=(conditions,), method='Radau', atol=1e-24
)
print(sol)

fig, ax = plt.subplots(2, 1, sharex='all')
for i, specie in enumerate(model.species):
    if 's' in specie:
        ax[1].loglog(t_values, sol.y[i], label=specie)
    else:
        ax[0].loglog(t_values, sol.y[i], label=specie)
ax[0].set_title('T$_{surface}$ = %.1f C, p = %.1f mbar' % (
        conditions['surface_temperature_kelvin']-273.15, pressure_mbar))
ax[0].set(xlim=(t_values[0], t_values[-1]), ylim=(1, 2*density_cm3), ylabel='density [cm$^{-3}$]')
ax[1].set(ylim=(1e-10, 2), xlabel='time [s]', ylabel='surface coverage [-]')
ax[0].legend(loc=(1.01, 0.0)), ax[1].legend(loc=(1.01, 0.0))
[axi.grid(True) for axi in ax]
fig.tight_layout()
plt.show()
