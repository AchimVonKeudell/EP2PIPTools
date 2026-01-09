import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from kinetic_model import PlasmaKineticModel

file = f"Input\\reaction_sets\\NH3-plasma-surface-reactions-metal.csv"
print(f'{file=}')
with open(file, 'r') as file_data:
    reactions_raw_string = file_data.readlines()


conditions = {
    'pressure_mbar': 1000,
    'surface_temperature_kelvin': 300,
    'gas_temperature_kelvin': 300,
    'volume_to_surface_area_ratio_cm': 1e-4,
    'surface_site_density_cm2': 1e15,
    'characteristic_diffusion_length_cm': 1.0 / 2,
    # 'electron_temperature_ev': 2.0,
    # 'electron_density_cm3': 0
}
conditions['density_cm3'] = 1e-6 * (
    1e2 * (conditions['pressure_mbar']) /
    (1.38e-23 * conditions['gas_temperature_kelvin'])
)

model = PlasmaKineticModel(reactions_raw_string, conditions)

initial_concentrations = {
    's': 0.0, 'Hs': 1., 'H2s': 0., 'Ns': 0., 'N2s': 0.0, 'NHs': 0., 'NH2s': 0., 'NH3s': 0.,
    'N2': 1. * conditions['density_cm3'], 'N': 0., 'H2': 0.0 * conditions['density_cm3'], 'H': 0,
    'NH': 0., 'NH2': 0., 'NH3': 0.}

initial_densities = [
    initial_concentrations[specie] for specie in model.species
]   # [1/cm3] re-configure to list following the order of model.species

reactor_volume_L = 1
total_flow_rate_L_min = 1
residence_time_s = reactor_volume_L / (total_flow_rate_L_min / 60)


def calculate_gas_density_in(influx_composition):
    y = np.zeros_like(model.species, dtype=np.float64)
    for species_index, species_label in enumerate(model.species):
        if species_label in influx_composition:
            y[species_index] = influx_composition[species_label] * conditions['density_cm3']
    return y

pulse_duration = 10e-6
repetition_frequency = 1e3

def ode(t, y):
    """
        dn/dt = reaction_vector + flux_in - flux_out
    """

    # composition = {'N2': 0.90, 'H2': 0.10, 'N': 0.0, 'H': 0.}
    # if (t % (1/repetition_frequency)) <= pulse_duration:
    composition = {'N2': 0.99, 'H2': 0., 'N': 0.01, 'H': 0.}

    gas_phase_density_in = calculate_gas_density_in(composition)

    flux_balance = np.zeros_like(y)
    for species_index, species_label in enumerate(model.species):
        if 's' in species_label:
            pass
        else:
            flux_balance[species_index] = (gas_phase_density_in - y)[species_index] / residence_time_s

    reactions_vector = model.ode(t, y)
    return reactions_vector + flux_balance


t_values = np.logspace(-11, 1, 200)
# t_values = np.linspace(0, 1)
sol = solve_ivp(
    ode, t_span=(t_values[0], t_values[-1]), t_eval=t_values,
    y0=initial_densities, method='Radau', atol=1e-21
)
print(sol)

fig, ax = plt.subplots(2, 1, sharex='all')
for i, species_label in enumerate(model.species):
    if 's' in species_label:
        ax[1].loglog(t_values, sol.y[i], label=species_label)
    else:
        ax[0].loglog(t_values, 1e2 * sol.y[i] / conditions['density_cm3'], label=species_label)

ax[0].set_title('T$_{surface}$ = %.1f C, p = %.1f mbar' % (
        conditions['surface_temperature_kelvin']-273.15, conditions['pressure_mbar']))
ax[0].set(xlim=(t_values[0], t_values[-1]),
          ylim=(1e-6, 2e2), ylabel='molar fraction / %'
          # ylim=(1e10, 2*conditions['density_cm3']), ylabel='density [cm$^{-3}$]'
          )
ax[1].set(ylim=(1e-10, 2), xlabel='time [s]', ylabel='surface coverage [-]')
ax[0].legend(loc=(1.01, 0.0)), ax[1].legend(loc=(1.01, 0.0))
[axi.grid(True) for axi in ax]
fig.tight_layout()
plt.show()
