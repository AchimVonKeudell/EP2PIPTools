import numpy as np
import matplotlib.pyplot as plt
from nh3_model.combining_surface_and_volume_kinetics import solve_ode, Conditions, VolumeRateConstants, SurfaceReactions
from nh3_model.GasKinetics.ElectronImpactReactions import get_look_up_table


flow_rates = [200, 5, 5]  # He, N2, H2
pressure_torr = 25
surface_temperature_kelvin = 273 + 50


plasma_conditions = Conditions(
    pressure_atm=pressure_torr/760,
    electron_density_cm3=1e9,
    gas_temperature_kelvin=surface_temperature_kelvin, vibrational_temperature_n2_kelvin=surface_temperature_kelvin,
    electron_temperature_ev=1.0,
    active_surface_site_density_per_cm2=1e15
)

electron_energy_distribution = get_look_up_table('maxwellian')
gas_rate_coefficients = VolumeRateConstants.from_conditions(
    rate_constant_data_file='Input\\gas_phase_reactions.dat',
    gas_temperature_kelvin=surface_temperature_kelvin,
    electron_temperature_ev=electron_energy_distribution.electron_temperature(
        plasma_conditions.electron_temperature_ev, 0.),
    electron_impact_reaction_rates=electron_energy_distribution.get_rates(
        plasma_conditions.electron_temperature_ev, 0.)
)

time_values = np.logspace(-6, 4, 100)
fig, ax = plt.subplots(1, 1, sharex='all')

for surface, line_style in {'metal': '-'}.items():
    surface_rate_coefficients = SurfaceReactions.from_conditions(
        sticking_coefficients_data_file=f'Input\\sticking-probabilities_{surface}.dat',
        energy_barriers_data_file=f'Input\\surface-energy-barrier_{surface}.dat',
        gas_temperature_kelvin=surface_temperature_kelvin
    )

    solution = solve_ode(
        time_values=time_values,
        initial_values=[0, 0.5/200, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        conditions=plasma_conditions, gas_kinetic_rates=gas_rate_coefficients,
        surface_kinetic_rates=surface_rate_coefficients
    )
    for label, i in {'N': 8, 'H': 9, 'NH': 10, 'NH2': 11}.items():
        ax.loglog(time_values, solution[i], label=label, ls=line_style)

ax.grid(True)
ax.set(xlim=(1e-6, time_values[-1]), xlabel='times [s]',
          ylim=(1e3/plasma_conditions.active_surface_site_density_per_cm2, 2))
ax.legend(loc=3)
plt.show()