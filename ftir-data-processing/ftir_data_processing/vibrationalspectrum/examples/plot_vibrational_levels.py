from matplotlib import pyplot as plt

from spectrumsimulator.distributions import TreanorTriatomicLinearFermiResonantFractionalDistribution
from spectrumsimulator.constants import *
from spectrumsimulator.database import models

plt.style.use('paper')
fig, ax1 = plt.subplots(figsize=(5, 4))
plt_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

distribution = TreanorTriatomicLinearFermiResonantFractionalDistribution(
    rotational_statistical_weight=CO2_statistical_weights[0],
    molecular_constants=CO2_molecular_properties,
    vibrational_constants=CO2_vibrational_constants[0],
    rotational_constants=CO2_rotational_constants[0],
    fermi_resonance_constants=CO2_resonance_constants[0]
)

for v1 in range(1, 8):
    for ranking_index in range(1, v1+2):
        quanta = models.TriatomicLinearFermiResonantGlobalQuanta(
            vibrational_quantum_1=v1,
            vibrational_quantum_2=0,
            vibrational_quantum_3=0,
            vibrational_angular_momentum=0,
            vibrational_ranking_index=ranking_index
        )
        energy = distribution._get_unperturbed_vibrational_energy(quanta)
        energy_fermi = distribution.get_vibrational_energy(quanta)

        ax1.plot([0, 1], [energy_fermi, energy_fermi], color=plt_colors[0])
        ax1.plot([0, 1], [energy, energy], ':', color=plt_colors[3])

for v2 in range(1, 16):
    quanta = models.TriatomicLinearFermiResonantGlobalQuanta(
        vibrational_quantum_1=0,
        vibrational_quantum_2=v2,
        vibrational_quantum_3=0,
        vibrational_angular_momentum=v2,
        vibrational_ranking_index=1
    )
    energy = distribution.get_vibrational_energy(quanta)

    ax1.plot([1.5, 2.5], [energy, energy], color=plt_colors[1])

for v3 in range(1, 5):
    quanta = models.TriatomicLinearFermiResonantGlobalQuanta(
        vibrational_quantum_1=0,
        vibrational_quantum_2=0,
        vibrational_quantum_3=v3,
        vibrational_angular_momentum=0,
        vibrational_ranking_index=1
    )
    energy = distribution.get_vibrational_energy(quanta)

    ax1.plot([3, 4], [energy, energy], color=plt_colors[2])

ax1.set_xlabel("")
ax1.set_ylabel(r"Energy (cm$^{-1}$)")
ax1.set_xlim(-0.25, 4.25)
ax1.set_ylim(0, 10000)
ax1.set_xticks([0.5, 2, 3.5])
ax1.set_xticklabels([r"$\nu_1$", r"$\nu_2$", r"$\nu_3$"])
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 5))

ax2 = ax1.twinx()
ax2.set_ylabel(r"Energy (eV)")
ax2.set_ylim(0, 10000/8065.5)

plt.savefig('co2_vibrational_levels.pdf', dpi=150)
