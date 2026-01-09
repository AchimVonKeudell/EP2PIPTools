import matplotlib.pyplot as plt
from spectrumsimulator.simulators import MoleculeSimulator, IsotopologueSimulator, StructureSimulator, LineSimulator
from spectrumsimulator.database.repositories import MoleculeRepository, IsotopologueRepository, \
    StructureLinearMoleculeRepository, LineLinearMoleculeRepository, PartitionSumRepository
from spectrumsimulator.database.filters import MoleculeIdMoleculeQueryFilter, QueryFilterComposite, \
    MoleculeIdIsotopologueQueryFilter, MostAbundantIsotopologueQueryFilter, IsotopologueIdStructureQueryFilter, \
    WavenumberRangeStructureQueryFilter, StrongestStructureQueryFilter, StructureIdLineQueryFilter, \
    WavenumberRangeLineQueryFilter, TriatomicLinearFermiResonantVibrationalQuantum3LowerStructureQueryFilter
from spectrumsimulator.database.hitran import HitranDatabase
from spectrumsimulator.spectrum import Spectrum
from spectrumsimulator.conditions import *
from spectrumsimulator.distributions import EquilibriumDistribution, BoltzmannTriatomicLinearFermiResonantFractionalDistribution, TreanorTriatomicLinearFermiResonantFractionalDistribution
from spectrumsimulator.constants import *

# Read the custom HITRAN database file format and create repository objects for the database to retrieve the data
# from the database easily.
db = HitranDatabase("../hitran_database.sqlite")
repo_molecules = MoleculeRepository(db)
repo_isotopologues = IsotopologueRepository(db)
repo_structures = StructureLinearMoleculeRepository(db)
repo_lines = LineLinearMoleculeRepository(db)
repo_partition_sum = PartitionSumRepository(db)

# Define the wavenumber range for which to perform the simulation
wavenumber_start = 2260  #2260.7
wavenumber_end = 2262  #2260.9
wavenumber_resolution = 0.001

# Create a figure that can be used for plotting the simulation
plt = InteractiveTogglePlotFigure('paper')

total_spectrum = Spectrum(wavenumber_start, wavenumber_end, wavenumber_resolution)
total_spectrum.data_array = np.ones_like(total_spectrum.wavenumber_array)

import time

v3_max = 5
v3_simulators = []

for v3 in range(v3_max):
    molecule_filter = MoleculeIdMoleculeQueryFilter(2)
    molecule = repo_molecules.get_filtered(molecule_filter).all()[0]
    molecule_simulator = MoleculeSimulator(molecule)

    isotopologue_filter = QueryFilterComposite()
    isotopologue_filter.add(MoleculeIdIsotopologueQueryFilter(molecule.molecule_id))
    isotopologue_filter.add(MostAbundantIsotopologueQueryFilter(3))

    for isotopologue in repo_isotopologues.get_filtered(isotopologue_filter):
        partition_sum_set = repo_partition_sum.get_by_isotopologue(isotopologue)
        isotopologue_simulator = IsotopologueSimulator(isotopologue, partition_sum_set)

        structure_filter = QueryFilterComposite()
        structure_filter.add(IsotopologueIdStructureQueryFilter(isotopologue.molecule_id, isotopologue.isotopologue_order_num))
        structure_filter.add(WavenumberRangeStructureQueryFilter(wavenumber_start, wavenumber_end))
        structure_filter.add(TriatomicLinearFermiResonantVibrationalQuantum3LowerStructureQueryFilter(v3))
        structure_filter.add(StrongestStructureQueryFilter(20))

        for structure in repo_structures.get_filtered(structure_filter):
            structure_simulator = StructureSimulator(structure)

            line_filter = QueryFilterComposite()
            line_filter.add(StructureIdLineQueryFilter(structure.structure_id))
            line_filter.add(WavenumberRangeLineQueryFilter(wavenumber_start, wavenumber_end))

            for line in repo_lines.get_filtered(line_filter):
                line_simulator = LineSimulator(line)

                structure_simulator.add_line_simulator(line_simulator)
            isotopologue_simulator.add_structure_simulator(structure_simulator)
        molecule_simulator.add_isotopologue_simulator(isotopologue_simulator)

    v3_simulators.append(molecule_simulator)

for v3 in range(v3_max):
    print(time.time())

    # Create a spectrum object to be able to simulate a spectrum
    spectrum_non_eq = Spectrum(wavenumber_start, wavenumber_end, wavenumber_resolution)

    # Define the conditions of the measurement
    conditions = TransmissionMeasurementGasConditions(
        pressure_atm=0.0068,
        temperature_K=491,
        path_length_cm=17
    )

    # Retrieve the conditions for the partial gas of the given molecule
    partial_conditions = conditions.get_species_conditions_from_mixing_ratio(0.68)

    # Set a Treanor distribution for the isotopologues
    for index in range(len(v3_simulators[v3]._isotopologue_simulators)):
        distribution = TreanorTriatomicLinearFermiResonantFractionalDistribution(
            rotational_statistical_weight=CO2_statistical_weights[index],
            molecular_constants=CO2_molecular_properties[index],
            vibrational_constants=CO2_vibrational_constants[index],
            rotational_constants=CO2_rotational_constants[index],
            fermi_resonance_constants=CO2_resonance_constants[index],
            v1_max=4,
            v2_max=6,
            v3_max=3,
            J_max=60
        )
        distribution.set_temperature_K(conditions.temperature_K, 517, 1500)
        v3_simulators[v3]._isotopologue_simulators[index].set_distribution(distribution)

    # Add the transmission data to the spectrum object, and define the thresholds for line strength inclusion and the
    # threshold at which the simulation of the line will stop
    v3_simulators[v3].add_transmission_to_spectrum(spectrum_non_eq, partial_conditions,
                                                      inclusion_strength_threshold=1e-25,
                                                      simulation_transmission_threshold=0.999)

    # Create a spectrum object to be able to simulate a spectrum
    spectrum_eq = Spectrum(wavenumber_start, wavenumber_end, wavenumber_resolution)

    # Define the conditions of the measurement
    conditions = TransmissionMeasurementGasConditions(
        pressure_atm=0.0068,
        temperature_K=300,
        path_length_cm=6
    )

    # Retrieve the conditions for the partial gas of the given molecule
    partial_conditions = conditions.get_species_conditions_from_mixing_ratio(0.68)

    # Set a Treanor distribution for the isotopologues
    for index in range(len(v3_simulators[v3]._isotopologue_simulators)):
        distribution = EquilibriumDistribution()
        distribution.set_temperature_K(conditions.temperature_K)
        v3_simulators[v3]._isotopologue_simulators[index].set_distribution(distribution)

    # Add the transmission data to the spectrum object, and define the thresholds for line strength inclusion and the
    # threshold at which the simulation of the line will stop
    v3_simulators[v3].add_transmission_to_spectrum(spectrum_eq, partial_conditions,
                                                    inclusion_strength_threshold=1e-25,
                                                    simulation_transmission_threshold=0.995)

    spectrum_non_eq.data_array *= spectrum_eq.data_array

    print(time.time())

    # Plot the transmission data for the given molecule
    plt.plot(spectrum_non_eq.wavenumber_array, spectrum_non_eq.data_array, label=r"$\nu_3$ = %i $\rightarrow$ %i" % (v3, v3+1), linewidth=2)

    # Calculate the total spectrum
    total_spectrum.data_array *= spectrum_non_eq.data_array

# Plot the transmission data for the given molecule
#plt.plot(total_spectrum.wavenumber_array, total_spectrum.data_array, label="Total", color=[0, 0, 0])

# Show the plot with all the molecules
plt.draw_git_describe()
plt.set_xtick_spacing(0.4)
plt.set_y_label("Transmission")
plt.set_x_label("Wavenumber (cm$^{-1}$)")
plt.set_xlim(min(spectrum_non_eq.wavenumber_array), max(spectrum_non_eq.wavenumber_array))
plt.tight_layout(margin=0.01)
plt.save("2260_2262_transmission_spectrum.jpg", dpi=150)
