import matplotlib.pyplot as plt
from spectrumsimulator.simulators import MoleculeSimulator, IsotopologueSimulator, StructureSimulator, LineSimulator
from spectrumsimulator.database.repositories import MoleculeRepository, IsotopologueRepository, \
    StructureLinearMoleculeRepository, LineLinearMoleculeRepository, PartitionSumRepository
from spectrumsimulator.database.filters import MoleculeIdMoleculeQueryFilter, QueryFilterComposite, \
    MoleculeIdIsotopologueQueryFilter, MostAbundantIsotopologueQueryFilter, IsotopologueIdStructureQueryFilter, \
    WavenumberRangeStructureQueryFilter, StrongestStructureQueryFilter, StructureIdLineQueryFilter, \
    WavenumberRangeLineQueryFilter
from spectrumsimulator.database.hitran import HitranDatabase
from spectrumsimulator.spectrum import Spectrum
from spectrumsimulator.conditions import *
from spectrumsimulator.distributions import Distribution, BoltzmannDiatomicFractionalDistribution, TreanorDiatomicFractionalDistribution
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
wavenumber_start = 1800
wavenumber_end = 2400
wavenumber_resolution = 0.01

molecule_filter = MoleculeIdMoleculeQueryFilter(5)          # 5: CO
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

# Define the conditions of the measurement
conditions = TransmissionMeasurementGasConditions(
    pressure_atm=1,
    temperature_K=800,          # gas temperature
    path_length_cm=0.25
)

# Retrieve the conditions for the partial gas of the given molecule
partial_conditions = conditions.get_species_conditions_from_mixing_ratio(1)


# Set a Treanor distribution for the isotopologues
for index, CO_isotopologue_dunham_matrix in enumerate(CO_dunham_matrices):
    distribution = TreanorDiatomicFractionalDistribution(
        rotational_statistical_weight=CO_statistical_weights[index],
        dunham_matrix=CO_isotopologue_dunham_matrix
    )
    distribution.set_temperature_K(conditions.temperature_K, 3000)          # rotational and vibrational temperatures
    molecule_simulator._isotopologue_simulators[index].set_distribution(distribution)

# Create a spectrum object to be able to simulate a spectrum
spectrum = Spectrum(wavenumber_start, wavenumber_end, wavenumber_resolution)

# Add the transmission data to the spectrum object, and define the thresholds for line strength inclusion and the
# threshold at which the simulation of the line will stop
molecule_simulator.add_transmission_to_spectrum(spectrum, partial_conditions, inclusion_strength_threshold=1e-25,
                                                simulation_transmission_threshold=0.9999)

# Plot the transmission data for the given molecule
plt.plot(spectrum.wavenumber_array, spectrum.data_array, label="Transmission")

# Show the plot with all the molecules
plt.show()
