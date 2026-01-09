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
from spectrumsimulator.distributions import *

# Read the custom HITRAN database file format and create repository objects for the database to retrieve the data
# from the database easily.
db = HitranDatabase("../hitran_database.sqlite")
repo_molecules = MoleculeRepository(db)
repo_isotopologues = IsotopologueRepository(db)
repo_structures = StructureLinearMoleculeRepository(db)
repo_lines = LineLinearMoleculeRepository(db)
repo_partition_sum = PartitionSumRepository(db)

# Define selection criteria for the molecules that should be included. This includes the id of the molecules that should
# be included and the amount of vibrational structures that should be included for each isotopologue of the molecule.
# The series [10, 5, 2] for the structure_count property means that for the most abundant isotopologue, 10 structures
# are included, while for the second most abundant isotopologue, 5 structures are included, and so forth. Also, the
# mixing ratio (fraction of the total gas) of each of the molecules is defined here.
selection_criteria = [
    {
        'molecule_id': 2,
        'structure_count': [25, 10, 5, 2, 1],
        'mixing_ratio': 0.5,
    }
]

# Define the wavenumber range for which to perform the simulation
wavenumber_start = 1800
wavenumber_end = 2400
wavenumber_resolution = 0.01

# Create some variables to store the simulators and keep track of the amount of absorption lines that are included
line_simulator_cache = []
line_count = 0

# Build a tree of nested simulators for molecules, isotopologues, structures and lines based on the selection criteria
# defined above.
for selection_criterion in selection_criteria:
    molecule_filter = MoleculeIdMoleculeQueryFilter(selection_criterion['molecule_id'])
    molecule = repo_molecules.get_filtered(molecule_filter).all()[0]
    molecule_simulator = MoleculeSimulator(molecule)

    isotopologue_filter = QueryFilterComposite()
    isotopologue_filter.add(MoleculeIdIsotopologueQueryFilter(molecule.molecule_id))
    isotopologue_filter.add(MostAbundantIsotopologueQueryFilter(len(selection_criterion['structure_count'])))

    for isotopologue in repo_isotopologues.get_filtered(isotopologue_filter):
        partition_sum_set = repo_partition_sum.get_by_isotopologue(isotopologue)
        isotopologue_simulator = IsotopologueSimulator(isotopologue, partition_sum_set)

        structure_filter = QueryFilterComposite()
        structure_filter.add(IsotopologueIdStructureQueryFilter(isotopologue.molecule_id, isotopologue.isotopologue_order_num))
        structure_filter.add(WavenumberRangeStructureQueryFilter(wavenumber_start, wavenumber_end))
        structure_filter.add(StrongestStructureQueryFilter(selection_criterion['structure_count'][isotopologue.isotopologue_order_num-1]))

        for structure in repo_structures.get_filtered(structure_filter):
            structure_simulator = StructureSimulator(structure)

            line_filter = QueryFilterComposite()
            line_filter.add(StructureIdLineQueryFilter(structure.structure_id))
            line_filter.add(WavenumberRangeLineQueryFilter(wavenumber_start, wavenumber_end))

            for line in repo_lines.get_filtered(line_filter):
                line_simulator = LineSimulator(line)
                line_count += 1

                structure_simulator.add_line_simulator(line_simulator)
            isotopologue_simulator.add_structure_simulator(structure_simulator)
        molecule_simulator.add_isotopologue_simulator(isotopologue_simulator)
    line_simulator_cache.append(molecule_simulator)

print("%i lines selected" % line_count)

# Define the conditions of the measurement
conditions = TransmissionMeasurementGasConditions(
    pressure_atm=1,
    temperature_K=296,          # gas temperature
    path_length_cm=1
)

# Define the distribution for the energy levels
# This type instructs the simulators to rely on HITRAN data
distribution = EquilibriumDistribution()
distribution.set_temperature_K(conditions.temperature_K)        # set equilibrium temperature (for cross sections)


# Loop through each molecule simulator separately so that the transmission spectrum of each molecule can be calculated
# individually
for criteria_index, molecule_simulator in enumerate(line_simulator_cache):
    # Set the distribution for the molecule simulator
    molecule_simulator.set_distribution(distribution)

    # Retrieve the conditions for the partial gas of the given molecule
    partial_conditions = conditions.get_species_conditions_from_mixing_ratio(selection_criteria[criteria_index]['mixing_ratio'])

    # Create a spectrum object to be able to simulate a spectrum
    spectrum = Spectrum(wavenumber_start, wavenumber_end, wavenumber_resolution)

    # Add the transmission data to the spectrum object, and define the thresholds for line strength inclusion and the
    # threshold at which the simulation of the line will stop
    molecule_simulator.add_transmission_to_spectrum(spectrum, partial_conditions, inclusion_strength_threshold=1e-25, simulation_transmission_threshold=0.9999)

    # Plot the transmission data for the given molecule
    plt.plot(spectrum.wavenumber_array, spectrum.data_array, label="Cross section")

# Show the plot with all the molecules
plt.show()
