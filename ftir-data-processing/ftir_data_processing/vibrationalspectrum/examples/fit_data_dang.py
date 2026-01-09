import numpy as np
from lmfit import Model, Parameters, Parameter
import matplotlib.pyplot as plt
from ..spectrumsimulator.simulators import MoleculeSimulator, IsotopologueSimulator, StructureSimulator, LineSimulator
from ..spectrumsimulator.database.repositories import MoleculeRepository, IsotopologueRepository, \
    StructureLinearMoleculeRepository, LineLinearMoleculeRepository, PartitionSumRepository
from ..spectrumsimulator.database.filters import MoleculeIdMoleculeQueryFilter, QueryFilterComposite, \
    MoleculeIdIsotopologueQueryFilter, MostAbundantIsotopologueQueryFilter, IsotopologueIdStructureQueryFilter, \
    WavenumberRangeStructureQueryFilter, StrongestStructureQueryFilter, StructureIdLineQueryFilter, \
    WavenumberRangeLineQueryFilter, TriatomicLinearFermiResonantVibrationalQuantum3LowerStructureQueryFilter
from ..spectrumsimulator.database.hitran import HitranDatabase
from ..spectrumsimulator.spectrum import Spectrum
from ..spectrumsimulator.conditions import *
from ..spectrumsimulator.distributions import TreanorTriatomicLinearFermiResonantFractionalDistribution
from ..spectrumsimulator.constants import *

####
# LOAD HITRAN DATABASE
####

# Read the custom HITRAN database file format and create repository objects for the database to retrieve the data
# from the database easily.
db = HitranDatabase("../hitran_database.sqlite")
repo_molecules = MoleculeRepository(db)
repo_isotopologues = IsotopologueRepository(db)
repo_structures = StructureLinearMoleculeRepository(db)
repo_lines = LineLinearMoleculeRepository(db)
repo_partition_sum = PartitionSumRepository(db)


####
# DEFINE THE WAVENUMBER RANGE OF THE DATA AND FIT
####

# Create a spectrum object to be able to load the data into it
x, y = np.loadtxt("data_dang_1982.csv", delimiter=",", unpack=True, skiprows=1)
wavenumber_start = min(x)
wavenumber_end = max(x)
wavenumber_resolution = (wavenumber_end - wavenumber_start) / len(x)
spectrum_dang = Spectrum(wavenumber_start, wavenumber_end, wavenumber_resolution)
spectrum_dang.data_array = np.interp(spectrum_dang.wavenumber_array, x, y)


####
# SELECT RELEVANT LINES FOR CO2
####

molecule_filter = MoleculeIdMoleculeQueryFilter(2)
molecule = repo_molecules.get_filtered(molecule_filter).all()[0]
molecule_simulator = MoleculeSimulator(molecule)

isotopologue_filter = QueryFilterComposite()
isotopologue_filter.add(MoleculeIdIsotopologueQueryFilter(molecule.molecule_id))
isotopologue_filter.add(MostAbundantIsotopologueQueryFilter(3))

for isotopologue in repo_isotopologues.get_filtered(isotopologue_filter):
    partition_sum_set = repo_partition_sum.get_by_isotopologue(isotopologue)
    isotopologue_simulator = IsotopologueSimulator(isotopologue, partition_sum_set)

    for v3 in range(5):
        structure_filter = QueryFilterComposite()
        structure_filter.add(IsotopologueIdStructureQueryFilter(isotopologue.molecule_id, isotopologue.isotopologue_order_num))
        structure_filter.add(WavenumberRangeStructureQueryFilter(wavenumber_start, wavenumber_end))
        structure_filter.add(TriatomicLinearFermiResonantVibrationalQuantum3LowerStructureQueryFilter(v3))
        structure_filter.add(StrongestStructureQueryFilter(15))

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


####
# DEFINE FITTING FUNCTION
####

def create_spectrum(wavenumbers, pressure_atm, path_length_cm, fraction_CO2, temperature_rot_K, temperature_vib_12_K, temperature_vib_3_K):
    # Create a spectrum object to be able to simulate a spectrum
    spectrum_non_eq = Spectrum(wavenumbers, np.ones_like(wavenumbers))

    # Define the conditions of the measurement
    conditions = TransmissionMeasurementGasConditions(
        pressure_atm=pressure_atm,
        temperature_K=temperature_rot_K,
        path_length_cm=path_length_cm
    )

    # Retrieve the conditions for the partial gas of the given molecule
    partial_conditions = conditions.get_species_conditions_from_mixing_ratio(fraction_CO2)

    # Set a Treanor distribution for the isotopologues
    for index in range(len(molecule_simulator._isotopologue_simulators)):
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
        distribution.set_temperature_K(conditions.temperature_K, temperature_vib_12_K, temperature_vib_3_K)
        molecule_simulator._isotopologue_simulators[index].set_distribution(distribution)

    # Add the transmission data to the spectrum object, and define the thresholds for line strength inclusion and the
    # threshold at which the simulation of the line will stop
    spectrum_non_eq.data_array *= molecule_simulator.add_transmission_to_spectrum(
        spectrum_non_eq, partial_conditions,
        inclusion_strength_threshold=1e-25, simulation_transmission_threshold=0.9995)

    return spectrum_non_eq.data_array


####
# PERFORM FIT AND PRINT RESULT TO CONSOLE
####

model = Model(create_spectrum, independent_vars=['wavenumbers'])
params = Parameters()
params.add("pressure_atm", 0.02066, vary=False)
params.add("path_length_cm", 10, vary=False)
params.add("fraction_CO2", 0.05, min=0, max=1)
params.add("temperature_rot_K", 400)
params.add("temperature_vib_12_K", 400)
params.add("temperature_vib_3_K", 800)

result = model.fit(spectrum_dang.data_array, params, wavenumbers=spectrum_dang.wavenumber_array)
print(result.fit_report())


####
# PLOT RESULT
####


# Plot the transmission data for the given molecule
plt.plot(x, y, 'k.', label="Data")

# Plot the fit and its initial guess
plt.plot(spectrum_dang.wavenumber_array, result.best_fit, '-', label="Fit")
plt.plot(spectrum_dang.wavenumber_array, result.residual, label="Fit")

# Show the plot with all the molecules
plt.gca().set_xlabel(r"Wavenumbers (cm$^{-1}$)")
plt.gca().set_ylabel(r"Transmittance")
plt.gca().set_xlim(2284.235, 2284.615)

