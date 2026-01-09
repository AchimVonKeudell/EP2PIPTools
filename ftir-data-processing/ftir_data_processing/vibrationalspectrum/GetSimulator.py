import os.path
from .spectrumsimulator.simulators import MoleculeSimulator, IsotopologueSimulator, StructureSimulator, LineSimulator
from .spectrumsimulator.database.repositories import MoleculeRepository, IsotopologueRepository, PartitionSumRepository, \
    StructureLinearMoleculeRepository, LineLinearMoleculeRepository, \
    StructureNonLinearMoleculeRepository, LineNonLinearMoleculeRepository, \
    StructurePyramidalTetratomicMoleculeRepository, LinePyramidalTetratomicMoleculeRepository
from .spectrumsimulator.database.filters import MoleculeIdMoleculeQueryFilter, QueryFilterComposite, \
    MoleculeIdIsotopologueQueryFilter, MostAbundantIsotopologueQueryFilter, IsotopologueIdStructureQueryFilter, \
    WavenumberRangeStructureQueryFilter, StrongestStructureQueryFilter, StructureIdLineQueryFilter,\
    WavenumberRangeLineQueryFilter, TriatomicLinearFermiResonantVibrationalQuantum3LowerStructureQueryFilter, \
    DiatomicVibrationalQuantumLowerStructureQueryFilter, PyramidalTetratomicStructureQueryFilter
from .spectrumsimulator.database.hitran import HitranDatabase


def get_simulator(hitran_database_file, wavenumber_start, wavenumber_end, molecule_id, structure_filters,
                  molecule_type='standard', *, number_of_most_abundant_isotopologues=3):
    """
    Get a simulator object that includes absorption lines for the given molecule between the start and stop wavenumber.
    As an optional parameter, a nested array of filters can be provided that filter the structure list for a given
    combination of the isotopologue order number and the v quantum.
    :param hitran_database_file:
    :param wavenumber_start:
    :param wavenumber_end:
    
    :param molecule_id:
    :param structure_filters:
    :param molecule_type:
    :param number_of_most_abundant_isotopologues:
    :return:
    """
    if not os.path.exists(hitran_database_file):
        raise FileExistsError(f'{hitran_database_file=} does not exist')

    def structure_and_line_repository(_db):
        if molecule_type == 'standard':     # i.e. linear molecule
            """The standard is meant for a 'linear' molecule. This mode must also be used for simulating ro-vibrational
            equilibrium spectra
            """
            return StructureLinearMoleculeRepository(_db), LineLinearMoleculeRepository(_db)
        elif molecule_type == 'triatomic non linear':
            return StructureNonLinearMoleculeRepository(_db), LineNonLinearMoleculeRepository(_db)
        elif molecule_type == 'pyramidal tetratomic':
            return StructurePyramidalTetratomicMoleculeRepository(_db), LinePyramidalTetratomicMoleculeRepository(_db)
        raise ValueError(f'{molecule_type=} not known')

    # Read the custom HITRAN database file format and create repository objects for the database to retrieve the data
    # from the database easily.
    db = HitranDatabase(hitran_database_file)
    repo_molecules = MoleculeRepository(db)
    repo_isotopologues = IsotopologueRepository(db)
    repo_structures, repo_lines = structure_and_line_repository(db)
    repo_partition_sum = PartitionSumRepository(db)

    # Select relevant lines for the molecule
    molecule_filter = MoleculeIdMoleculeQueryFilter(molecule_id)
    molecule = repo_molecules.get_filtered(molecule_filter).all()[0]
    molecule_simulator = MoleculeSimulator(molecule)

    isotopologue_filter = QueryFilterComposite()
    isotopologue_filter.add(MoleculeIdIsotopologueQueryFilter(molecule.molecule_id))
    isotopologue_filter.add(MostAbundantIsotopologueQueryFilter(number_of_most_abundant_isotopologues))

    for isotopologue in repo_isotopologues.get_filtered(isotopologue_filter):
        partition_sum_set = repo_partition_sum.get_by_isotopologue(isotopologue)
        isotopologue_simulator = IsotopologueSimulator(isotopologue, partition_sum_set)
        isotopologue_index = isotopologue.isotopologue_order_num-1

        for v in range(len(structure_filters[isotopologue_index])):
            structure_filter = QueryFilterComposite()
            structure_filter.add(IsotopologueIdStructureQueryFilter(isotopologue.molecule_id, isotopologue.isotopologue_order_num))
            structure_filter.add(WavenumberRangeStructureQueryFilter(wavenumber_start, wavenumber_end))
            for v_structure_filter in structure_filters[isotopologue_index][v]:
                structure_filter.add(v_structure_filter)

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
    return molecule_simulator


def get_HCN_simulator(hitran_database_file, wavenumber_start, wavenumber_end, v_structure_count_list=[[100], [50], [50]]):
    filters = []
    for i, v_isotopologue_count_list in enumerate(v_structure_count_list):
        isotopologue_filters = []
        for j, v_structure_count in enumerate(v_isotopologue_count_list):
            isotopologue_filters.append([
                StrongestStructureQueryFilter(v_structure_count)
            ])
        filters.append(isotopologue_filters)

    return get_simulator(hitran_database_file, wavenumber_start, wavenumber_end, 23, filters,
                         number_of_most_abundant_isotopologues=3)


def get_C2H2_simulator(hitran_database_file, wavenumber_start, wavenumber_end, v_structure_count_list=[[100], [50], [50]]):
    filters = []
    for i, v_isotopologue_count_list in enumerate(v_structure_count_list):
        isotopologue_filters = []
        for j, v_structure_count in enumerate(v_isotopologue_count_list):
            isotopologue_filters.append([
                StrongestStructureQueryFilter(v_structure_count)
            ])
        filters.append(isotopologue_filters)

    return get_simulator(hitran_database_file, wavenumber_start, wavenumber_end, 26, filters,
                         number_of_most_abundant_isotopologues=3)


def get_CH4_simulator(hitran_database_file, wavenumber_start, wavenumber_end, v_structure_count_list=[[100], [50], [50], [25]]):
    filters = []

    for i, v_isotopologue_count_list in enumerate(v_structure_count_list):
        isotopologue_filters = []
        for j, v_structure_count in enumerate(v_isotopologue_count_list):
            isotopologue_filters.append([
                StrongestStructureQueryFilter(v_structure_count)
            ])
        filters.append(isotopologue_filters)

    return get_simulator(hitran_database_file, wavenumber_start, wavenumber_end, 6, filters,
                         number_of_most_abundant_isotopologues=4)



def get_CO_simulator(hitran_database_file, wavenumber_start, wavenumber_end,
                     v_structure_count_list=[[25, 10, 5, 2, 1], [10, 5, 2, 1, 0], [5, 2, 1, 0, 0]]):
    """
    Get a simulator object that includes absorption lines for CO between the start and stop wavenumber. As an optional
    parameter, the structure count per vibrational level can be provided.
    :param wavenumber_start:
    :param wavenumber_end:
    :param v_structure_count_list:
    :return:
    """
    filters = []

    for i, v_isotopologue_count_list in enumerate(v_structure_count_list):
        isotopologue_filters = []
        for j, v_structure_count in enumerate(v_isotopologue_count_list):
            isotopologue_filters.append([
                DiatomicVibrationalQuantumLowerStructureQueryFilter(j),
                StrongestStructureQueryFilter(v_structure_count)
            ])
        filters.append(isotopologue_filters)

    return get_simulator(hitran_database_file, wavenumber_start, wavenumber_end, 5, filters)


def get_CO2_simulator(hitran_database_file, wavenumber_start, wavenumber_end,
                      v3_structure_count_list=[[100, 80, 50, 25, 10, 5, 2], [25, 15, 5, 2, 1, 0, 0], [10, 5, 2, 1, 0, 0, 0]]):
    """
    Get a simulator object that includes absorption lines for CO2 between the start and stop wavenumber. As an optional
    parameter, the structure count per v3 level can be provided.
    :param wavenumber_start:
    :param wavenumber_end:
    :param v3_structure_count_list:
    :return:
    """
    filters = []

    for i, v3_isotopologue_count_list in enumerate(v3_structure_count_list):
        isotopologue_filters = []
        for j, v3_structure_count in enumerate(v3_isotopologue_count_list):
            isotopologue_filters.append([
                TriatomicLinearFermiResonantVibrationalQuantum3LowerStructureQueryFilter(j),
                StrongestStructureQueryFilter(v3_structure_count)
            ])
        filters.append(isotopologue_filters)

    return get_simulator(hitran_database_file, wavenumber_start, wavenumber_end, 2, filters)


def get_H2O_simulator(hitran_database_file, wavenumber_start, wavenumber_end,
                     v_structure_count_list=[[1000], [1000], [1000]]):
    filters = []
    for v_isotopologue_count_list in v_structure_count_list:
        isotopologue_filters = []
        for j, v_structure_count_list in enumerate(v_isotopologue_count_list):
            isotopologue_filters.append([
                StrongestStructureQueryFilter(v_structure_count_list)
            ])
        filters.append(isotopologue_filters)
    return get_simulator(hitran_database_file, wavenumber_start, wavenumber_end, 1, filters,
                         molecule_type='triatomic non linear')


# TODO: Include structure filter?
def get_N2O_simulator(hitran_database_file, wavenumber_start, wavenumber_end,
                      v_structure_count_list=[[1000], [1000], [1000], [1000], [1000]]):
    filters = []

    for i, v_isotopologue_count_list in enumerate(v_structure_count_list):
        isotopologue_filters = []
        for j, v_structure_count in enumerate(v_isotopologue_count_list):
            isotopologue_filters.append([
                StrongestStructureQueryFilter(v_structure_count)
            ])
        filters.append(isotopologue_filters)

    return get_simulator(hitran_database_file, wavenumber_start, wavenumber_end, 4, filters,
                         number_of_most_abundant_isotopologues=5)

def get_NO_simulator(hitran_database_file, wavenumber_start, wavenumber_end, v_structure_count_list=[[100], [50], [50]]):
    filters = []

    for i, v_isotopologue_count_list in enumerate(v_structure_count_list):
        isotopologue_filters = []
        for j, v_structure_count in enumerate(v_isotopologue_count_list):
            isotopologue_filters.append([
                StrongestStructureQueryFilter(v_structure_count)
            ])
        filters.append(isotopologue_filters)

    return get_simulator(hitran_database_file, wavenumber_start, wavenumber_end, 8, filters,
                         number_of_most_abundant_isotopologues=3)


def get_NO2_simulator(hitran_database_file, wavenumber_start, wavenumber_end, v_structure_count_list=[[100], [50]]):
    filters = []

    for i, v_isotopologue_count_list in enumerate(v_structure_count_list):
        isotopologue_filters = []
        for j, v_structure_count in enumerate(v_isotopologue_count_list):
            isotopologue_filters.append([
                StrongestStructureQueryFilter(v_structure_count)
            ])
        filters.append(isotopologue_filters)

    return get_simulator(hitran_database_file, wavenumber_start, wavenumber_end, 10, filters,
                         number_of_most_abundant_isotopologues=2)

# TODO: Include a structure filter? As is done with get_CO2_simulator and get_CO_simulator
def get_NH3_simulator(hitran_database_file, wavenumber_start, wavenumber_end,
                      v_structure_count_list=[[100], [50]]):
    filters = []
    for v_isotopologue_count_list in v_structure_count_list:
        isotopologue_filters = []
        for j, v_structure_count_list in enumerate(v_isotopologue_count_list):
            isotopologue_filters.append([
                StrongestStructureQueryFilter(v_structure_count_list)
            ])
        filters.append(isotopologue_filters)
    return get_simulator(hitran_database_file, wavenumber_start, wavenumber_end, 11, filters,
                         molecule_type='pyramidal tetratomic')


def get_O3_simulator(hitran_database_file, wavenumber_start, wavenumber_end,
                     v_structure_count_list=[[1000], [1000], [1000]]):
    filters = []
    for v_isotopologue_count_list in v_structure_count_list:
        isotopologue_filters = []
        for j, v_structure_count_list in enumerate(v_isotopologue_count_list):
            isotopologue_filters.append([
                StrongestStructureQueryFilter(v_structure_count_list)
            ])
        filters.append(isotopologue_filters)
    return get_simulator(hitran_database_file, wavenumber_start, wavenumber_end, 3, filters,
                         molecule_type='triatomic non linear')
