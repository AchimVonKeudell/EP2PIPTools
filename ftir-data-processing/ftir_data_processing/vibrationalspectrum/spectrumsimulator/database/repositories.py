import sqlite3
from pypika import Query, Table, Order, Tables, JoinType
from . import models, queries, results, filters


class MoleculeRepositoryInterface:
    def save(self, molecule: models.Molecule, commit=True):
        pass

    def get_by_id(self, molecule_id):
        pass

    def commit(self):
        pass


class IsotopologueRepositoryInterface:
    def save(self, isotopologue: models.Isotopologue, commit=True):
        pass

    def get_by_molecule_id(self, molecule_id):
        pass

    def get_by_isotopologue_order_num(self, molecule_id, isotopologue_order_num):
        pass

    def commit(self):
        pass


class StructureRepositoryInterface:
    def save(self, structure: models.Structure, commit=True):
        pass

    def get_by_id(self, structure_id):
        pass

    def get_by_isotopologue_order_num(self, molecule_id, isotopologue_order_num):
        pass

    def commit(self):
        pass


class LineRepositoryInterface:
    def save(self, line: models.Line, molecule_id, commit=True):
        pass

    def get_by_id(self, line_id):
        pass

    def get_by_structure_id(self, structure_id):
        pass

    def commit(self):
        pass


class PartitionSumRepositoryInterface:
    def get_by_isotopologue(self, isotopologue: models.Isotopologue):
        pass


class MoleculeRepository(MoleculeRepositoryInterface):
    """
    Object to abstract the database communication from the rest of the code. This object provides methods to interact
    with molecule records in the database.
    """
    _database = None

    def __init__(self, database: sqlite3.Connection):
        self._database = database
        self._molecule_query = Query().from_('molecules').select('molecule_id', 'molecule_name')

    def save(self, molecule: models.Molecule, commit=True):
        """
        Method to save a molecule object to the database
        :param molecule: Molecule object
        :param commit: Boolean indicating of a commit should be done immediately after the save operation. If this is
        set to False, the commit method has to be called to commit the changes.
        :return: Inserted molecule id
        """
        cursor = self._database.cursor()
        cursor.execute(queries.molecule_insert, (molecule.molecule_id, molecule.molecule_name, ))

        if commit:
            self._database.commit()

        molecule.molecule_id = cursor.lastrowid

        return molecule.molecule_id

    def get_all(self):
        """
        Method to retrieve all molecules from the database
        :return: 
        """
        cursor = self._database.cursor()
        cursor.execute(str(self._molecule_query))

        return results.MoleculeResultSet(cursor)

    def get_by_id(self, molecule_id):
        """
        Method to retrieve a molecule object from the database with a given id
        :param molecule_id: The id of the molecule to retrieve
        :return: 
        """
        query_filter = filters.MoleculeIdMoleculeQueryFilter(molecule_id)
        result = self.get_filtered(query_filter)

        try:
            return result.__next__()
        except StopIteration:
            raise IndexError

    def get_filtered(self, query_filter: filters.QueryFilterInterface):
        """
        Method to retrieve molecules with a filter applied to the results.
        :param query_filter: A filter object that filters the results
        :return: 
        """
        query = query_filter.filter_query(self._molecule_query)
        cursor = self._database.cursor()
        cursor.execute(str(query))

        return results.MoleculeResultSet(cursor)

    def commit(self):
        """
        Write pending changes to the database file. Typically, this is done automatically when calling save(), but this
        can be disabled by passing commit=False to the save method. By not committing after each save, database
        throughput can be maximized.
        :return: 
        """
        self._database.commit()


class IsotopologueRepository(IsotopologueRepositoryInterface):
    """
    Object to abstract the database communication from the rest of the code. This object provides methods to interact
    with isotopologue records in the database.
    """
    _database = None

    def __init__(self, database: sqlite3.Connection):
        self._database = database
        self._isotopologue_query = Query().from_('isotopologues').select(
            'molecule_id',
            'isotopologue_order_num',
            'afgl_code',
            'abundance',
            'molar_mass'
        )

    def save(self, isotopologue: models.Isotopologue, commit=True):
        """
        Method to save a isotopologue object to the database
        :param isotopologue: Isotopologue object
        :param commit: Boolean indicating of a commit should be done immediately after the save operation. If this is
        set to False, the commit method has to be called to commit the changes.
        """
        cursor = self._database.cursor()
        cursor.execute(queries.isotopologue_insert, (
            isotopologue.molecule_id, isotopologue.isotopologue_order_num, isotopologue.afgl_code,
            isotopologue.abundance, isotopologue.molar_mass,
        ))

        if commit:
            self._database.commit()

    def get_by_molecule_id(self, molecule_id):
        """
        Method to retrieve all isotopologues for a given molecule
        :param molecule_id: The id of the molecule for which to retrieve all isotopologues
        :return: 
        """
        query_filter = filters.MoleculeIdIsotopologueQueryFilter(molecule_id)
        return self.get_filtered(query_filter)

    def get_by_isotopologue_order_num(self, molecule_id, isotopologue_order_num):
        """
        Method to retrieve an isotopologue object from the database with a given id
        :param molecule_id: The id of the molecule to which the isotopologue belongs
        :param isotopologue_order_num: The abundance order number of the isotopologue for the given molecule_id where 1 is the most abundant
        :return: 
        """
        query_filter = filters.IsotopologueIdIsotopologueQueryFilter(molecule_id, isotopologue_order_num)
        result = self.get_filtered(query_filter)

        try:
            return result.__next__()
        except StopIteration:
            raise IndexError

    def get_filtered(self, query_filter: filters.QueryFilterInterface):
        """
        Method to retrieve isotopologues with a filter applied to the results.
        :param query_filter: A filter object that filters the results
        :return: 
        """
        query = query_filter.filter_query(self._isotopologue_query)
        cursor = self._database.cursor()
        cursor.execute(str(query))

        return results.IsotopologueResultSet(cursor)

    def commit(self):
        """
        Write pending changes to the database file. Typically, this is done automatically when calling save(), but this
        can be disabled by passing commit=False to the save method. By not committing after each save, database
        throughput can be maximized.
        :return: 
        """
        self._database.commit()


class StructureLinearMoleculeRepository(StructureRepositoryInterface):
    """
    Object to abstract the database communication from the rest of the code. This object provides methods to interact
    with structure records in the database.
    """
    _database = None

    def __init__(self, database: sqlite3.Connection):
        structures_table, structures_diatomic_table, structures_triatomic_linear_fermi_resonant_table = \
            Tables('structures', 'structures_diatomic', 'structures_triatomic_linear_fermi_resonant',)
        self._database = database
        self._structure_query = (
            Query().from_(structures_table)
                .join(structures_diatomic_table, how=JoinType.left_outer)
                .on(structures_table.structure_id == structures_diatomic_table.structure_id)
                .join(structures_triatomic_linear_fermi_resonant_table, how=JoinType.left_outer)
                .on(structures_table.structure_id == structures_triatomic_linear_fermi_resonant_table.structure_id)
                .select(
                    structures_table.structure_id,
                    structures_table.molecule_id,
                    structures_table.isotopologue_order_num,
                    structures_table.global_quanta_lower,
                    structures_table.global_quanta_upper,
                    structures_table.max_line_strength_296K,
                    structures_diatomic_table.vibrational_quantum_lower,
                    structures_diatomic_table.vibrational_quantum_upper,
                    structures_triatomic_linear_fermi_resonant_table.vibrational_quantum_1_lower,
                    structures_triatomic_linear_fermi_resonant_table.vibrational_quantum_2_lower,
                    structures_triatomic_linear_fermi_resonant_table.vibrational_quantum_3_lower,
                    structures_triatomic_linear_fermi_resonant_table.vibrational_angular_momentum_lower,
                    structures_triatomic_linear_fermi_resonant_table.vibrational_ranking_index_lower,
                    structures_triatomic_linear_fermi_resonant_table.vibrational_quantum_1_upper,
                    structures_triatomic_linear_fermi_resonant_table.vibrational_quantum_2_upper,
                    structures_triatomic_linear_fermi_resonant_table.vibrational_quantum_3_upper,
                    structures_triatomic_linear_fermi_resonant_table.vibrational_angular_momentum_upper,
                    structures_triatomic_linear_fermi_resonant_table.vibrational_ranking_index_upper
            )
        )

    def save(self, structure: models.Structure, commit=True):
        """
        Method to save a structure object to the database
        :param structure: Structure object
        :param commit: Boolean indicating of a commit should be done immediately after the save operation. If this is
        set to False, the commit method has to be called to commit the changes.
        """
        cursor = self._database.cursor()
        cursor.execute(queries.structure_insert_query, (
            structure.molecule_id, structure.isotopologue_order_num, structure.global_quanta_upper.hitran_raw_quanta,
            structure.global_quanta_lower.hitran_raw_quanta, structure.max_line_strength_296K,
        ))

        # if the structure is for a diatomic molecule, additional parameters have to be stored
        if models.Structure.get_global_quanta_model_type_by_molecule_id(structure.molecule_id) == models.DiatomicGlobalQuanta:
            cursor.execute(queries.structure_diatomic_insert_query, (
                cursor.lastrowid,
                structure.global_quanta_lower.vibrational_quantum,
                structure.global_quanta_upper.vibrational_quantum
            ))
        elif models.Structure.get_global_quanta_model_type_by_molecule_id(structure.molecule_id) == models.TriatomicLinearFermiResonantGlobalQuanta:
            cursor.execute(queries.structure_triatomic_linear_fermi_resonance_insert_query, (
                cursor.lastrowid,
                structure.global_quanta_lower.vibrational_quantum_1,
                structure.global_quanta_lower.vibrational_quantum_2,
                structure.global_quanta_lower.vibrational_quantum_3,
                structure.global_quanta_lower.vibrational_angular_momentum,
                structure.global_quanta_lower.vibrational_ranking_index,
                structure.global_quanta_upper.vibrational_quantum_1,
                structure.global_quanta_upper.vibrational_quantum_2,
                structure.global_quanta_upper.vibrational_quantum_3,
                structure.global_quanta_upper.vibrational_angular_momentum,
                structure.global_quanta_upper.vibrational_ranking_index
            ))

        if commit:
            self._database.commit()

        structure.structure_id = cursor.lastrowid

        return structure.structure_id

    def get_by_id(self, structure_id):
        """
        Method to retrieve a structure object from the database with a given id
        :param structure_id: The id of the structure to retrieve
        :return: 
        """
        query_filter = filters.StructureIdStructureQueryFilter(structure_id)
        result = self.get_filtered(query_filter)

        try:
            return result.__next__()
        except StopIteration:
            raise IndexError

    def get_by_isotopologue_order_num(self, molecule_id, isotopologue_order_num):
        """
        Method to retrieve all structures for a given isotopologue
        :param molecule_id: The id of the molecule to which the isotopologue belongs
        :param isotopologue_order_num: The abundance order number of the isotopologue for the given molecule_id where 1 is the most abundant
        :return: 
        """
        query_filter = filters.IsotopologueIdStructureQueryFilter(molecule_id, isotopologue_order_num)
        return self.get_filtered(query_filter)

    def get_filtered(self, query_filter: filters.QueryFilterInterface):
        """
        Method to retrieve structures with a filter applied to the results.
        :param query_filter: A filter object that filters the results
        :return: 
        """
        query = query_filter.filter_query(self._structure_query)
        cursor = self._database.cursor()
        cursor.execute(str(query))

        return results.StructureResultSet(cursor)

    def commit(self):
        """
        Write pending changes to the database file. Typically, this is done automatically when calling save(), but this
        can be disabled by passing commit=False to the save method. By not committing after each save, database
        throughput can be maximized.
        :return: 
        """
        self._database.commit()


class StructureNonLinearMoleculeRepository(StructureLinearMoleculeRepository):
    """
    Object to abstract the database communication from the rest of the code. This object provides methods to interact
    with structure records in the database.
    """
    def __init__(self, database: sqlite3.Connection):
        structures_table, structures_non_linear_table = Tables('structures', 'structures_triatomic_non_linear',)
        self._database = database
        self._structure_query = (
            Query().from_(structures_table)
                .join(structures_non_linear_table, how=JoinType.left_outer)
                .on(structures_table.structure_id == structures_non_linear_table.structure_id)
                .select(
                    structures_table.structure_id,
                    structures_table.molecule_id,
                    structures_table.isotopologue_order_num,
                    structures_table.global_quanta_lower,
                    structures_table.global_quanta_upper,
                    structures_table.max_line_strength_296K,
                    structures_non_linear_table.vibrational_quantum_1_lower,
                    structures_non_linear_table.vibrational_quantum_2_lower,
                    structures_non_linear_table.vibrational_quantum_3_lower,
                    structures_non_linear_table.vibrational_quantum_1_upper,
                    structures_non_linear_table.vibrational_quantum_2_upper,
                    structures_non_linear_table.vibrational_quantum_3_upper
            )
        )

    def save(self, structure: models.Structure, commit=True):
        """
        Method to save a structure object to the database
        :param structure: Structure object
        :param commit: Boolean indicating of a commit should be done immediately after the save operation. If this is
        set to False, the commit method has to be called to commit the changes.
        """
        cursor = self._database.cursor()
        cursor.execute(queries.structure_insert_query, (
            structure.molecule_id, structure.isotopologue_order_num, structure.global_quanta_upper.hitran_raw_quanta,
            structure.global_quanta_lower.hitran_raw_quanta, structure.max_line_strength_296K,
        ))

        # if the structure is for a diatomic molecule, additional parameters have to be stored
        if models.Structure.get_global_quanta_model_type_by_molecule_id(structure.molecule_id) == models.TriatomicGlobalQuanta:
            cursor.execute(queries.structure_diatomic_insert_query, (
                cursor.lastrowid,
                structure.global_quanta_lower.vibrational_quantum_1,
                structure.global_quanta_lower.vibrational_quantum_2,
                structure.global_quanta_lower.vibrational_quantum_3,
                structure.global_quanta_upper.vibrational_quantum_1,
                structure.global_quanta_upper.vibrational_quantum_2,
                structure.global_quanta_upper.vibrational_quantum_3,
            ))

        if commit:
            self._database.commit()

        structure.structure_id = cursor.lastrowid

        return structure.structure_id


class StructurePyramidalTetratomicMoleculeRepository(StructureLinearMoleculeRepository):

    def __init__(self, database: sqlite3.Connection):
        structures_table, structures_pyramidal_tetratomic_table = Tables('structures', 'structures_pyramidal_tetratomic')
        self._database = database
        self._structure_query = (
            Query().from_(structures_table)
                .join(structures_pyramidal_tetratomic_table, how=JoinType.left_outer)
                .on(structures_table.structure_id == structures_pyramidal_tetratomic_table.structure_id)
                .select(
                structures_table.structure_id,
                structures_table.molecule_id,
                structures_table.isotopologue_order_num,
                structures_table.global_quanta_lower,
                structures_table.global_quanta_upper,
                structures_table.max_line_strength_296K,
                structures_pyramidal_tetratomic_table.vibrational_quantum_1_lower,
                structures_pyramidal_tetratomic_table.vibrational_quantum_2_lower,
                structures_pyramidal_tetratomic_table.vibrational_quantum_3_lower,
                structures_pyramidal_tetratomic_table.vibrational_quantum_4_lower,
                structures_pyramidal_tetratomic_table.vibrational_angular_momentum_3_lower,
                structures_pyramidal_tetratomic_table.vibrational_angular_momentum_4_lower,
                structures_pyramidal_tetratomic_table.vibrational_angular_momentum_lower,
                structures_pyramidal_tetratomic_table.vibrational_symmetry_lower,
                structures_pyramidal_tetratomic_table.vibrational_quantum_1_upper,
                structures_pyramidal_tetratomic_table.vibrational_quantum_2_upper,
                structures_pyramidal_tetratomic_table.vibrational_quantum_3_upper,
                structures_pyramidal_tetratomic_table.vibrational_quantum_4_upper,
                structures_pyramidal_tetratomic_table.vibrational_angular_momentum_3_upper,
                structures_pyramidal_tetratomic_table.vibrational_angular_momentum_4_upper,
                structures_pyramidal_tetratomic_table.vibrational_angular_momentum_upper,
                structures_pyramidal_tetratomic_table.vibrational_symmetry_upper
            )
        )

    def save(self, structure: models.Structure, commit=True):
        """
        Method to save a structure object to the database
        :param structure: Structure object
        :param commit: Boolean indicating of a commit should be done immediately after the save operation. If this is
        set to False, the commit method has to be called to commit the changes.
        """
        cursor = self._database.cursor()
        cursor.execute(queries.structure_insert_query, (
            structure.molecule_id, structure.isotopologue_order_num, structure.global_quanta_upper.hitran_raw_quanta,
            structure.global_quanta_lower.hitran_raw_quanta, structure.max_line_strength_296K,
        ))

        # if the structure is for a diatomic molecule, additional parameters have to be stored
        if models.Structure.get_global_quanta_model_type_by_molecule_id(
                structure.molecule_id) == models.PyramidalTetratomicGlobalQuanta:
            cursor.execute(queries.structure_pyramidal_tetratomic_insert_query, (
                cursor.lastrowid,
                structure.global_quanta_lower.vibrational_quantum_1,
                structure.global_quanta_lower.vibrational_quantum_2,
                structure.global_quanta_lower.vibrational_quantum_3,
                structure.global_quanta_lower.vibrational_quantum_4,
                structure.global_quanta_lower.vibrational_angular_momentum_3,
                structure.global_quanta_lower.vibrational_angular_momentum_4,
                structure.global_quanta_lower.vibrational_angular_momentum,
                structure.global_quanta_lower.vibrational_symmetry,
                structure.global_quanta_upper.vibrational_quantum_1,
                structure.global_quanta_upper.vibrational_quantum_2,
                structure.global_quanta_upper.vibrational_quantum_3,
                structure.global_quanta_upper.vibrational_quantum_4,
                structure.global_quanta_upper.vibrational_angular_momentum_3,
                structure.global_quanta_upper.vibrational_angular_momentum_4,
                structure.global_quanta_upper.vibrational_angular_momentum,
                structure.global_quanta_upper.vibrational_symmetry
            ))
        if commit:
            self._database.commit()

        structure.structure_id = cursor.lastrowid

        return structure.structure_id


class StructurePentatomicMoleculeRepository(StructureLinearMoleculeRepository):

    def __init__(self, database: sqlite3.Connection):
        structures_table, structures = Tables('structures', 'structures_pentatomic')
        self._database = database
        self._structure_query = (
            Query().from_(structures_table)
                .join(structures, how=JoinType.left_outer)
                .on(structures_table.structure_id == structures.structure_id)
                .select(
                structures_table.structure_id,
                structures_table.molecule_id,
                structures_table.isotopologue_order_num,
                structures_table.global_quanta_lower,
                structures_table.global_quanta_upper,
                structures_table.max_line_strength_296K,
                structures.vibrational_quantum_1_lower,
                structures.vibrational_quantum_2_lower,
                structures.vibrational_quantum_3_lower,
                structures.vibrational_quantum_4_lower,
                structures.vibrational_quantum_5_lower,
                structures.vibrational_quantum_6_lower,
                structures.vibrational_symmetry_lower,
                structures.vibrational_quantum_multiplicity_index_lower,
                structures.vibrational_quantum_1_upper,
                structures.vibrational_quantum_2_upper,
                structures.vibrational_quantum_3_upper,
                structures.vibrational_quantum_4_upper,
                structures.vibrational_quantum_5_upper,
                structures.vibrational_quantum_6_upper,
                structures.vibrational_symmetry_upper,
                structures.vibrational_quantum_multiplicity_index_upper
            )
        )

    def save(self, structure: models.Structure, commit=True):
        """
        Method to save a structure object to the database
        :param structure: Structure object
        :param commit: Boolean indicating of a commit should be done immediately after the save operation. If this is
        set to False, the commit method has to be called to commit the changes.
        """
        cursor = self._database.cursor()
        cursor.execute(queries.structure_insert_query, (
            structure.molecule_id, structure.isotopologue_order_num, structure.global_quanta_upper.hitran_raw_quanta,
            structure.global_quanta_lower.hitran_raw_quanta, structure.max_line_strength_296K,
        ))

        # if the structure is for a diatomic molecule, additional parameters have to be stored
        if models.Structure.get_global_quanta_model_type_by_molecule_id(
                structure.molecule_id) == models.PentatomicGlobalQuanta:
            cursor.execute(queries.structure_pentatomic_insert_query, (
                cursor.lastrowid,
                structure.vibrational_quantum_1_lower,
                structure.vibrational_quantum_2_lower,
                structure.vibrational_quantum_3_lower,
                structure.vibrational_quantum_4_lower,
                structure.vibrational_quantum_5_lower,
                structure.vibrational_quantum_6_lower,
                structure.vibrational_symmetry_lower,
                structure.vibrational_quantum_multiplicity_index_lower,
                structure.vibrational_quantum_1_upper,
                structure.vibrational_quantum_2_upper,
                structure.vibrational_quantum_3_upper,
                structure.vibrational_quantum_4_upper,
                structure.vibrational_quantum_5_upper,
                structure.vibrational_quantum_6_upper,
                structure.vibrational_symmetry_upper,
                structure.vibrational_quantum_multiplicity_index_upper
            ))
        if commit:
            self._database.commit()

        structure.structure_id = cursor.lastrowid

        return structure.structure_id


class LineLinearMoleculeRepository(LineRepositoryInterface):
    """
    Object to abstract the database communication from the rest of the code. This object provides methods to interact
    with absorption line records in the database.
    """
    _database = None

    def __init__(self, database: sqlite3.Connection):
        lines_table, lines_diatomic_or_linear_table, structures_table = \
            Tables('lines', 'lines_diatomic_or_linear','structures')

        self._database = database
        self._line_query = Query().from_(lines_table) \
            .join(structures_table, how=JoinType.inner) \
            .on(structures_table.structure_id == lines_table.structure_id) \
            .join(lines_diatomic_or_linear_table, how=JoinType.left_outer) \
            .on(lines_table.line_id == lines_diatomic_or_linear_table.line_id) \
            .select(
                structures_table.molecule_id,
                lines_table.line_id,
                lines_table.structure_id,
                lines_table.wavenumber_vacuum,
                lines_table.line_strength_296K,
                lines_table.einstein_A,
                lines_table.broadening_hw_air,
                lines_table.broadening_hw_self,
                lines_table.energy_lower_state,
                lines_table.broadening_temp_coefficient,
                lines_table.pressure_shift,
                lines_table.local_quanta_lower,
                lines_table.local_quanta_upper,
                lines_table.rotational_statistical_weight_lower,
                lines_table.rotational_statistical_weight_upper,
                lines_diatomic_or_linear_table.J_lower,
                lines_diatomic_or_linear_table.symmetry_lower,
                lines_diatomic_or_linear_table.F_lower,
                lines_diatomic_or_linear_table.J_upper,
                lines_diatomic_or_linear_table.symmetry_upper,
                lines_diatomic_or_linear_table.F_upper,
            )

    def save(self, line: models.Line, molecule_id, commit=True):
        """
        Method to save an absorption line object to the database
        :param line: Line object
        :param molecule_id: The id of the molecule to which the line belongs. This is important to find out if an additional
        record needs to be stored for the processed local quanta. This can be done for certain molecule ids, while this
        is not done for other molecules.
        :param commit: Boolean indicating of a commit should be done immediately after the save operation. If this is
        set to False, the commit method has to be called to commit the changes.
        """
        cursor = self._database.cursor()
        cursor.execute(queries.line_insert_query, (
            line.structure_id, line.wavenumber_vacuum, line.line_strength_296K, line.einstein_A, line.broadening_hw_air,
            line.broadening_hw_self, line.energy_lower_state, line.broadening_temp_coefficient, line.pressure_shift,
            line.local_quanta_lower.hitran_raw_quanta, line.local_quanta_upper.hitran_raw_quanta,
            line.rotational_statistical_weight_lower, line.rotational_statistical_weight_upper,
        ))

        # if the structure is for a diatomic or linear molecule, additional parameters have to be stored
        if models.Line.get_local_quanta_model_type_by_molecule_id(molecule_id) == models.DiatomicOrLinearLocalQuanta:
            cursor.execute(queries.line_diatomic_or_linear_insert_query, (
                cursor.lastrowid,
                line.local_quanta_lower.J,
                line.local_quanta_lower.symmetry,
                line.local_quanta_lower.F,
                line.local_quanta_upper.J,
                line.local_quanta_upper.symmetry,
                line.local_quanta_upper.F,
            ))
        elif models.Line.get_local_quanta_model_type_by_molecule_id(molecule_id) == models.PyramidalTetratomicLocalQuanta:
            cursor.execute(queries.line_pyramidal_tetratomic_insert_query, (
                cursor.lastrowid,
                line.local_quanta_lower.J,
                line.local_quanta_lower.K,
                line.local_quanta_lower.inversion_symmetry,
                line.local_quanta_lower.symmetry_rotational,
                line.local_quanta_lower.symmetry_total,
                line.local_quanta_upper.J,
                line.local_quanta_upper.K,
                line.local_quanta_upper.inversion_symmetry,
                line.local_quanta_upper.symmetry_rotational,
                line.local_quanta_upper.symmetry_total
            ))

        if commit:
            self._database.commit()

        line.line_id = cursor.lastrowid

        return line.line_id

    def get_by_id(self, line_id):
        """
        Method to retrieve a line object from the database with a given line id
        :param line_id: The id of the line to retrieve
        :return: 
        """
        query_filter = filters.LineIdLineQueryFilter(line_id)
        result = self.get_filtered(query_filter)

        try:
            return result.__next__()
        except StopIteration:
            raise IndexError

    def get_by_structure_id(self, structure_id):
        """
        Method to retrieve all lines for a given (vibrational) structure
        :param structure_id: The id of the structure to which the line belongs
        :return: 
        """
        query_filter = filters.StructureIdLineQueryFilter(structure_id)
        return self.get_filtered(query_filter)

    def get_filtered(self, query_filter: filters.QueryFilterInterface):
        """
        Method to retrieve lines with a filter applied to the results.
        :param query_filter: A filter object that filters the results
        :return: 
        """
        query = query_filter.filter_query(self._line_query)
        cursor = self._database.cursor()
        cursor.execute(str(query))

        return results.LineResultSet(cursor)

    def commit(self):
        """
        Write pending changes to the database file. Typically, this is done automatically when calling save(), but this
        can be disabled by passing commit=False to the save method. By not committing after each save, database
        throughput can be maximized.
        :return: 
        """
        self._database.commit()


class LineNonLinearMoleculeRepository(LineLinearMoleculeRepository):
    def __init__(self, database: sqlite3.Connection):
        lines_table, lines_triatomic_non_linear_table, structures_table = \
            Tables('lines', 'lines_triatomic_non_linear', 'structures')

        self._database = database
        self._line_query = Query().from_(lines_table) \
            .join(structures_table, how=JoinType.inner) \
            .on(structures_table.structure_id == lines_table.structure_id) \
            .join(lines_triatomic_non_linear_table, how=JoinType.left_outer) \
            .on(lines_table.line_id == lines_triatomic_non_linear_table.line_id) \
            .select(
            structures_table.molecule_id,
            lines_table.line_id,
            lines_table.structure_id,
            lines_table.wavenumber_vacuum,
            lines_table.line_strength_296K,
            lines_table.einstein_A,
            lines_table.broadening_hw_air,
            lines_table.broadening_hw_self,
            lines_table.energy_lower_state,
            lines_table.broadening_temp_coefficient,
            lines_table.pressure_shift,
            lines_table.local_quanta_lower,
            lines_table.local_quanta_upper,
            lines_table.rotational_statistical_weight_lower,
            lines_table.rotational_statistical_weight_upper,
            lines_triatomic_non_linear_table.J_lower,
            lines_triatomic_non_linear_table.K_lower,
            lines_triatomic_non_linear_table.J_upper,
            lines_triatomic_non_linear_table.K_upper
        )

    def save(self, line: models.Line, molecule_id, commit=True):
        """
        Method to save an absorption line object to the database
        :param line: Line object
        :param molecule_id: The id of the molecule to which the line belongs. This is important to find out if an additional
        record needs to be stored for the processed local quanta. This can be done for certain molecule ids, while this
        is not done for other molecules.
        :param commit: Boolean indicating of a commit should be done immediately after the save operation. If this is
        set to False, the commit method has to be called to commit the changes.
        """
        cursor = self._database.cursor()
        cursor.execute(queries.line_insert_query, (
            line.structure_id, line.wavenumber_vacuum, line.line_strength_296K, line.einstein_A, line.broadening_hw_air,
            line.broadening_hw_self, line.energy_lower_state, line.broadening_temp_coefficient, line.pressure_shift,
            line.local_quanta_lower.hitran_raw_quanta, line.local_quanta_upper.hitran_raw_quanta,
            line.rotational_statistical_weight_lower, line.rotational_statistical_weight_upper,
        ))

        if models.Line.get_local_quanta_model_type_by_molecule_id(
                molecule_id) == models.TriatomicNonLinearLocalQuanta:
            cursor.execute(queries.line_triatomic_non_linear_insert_query, (
                cursor.lastrowid,
                line.local_quanta_lower.J,
                line.local_quanta_lower.K,
                line.local_quanta_upper.J,
                line.local_quanta_upper.K
            ))

        if commit:
            self._database.commit()

        line.line_id = cursor.lastrowid

        return line.line_id


class LinePyramidalTetratomicMoleculeRepository(LineLinearMoleculeRepository):
    def __init__(self, database: sqlite3.Connection):
        lines_table, lines_pyramidal_tetratomic_table, structures_table = \
            Tables('lines', 'lines_pyramidal_tetratomic', 'structures')

        self._database = database
        self._line_query = Query().from_(lines_table) \
            .join(structures_table, how=JoinType.inner) \
            .on(structures_table.structure_id == lines_table.structure_id) \
            .join(lines_pyramidal_tetratomic_table, how=JoinType.left_outer) \
            .on(lines_table.line_id == lines_pyramidal_tetratomic_table.line_id) \
            .select(
            structures_table.molecule_id,
            lines_table.line_id,
            lines_table.structure_id,
            lines_table.wavenumber_vacuum,
            lines_table.line_strength_296K,
            lines_table.einstein_A,
            lines_table.broadening_hw_air,
            lines_table.broadening_hw_self,
            lines_table.energy_lower_state,
            lines_table.broadening_temp_coefficient,
            lines_table.pressure_shift,
            lines_table.local_quanta_lower,
            lines_table.local_quanta_upper,
            lines_table.rotational_statistical_weight_lower,
            lines_table.rotational_statistical_weight_upper,
            lines_pyramidal_tetratomic_table.J_lower,
            lines_pyramidal_tetratomic_table.K_lower,
            lines_pyramidal_tetratomic_table.inversion_symmetry_lower,
            lines_pyramidal_tetratomic_table.symmetry_rotational_lower,
            lines_pyramidal_tetratomic_table.symmetry_total_lower,
            lines_pyramidal_tetratomic_table.J_upper,
            lines_pyramidal_tetratomic_table.K_upper,
            lines_pyramidal_tetratomic_table.inversion_symmetry_upper,
            lines_pyramidal_tetratomic_table.symmetry_rotational_upper,
            lines_pyramidal_tetratomic_table.symmetry_total_upper,
        )

    def save(self, line: models.Line, molecule_id, commit=True):
        """
        Method to save an absorption line object to the database
        :param line: Line object
        :param molecule_id: The id of the molecule to which the line belongs. This is important to find out if an additional
        record needs to be stored for the processed local quanta. This can be done for certain molecule ids, while this
        is not done for other molecules.
        :param commit: Boolean indicating of a commit should be done immediately after the save operation. If this is
        set to False, the commit method has to be called to commit the changes.
        """
        cursor = self._database.cursor()
        cursor.execute(queries.line_insert_query, (
            line.structure_id, line.wavenumber_vacuum, line.line_strength_296K, line.einstein_A, line.broadening_hw_air,
            line.broadening_hw_self, line.energy_lower_state, line.broadening_temp_coefficient, line.pressure_shift,
            line.local_quanta_lower.hitran_raw_quanta, line.local_quanta_upper.hitran_raw_quanta,
            line.rotational_statistical_weight_lower, line.rotational_statistical_weight_upper,
        ))

        if models.Line.get_local_quanta_model_type_by_molecule_id(
                molecule_id) == models.PyramidalTetratomicLocalQuanta:
            cursor.execute(queries.line_pyramidal_tetratomic_insert_query, (
                cursor.lastrowid,
                line.local_quanta_lower.J,
                line.local_quanta_lower.K,
                line.local_quanta_lower.inversion_symmetry,
                line.local_quanta_lower.symmetry_rotational,
                line.local_quanta_lower.symmetry_total,
                line.local_quanta_upper.J,
                line.local_quanta_upper.K,
                line.local_quanta_upper.inversion_symmetry,
                line.local_quanta_upper.symmetry_rotational,
                line.local_quanta_upper.symmetry_total
            ))

        if commit:
            self._database.commit()

        line.line_id = cursor.lastrowid

        return line.line_id


class PartitionSumRepository(PartitionSumRepositoryInterface):
    """
    Object to abstract the database communication from the rest of the code. This object provides methods to query
    partition sum records in the database.
    """
    _database = None

    def __init__(self, database: sqlite3.Connection):
        self._database = database
        self._partition_sum_table = Table('partition_sums')

    def get_by_isotopologue(self, isotopologue: models.Isotopologue):
        """
        Retrieve partition sum data for the full available temperature range in the database for a given isotopologue
        :param isotopologue: The isotopologue for which the partition sums should be retrieved
        :return: 
        """
        query = Query()\
            .from_(self._partition_sum_table).select('temperature', 'partition_sum')\
            .where(self._partition_sum_table.isotopologue_order_num == isotopologue.isotopologue_order_num)\
            .where(self._partition_sum_table.molecule_id == isotopologue.molecule_id)

        cursor = self._database.cursor()
        cursor.execute(str(query))

        rows = cursor.fetchall()
        temperatures, partition_sums = zip(*rows)

        return results.PartitionSumSet(temperatures, partition_sums)
