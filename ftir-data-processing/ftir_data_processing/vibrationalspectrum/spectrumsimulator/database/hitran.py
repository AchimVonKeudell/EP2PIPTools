import os
import glob
import struct
import sqlite3
import coloredlogs
import logging

from . import queries, repositories, models, readers, exomol

# Create a logger object.
logger = logging.getLogger('spectrum-database-generator')
logger.setLevel(logging.INFO)

coloredlogs.install(level='DEBUG', logger=logger)


class HITRANValueError(ValueError):
    pass


class HitranDatabaseGenerator:
    def __init__(self, database: 'HitranDatabase'):
        """
        Initialize a database generator object to create and fill a sqlite database from HITRAN/HITEMP source files. 
        The database connection is setup up and a cursor object is created for the creation of tables and the insertion
        of records.
        :param database: An sqlite3 connection object referring to an sqlite database.
        """
        self._database = database
        self._cursor = database.cursor()

        self._molecule_repository = repositories.MoleculeRepository(self._database)
        self._isotopologue_repository = repositories.IsotopologueRepository(self._database)
        self._structure_repository = self._get_structure_repository()
        self._line_repository = self._get_line_repository()

    def _get_structure_repository(self) -> repositories.StructureRepositoryInterface:
        raise NotImplementedError

    def _get_line_repository(self) -> repositories.LineRepositoryInterface:
        raise NotImplementedError

    def create_tables(self):
        """
        Create tables for molecules, isotopologues, (vibrational) structures and absorption lines
        :return: 
        """
        raise NotImplementedError

    def process_molecule_line(self, molecule_line):
        """
        Process a line in the molparam.txt file referring to a molecule. Extract the molecule id and name from the line
        and store it in the database.
        :param molecule_line: 
        :return: 
        """
        # save the molecule to the database
        return self.insert_molecule(models.Molecule(
            molecule_id=int(molecule_line[8:-2]),
            molecule_name=molecule_line[0:6].strip()
        ))

    def insert_molecule(self, molecule):
        """
        Insert a molecule into the database
        :param molecule: Molecule object
        """
        return self._molecule_repository.save(molecule, commit=False)

    def process_isotopologue_line(self, molecule_id, isotopologue_order_num, isotopologue_line):
        """
        Process a line belonging to an isotopologue and insert the data into the database
        :param molecule_id: The database of the current molecule
        :param isotopologue_order_num: The order number of the isotopologue in order of abundance (DESC)
        :param isotopologue_line: The line from the file that specifies the isotopologues parameters
        :return: 
        """
        return self.insert_isotopologue(models.Isotopologue(
            molecule_id=molecule_id,
            isotopologue_order_num=isotopologue_order_num,
            afgl_code=isotopologue_line[0:12].strip(),
            abundance=float(isotopologue_line[12:25]),
            molar_mass=float(isotopologue_line[44:58])
        ))

    def insert_isotopologue(self, isotopologue):
        """
        Insert an isotopologue into the database.
        :param isotopologue: Isotopologue object
        """
        return self._isotopologue_repository.save(isotopologue, commit=False)

    def process_line_line(self, structure_index_list, line_data):
        """
        Process a line in a line parameter file related to lines and structures
        :param structure_index_list: A cache containing all already inserted structures and their unique ids
        :param line_data: The line from the file to be processed
        :return: 
        """

        try:
            # process the global upper and lower quantum into seperate numbers depending on the type of molecule
            global_quanta_model_type = models.Structure.get_global_quanta_model_type_by_molecule_id(line_data.molecule_id)
            global_quanta_lower = self.process_global_quanta(global_quanta_model_type, line_data.global_quanta_lower)
            global_quanta_upper = self.process_global_quanta(global_quanta_model_type, line_data.global_quanta_upper)

            # process the global upper and lower quantum into seperate numbers depending on the type of molecule
            local_quanta_model_type = models.Line.get_local_quanta_model_type_by_molecule_id(line_data.molecule_id)
            local_quanta_lower, local_quanta_upper = self.process_local_quanta(local_quanta_model_type, line_data.local_quanta_lower, line_data.local_quanta_upper)

            # define a string that uniquely defines the structure
            structure_string_id = "%i_%s_%s" % (line_data.isotopologue_order_num, line_data.global_quanta_lower, line_data.global_quanta_upper)

            # check if the structure is already inserted
            if structure_string_id not in structure_index_list:
                # insert the structure
                structure_id = self.insert_structure(models.Structure(
                    molecule_id=line_data.molecule_id,
                    isotopologue_order_num=line_data.isotopologue_order_num,
                    global_quanta_lower=global_quanta_lower,
                    global_quanta_upper=global_quanta_upper
                ))
                structure_index_list[structure_string_id] = structure_id
            else:
                # get the structure id from the structure cache
                structure_id = structure_index_list[structure_string_id]

            # insert the line
            self.insert_line(models.Line(
                structure_id=structure_id,
                wavenumber_vacuum=line_data.wavenumber_vacuum,
                line_strength_296K=line_data.line_strength_296K,
                einstein_A=line_data.einstein_A,
                broadening_hw_air=line_data.broadening_hw_air,
                broadening_hw_self=line_data.broadening_hw_self,
                energy_lower_state=line_data.energy_lower_state,
                broadening_temp_coefficient=line_data.broadening_temp_coefficient,
                pressure_shift=line_data.pressure_shift,
                local_quanta_lower=local_quanta_lower,
                local_quanta_upper=local_quanta_upper,
                rotational_statistical_weight_lower=line_data.rotational_statistical_weight_lower,
                rotational_statistical_weight_upper=line_data.rotational_statistical_weight_upper,
            ), line_data.molecule_id)
        except HITRANValueError:
            return

    def process_global_quanta(self, global_quanta_model_type, hitran_raw_quanta) -> models.GlobalQuanta:
        """
        Extract global quantum numbers from a 15 character string in HITRAN format containing the relevant quantum numbers
        :param global_quanta_model_type:
        :param hitran_raw_quanta:
        :return: 
        """
        raise NotImplementedError

    def process_local_quanta(self, local_quanta_model_type, hitran_raw_quanta_lower, hitran_raw_quanta_upper) -> iter:
        """
        Extract local quantum numbers from two 15 character string in HITRAN format containing the relevant quantum
        numbers for the upper and lower state. Some information needed to reconstruct the upper rotational quantum
        number is stored in the lower local quanta. Therefore, this function needs both quanta as an input.
        :param local_quanta_model_type:
        :param hitran_raw_quanta_lower:
        :param hitran_raw_quanta_upper:
        :return:
        """
        raise NotImplementedError

    def insert_structure(self, structure):
        """
        Insert an structure into the database.
        :param structure: Structure object
        """
        return self._structure_repository.save(structure, commit=False)

    def insert_line(self, line, molecule_id):
        """
        Insert an line into the database.
        :param line: Line object
        :param molecule_id: The id of the molecule to which the line belongs
        """
        return self._line_repository.save(line, molecule_id, commit=False)

    def generate_database(self, hitran_location, hitran_parameter_mapper):
        """
        Generate an sqlite database from HITRAN source files. An existing sqlite database will be overwritten.
        :param hitran_location: 
        :return: 
        """
        self.create_tables()
        self.fill_tables(hitran_location, hitran_parameter_mapper)

    def fill_tables(self, hitran_location, hitran_parameter_mapper):
        """
        Fill the tables with data on molecules, isotopologues, structures and absorption lines
        :param hitran_location: 
        :return: 
        """
        # path of the file containing metadata on molecules and isotopologues
        molecular_parameter_file = r'%s\molparam.txt' % hitran_location

        molecule_counter = 0
        isotopologue_counter = 0

        try:
            with open(molecular_parameter_file) as molecular_parameter_file_handle:
                # skip first line of the file, since it is a header
                next(molecular_parameter_file_handle)

                # read the first line which contains molecular information
                file_line = next(molecular_parameter_file_handle)

                # keep reading lines from the file. the script will break from the loop with a StopIteration of a next
                # line is not available anymore.
                while True:
                    # process the molecule
                    molecule_id = self.process_molecule_line(file_line)

                    # read the line having data on the first isotopologue linked to the current molecule
                    file_line = next(molecular_parameter_file_handle)
                    # start at isotopologue order number 1
                    isotopologue_order_num = 1

                    # keep reading isotopologues until we reach the next molecule
                    while len(file_line) == 60:
                        self.process_isotopologue_line(molecule_id, isotopologue_order_num, file_line)

                        isotopologue_counter += 1
                        isotopologue_order_num += 1
                        file_line = next(molecular_parameter_file_handle)

                    # process the individual lines for a given molecule
                    self.process_molecule_line_data(molecule_id, hitran_location, hitran_parameter_mapper)

                    molecule_counter += 1
        except StopIteration:
            pass
        finally:
            logger.info("%i molecules and %i isotopologues imported" % (molecule_counter, isotopologue_counter))

            self.import_partition_sums(hitran_location)
            logger.info("Partition sums imported")

            self.update_max_line_strength()
            logger.info("Maximum line strengths calculated per structure")

            self._database.commit()
            logger.info("Changes written to database")

    def process_molecule_line_data(self, molecule_id, hitran_location, hitran_parameter_mapper):
        """
        Process all *.par files for a given molecule_id
        :param molecule_id: The id of the molecule for which to process absorption line parameter files
        :param hitran_location: The location of the hitran input files
        :return: 
        """
        # find all absorption line files related to the given molecule_id
        absorption_line_files = glob.glob("%s/%02i_*.par" % (hitran_location, molecule_id))
        # this variable is used to store all already inserted structures with their ids. This is used to link subsequent
        # absorption lines to the correct structures instead of making a new structure for each individual line
        structure_index_list = {}

        # loop through the input files for the given molecule
        for absorption_line_file in absorption_line_files:
            reader = readers.HITRANReader(absorption_line_file, mapper=hitran_parameter_mapper)
            for line_data in reader:
                self.process_line_line(structure_index_list, line_data)

            logger.info("Line parameter file %s imported" % absorption_line_file)

    def import_partition_sums(self, hitran_location):
        """
        Imports partition sums from the parsum.dat for a temperature range between 70K and 3000K and links them to the
        correct isotopologues.
        :param hitran_location: 
        :return: 
        """
        # path of the file containing partition sums for all isotopologues
        partition_sum_file = r'%s/parsum.dat' % hitran_location
        isotopologues_index = {}

        with open(partition_sum_file) as partition_sum_file_handle:
            header_line = next(partition_sum_file_handle)
            header_line_columns = header_line.split()

            column_data = {}

            # loop through all columns except the first one and find the molecule id and isotopologue order num for
            # each column from the molecule name and AFGL code
            for column_index in range(1, len(header_line_columns)):
                molecule_name, afgl_code = header_line_columns[column_index].split("_")

                self._cursor.execute(queries.isotopologue_select_by_afgl, (molecule_name, afgl_code,))
                isotopologue_row = self._cursor.fetchone()

                if isotopologue_row is not None:
                    column_data[column_index] = {'molecule_id': isotopologue_row[0], 'isotopologue_order_num': isotopologue_row[1]}

            # loop through all data and store it in the database
            for partition_sum_line in partition_sum_file_handle:
                partition_sum_line_data = partition_sum_line.split()
                temperature = partition_sum_line_data[0]

                for column_index in range(1, len(partition_sum_line_data)):
                    if column_index in column_data:
                        molecule_id = column_data[column_index]['molecule_id']
                        isotopologue_order_num = column_data[column_index]['isotopologue_order_num']
                        partition_sum = partition_sum_line_data[column_index]

                        self._cursor.execute(queries.partition_sum_insert_query, (molecule_id, isotopologue_order_num, temperature, partition_sum))

    def update_max_line_strength(self):
        """
        Calculates the maximum line strength per structure and stores it in the database
        :return: 
        """
        self._cursor.execute(queries.structure_calculate_statistics_query)


class HitranLinearMoleculesDatabaseGenerator(HitranDatabaseGenerator):

    def _get_structure_repository(self):
        return repositories.StructureLinearMoleculeRepository(self._database)

    def _get_line_repository(self):
        return repositories.LineLinearMoleculeRepository(self._database)

    def create_tables(self):
        """
        Create tables for molecules, isotopologues, (vibrational) structures and absorption lines
        :return:
        """
        self._cursor.executescript(queries.molecule_table_create)
        self._cursor.executescript(queries.isotopologue_table_create)
        self._cursor.executescript(queries.structure_table_create)
        self._cursor.executescript(queries.structure_diatomic_table_create)
        self._cursor.executescript(queries.structure_triatomic_linear_fermi_resonant_table_create)
        self._cursor.executescript(queries.line_table_create)
        self._cursor.executescript(queries.line_diatomic_or_linear_table_create)
        self._cursor.executescript(queries.partition_sums_table_create)

        logger.info("Created tables")

    def process_global_quanta(self, global_quanta_model_type, hitran_raw_quanta):
        """
        Extract global quantum numbers from a 15 character string in HITRAN format containing the relevant quantum numbers
        :param global_quanta_model_type:
        :param hitran_raw_quanta:
        :return:
        """
        _hitran_raw_quanta_stripped = hitran_raw_quanta.strip()
        if _hitran_raw_quanta_stripped == "":
            # print(f'{self.__class__.__name__}: {hitran_raw_quanta=}')
            raise HITRANValueError("Global quanta is empty")

        if global_quanta_model_type == models.DiatomicGlobalQuanta:
            return global_quanta_model_type(
                vibrational_quantum=int(hitran_raw_quanta[13:15])
            )
        elif global_quanta_model_type == models.TriatomicGlobalQuanta:
            return global_quanta_model_type(
                vibrational_quantum_1=to_int(hitran_raw_quanta[9:11]),
                vibrational_quantum_2=to_int(hitran_raw_quanta[11:13]),
                vibrational_quantum_3=to_int(hitran_raw_quanta[13:15])
            )
        elif global_quanta_model_type == models.TriatomicLinearGlobalQuanta:
            def to_int(selected_hitran_raw_quanta):
                # TODO: why is there a ' a', ' b', ' c', or ' d' in v2 for N2O?
                if any([non_int in selected_hitran_raw_quanta for non_int in ['a', 'b', 'c', 'd', 'e']]):
                    return 0
                return int(selected_hitran_raw_quanta)
            return global_quanta_model_type(
                vibrational_quantum_1=to_int(hitran_raw_quanta[7:9]),
                vibrational_quantum_2=to_int(hitran_raw_quanta[9:11]),
                vibrational_quantum_3=to_int(hitran_raw_quanta[13:15]),
                vibrational_angular_momentum=to_int(hitran_raw_quanta[11:13])
            )
        elif global_quanta_model_type == models.TriatomicLinearFermiResonantGlobalQuanta:
            quanta = global_quanta_model_type(
                vibrational_quantum_1=int(hitran_raw_quanta[6:8]),
                vibrational_quantum_2=int(hitran_raw_quanta[8:10]),
                vibrational_quantum_3=int(hitran_raw_quanta[12:14]),
                vibrational_angular_momentum=int(hitran_raw_quanta[10:12]),
                vibrational_ranking_index=int(hitran_raw_quanta[14:15])
            )

            # if the vibrational ranking index is 0, this means either the upper of lower state is unknown in the HITRAN
            # and we should ignore the line, because we cannot reconstruct upper and lower state densities properly
            if quanta.vibrational_ranking_index == 0:
                raise HITRANValueError("Global state unknown, should be ignored.")
            return quanta
        return models.GlobalQuanta(hitran_raw_quanta=hitran_raw_quanta)

    def process_local_quanta(self, local_quanta_model_type, hitran_raw_quanta_lower, hitran_raw_quanta_upper):
        """
        Extract local quantum numbers from two 15 character string in HITRAN format containing the relevant quantum
        numbers for the upper and lower state. Some information needed to reconstruct the upper rotational quantum
        number is stored in the lower local quanta. Therefore, this function needs both quanta as an input.
        :param local_quanta_model_type:
        :param hitran_raw_quanta_lower:
        :param hitran_raw_quanta_upper:
        :return:
        """
        if local_quanta_model_type == models.DiatomicOrLinearLocalQuanta:
            branch = hitran_raw_quanta_lower[5:6]
            J_lower = int(hitran_raw_quanta_lower[6:9])

            if branch == 'P':
                J_upper = J_lower - 1
            elif branch == 'Q':
                J_upper = J_lower
            elif branch == 'R':
                J_upper = J_lower + 1
            else:
                raise Exception("Unknown branch %s found in HITRAN database" % branch)

            local_quanta_lower = local_quanta_model_type(
                J=J_lower,
                symmetry=hitran_raw_quanta_lower[9:10].strip(),
                F=hitran_raw_quanta_lower[10:15].strip()
            )
            local_quanta_upper = local_quanta_model_type(
                J=J_upper,
                symmetry=hitran_raw_quanta_upper[9:10].strip(),
                F=hitran_raw_quanta_upper[10:15].strip()
            )
            return local_quanta_lower, local_quanta_upper
        return models.LocalQuanta(hitran_raw_quanta=hitran_raw_quanta_lower), \
               models.LocalQuanta(hitran_raw_quanta=hitran_raw_quanta_upper)


class HitranTriatomicNonLinearMoleculesDatabaseGenerator(HitranDatabaseGenerator):

    def _get_structure_repository(self):
        return repositories.StructureNonLinearMoleculeRepository(self._database)

    def _get_line_repository(self):
        return repositories.LineNonLinearMoleculeRepository(self._database)

    def create_tables(self):
        """
        Create tables for molecules, isotopologues, (vibrational) structures and absorption lines
        :return:
        """
        self._cursor.executescript(queries.molecule_table_create)
        self._cursor.executescript(queries.isotopologue_table_create)
        self._cursor.executescript(queries.structure_table_create)
        self._cursor.executescript(queries.structure_triatomic_non_linear_table_create)
        self._cursor.executescript(queries.line_table_create)
        self._cursor.executescript(queries.line_triatomic_non_linear_table_create)
        self._cursor.executescript(queries.partition_sums_table_create)

        logger.info("Created tables")

    def process_global_quanta(self, global_quanta_model_type, hitran_raw_quanta):
        """
        Extract global quantum numbers from a 15 character string in HITRAN format containing the relevant quantum numbers
        :param global_quanta_model_type:
        :param hitran_raw_quanta:
        :return:
        """
        _hitran_raw_quanta_stripped = hitran_raw_quanta.strip()
        if _hitran_raw_quanta_stripped == "":
            # print(f'{self.__class__.__name__}: {hitran_raw_quanta=}')
            raise HITRANValueError("Global quanta is empty")

        if global_quanta_model_type == models.TriatomicGlobalQuanta:
            return global_quanta_model_type(
                vibrational_quantum_1=to_int(hitran_raw_quanta[9:11]),
                vibrational_quantum_2=to_int(hitran_raw_quanta[11:13]),
                vibrational_quantum_3=to_int(hitran_raw_quanta[13:15])
            )
        return models.GlobalQuanta(hitran_raw_quanta=hitran_raw_quanta)

    def process_local_quanta(self, local_quanta_model_type, hitran_raw_quanta_lower, hitran_raw_quanta_upper):
        """
        Extract local quantum numbers from two 15 character string in HITRAN format containing the relevant quantum
        numbers for the upper and lower state. Some information needed to reconstruct the upper rotational quantum
        number is stored in the lower local quanta. Therefore, this function needs both quanta as an input.
        :param local_quanta_model_type:
        :param hitran_raw_quanta_lower:
        :param hitran_raw_quanta_upper:
        :return:
        """
        if local_quanta_model_type == models.TriatomicNonLinearLocalQuanta:
            J_lower = int(hitran_raw_quanta_lower[0:3])
            K_lower = int(hitran_raw_quanta_lower[3:6])

            J_upper = int(hitran_raw_quanta_upper[:3])
            K_upper = int(hitran_raw_quanta_upper[3:6])

            return local_quanta_model_type(J=J_lower, K=K_lower), \
                   local_quanta_model_type(J=J_upper, K=K_upper)
        return models.LocalQuanta(hitran_raw_quanta=hitran_raw_quanta_lower), \
               models.LocalQuanta(hitran_raw_quanta=hitran_raw_quanta_upper)


class HitranPyramidalTetratomicDatabaseGenerator(HitranDatabaseGenerator):

    def _get_structure_repository(self):
        return repositories.StructurePyramidalTetratomicMoleculeRepository(self._database)

    def _get_line_repository(self):
        return repositories.LinePyramidalTetratomicMoleculeRepository(self._database)

    def create_tables(self):
        """
        Create tables for molecules, isotopologues, (vibrational) structures and absorption lines
        :return:
        """
        self._cursor.executescript(queries.molecule_table_create)
        self._cursor.executescript(queries.isotopologue_table_create)
        self._cursor.executescript(queries.structure_table_create)
        self._cursor.executescript(queries.structure_pyramidal_tetratomic_table_create)
        self._cursor.executescript(queries.line_table_create)
        self._cursor.executescript(queries.line_pyramidal_tetratomic_table_create)
        self._cursor.executescript(queries.partition_sums_table_create)

        logger.info("Created tables")

    def process_global_quanta(self, global_quanta_model_type, hitran_raw_quanta):
        """
        Extract global quantum numbers from a 15 character string in HITRAN format containing the relevant quantum numbers
        :param global_quanta_model_type:
        :param hitran_raw_quanta:
        :return:
        """
        _hitran_raw_quanta_stripped = hitran_raw_quanta.strip()
        if _hitran_raw_quanta_stripped == "":
            raise HITRANValueError("Global quanta is empty")

        if global_quanta_model_type == models.PyramidalTetratomicGlobalQuanta:
            _hitran_raw_quanta_split = _hitran_raw_quanta_stripped.split()
            if len(_hitran_raw_quanta_split) == 4:          # for 14NH3, e.g. global quanta = " 1001 01 1 E'   "
                return global_quanta_model_type(
                    vibrational_quantum_1=int(hitran_raw_quanta[1]),
                    vibrational_quantum_2=int(hitran_raw_quanta[2]),
                    vibrational_quantum_3=int(hitran_raw_quanta[3]),
                    vibrational_quantum_4=int(hitran_raw_quanta[4]),
                    vibrational_angular_momentum_3=int(hitran_raw_quanta[6]),
                    vibrational_angular_momentum_4=int(hitran_raw_quanta[7]),
                    vibrational_angular_momentum=int(hitran_raw_quanta[9]),
                    vibrational_symmetry=hitran_raw_quanta[10:].strip()
                )
            elif len(_hitran_raw_quanta_split) == 5:        # for 15NH3, e.g. global quanta = "      0 0 1 1 s "
                return global_quanta_model_type(
                    vibrational_quantum_1=int(_hitran_raw_quanta_split[0]),
                    vibrational_quantum_2=int(_hitran_raw_quanta_split[1]),
                    vibrational_quantum_3=int(_hitran_raw_quanta_split[2]),
                    vibrational_quantum_4=int(_hitran_raw_quanta_split[3]),
                    vibrational_angular_momentum_3=None,
                    vibrational_angular_momentum_4=None,
                    vibrational_angular_momentum=None,
                    vibrational_symmetry=_hitran_raw_quanta_split[4]
                )
            else:
                raise HITRANValueError("Global state unknwon, should be ignored")
        return models.GlobalQuanta(hitran_raw_quanta=hitran_raw_quanta)

    def process_local_quanta(self, local_quanta_model_type, hitran_raw_quanta_lower, hitran_raw_quanta_upper):
        """
        Extract local quantum numbers from two 15 character string in HITRAN format containing the relevant quantum
        numbers for the upper and lower state. Some information needed to reconstruct the upper rotational quantum
        number is stored in the lower local quanta. Therefore, this function needs both quanta as an input.
        :param local_quanta_model_type:
        :param hitran_raw_quanta_lower:
        :param hitran_raw_quanta_upper:
        :return:
        """
        if local_quanta_model_type == models.PyramidalTetratomicLocalQuanta:
            if len(hitran_raw_quanta_lower.split()) == 2:
                return [
                    local_quanta_model_type(
                        J=int(hitran_raw_quanta[:3]),
                        K=int(hitran_raw_quanta[3:6]),
                        inversion_symmetry=hitran_raw_quanta[6:8].strip(),
                        symmetry_rotational='',
                        symmetry_total=''
                    ) for hitran_raw_quanta in (hitran_raw_quanta_lower, hitran_raw_quanta_upper)
                ]
            return [
                local_quanta_model_type(
                    J=int(hitran_raw_quanta[:2]),
                    K=int(hitran_raw_quanta[2:5]),
                    inversion_symmetry=hitran_raw_quanta[5:7].strip(),
                    symmetry_rotational=hitran_raw_quanta[8:11].strip(),
                    symmetry_total=hitran_raw_quanta[11:14].strip()
                ) for hitran_raw_quanta in (hitran_raw_quanta_lower, hitran_raw_quanta_upper)
            ]
        return models.LocalQuanta(hitran_raw_quanta=hitran_raw_quanta_lower), \
               models.LocalQuanta(hitran_raw_quanta=hitran_raw_quanta_upper)


class HitranDatabase(sqlite3.Connection):
    def __init__(self, database_path):
        self.database_path = database_path

        super().__init__(database_path)

    @staticmethod
    def generate_database_linear_molecules(database_path, hitran_directory,
                                           hitran_parameter_mapper=readers.HITRAN160ParameterMapper()):
        """
        Method to generate a new sqlite database from HITRAN source files. This returns the newly generated sqlite
        database.

        :param database_path: The path to the sqlite HITRAN database file
        :param hitran_directory:
        :param hitran_parameter_mapper:
        :return: object referring to the newly generated database
        """
        if not os.path.isdir(hitran_directory):
            raise Exception("Hitran directory '%s' does not exist." % hitran_directory)

        logger.info("Starting import of HITRAN files from '%s'" % hitran_directory)

        database = HitranDatabase(database_path)
        generator = HitranLinearMoleculesDatabaseGenerator(database)
        generator.generate_database(hitran_location=hitran_directory, hitran_parameter_mapper=hitran_parameter_mapper)
        return database

    @staticmethod
    def generate_database_pyramidal_tetratomic_molecules(database_path, hitran_directory,
                                                         hitran_parameter_mapper=readers.HITRAN160ParameterMapper()):
        """
        Method to generate a new sqlite database from HITRAN source files. This returns the newly generated sqlite
        database.

        :param database_path: The path to the sqlite HITRAN database file
        :param hitran_directory:
        param hitran_parameter_mapper:
        :return: hitran_directory: object referring to the newly generated database
        """
        if not os.path.isdir(hitran_directory):
            raise Exception("Hitran directory '%s' does not exist." % hitran_directory)

        logger.info("Starting import of HITRAN files from '%s'" % hitran_directory)

        database = HitranDatabase(database_path)
        generator = HitranPyramidalTetratomicDatabaseGenerator(database)
        generator.generate_database(hitran_location=hitran_directory, hitran_parameter_mapper=hitran_parameter_mapper)
        return database


    @staticmethod
    def generate_database_non_linear_molecules(database_path, hitran_directory,
                                               hitran_parameter_mapper=readers.HITRAN160ParameterMapper()):
        if not os.path.isdir(hitran_directory):
            raise Exception("Hitran directory '%s' does not exist." % hitran_directory)

        logger.info("Starting import of HITRAN files from '%s'" % hitran_directory)

        database = HitranDatabase(database_path)
        generator = HitranTriatomicNonLinearMoleculesDatabaseGenerator(database)
        generator.generate_database(hitran_location=hitran_directory, hitran_parameter_mapper=hitran_parameter_mapper)
        return database