import sqlite3
import numpy as np
from typing import List
from . import models


def record_to_molecule_model(record):
    """
    Takes a database record with all molecule fields selected and converts it to a molecule model object
    :param record: Array containing the values of all molecule database fields for a given row
    :return: 
    """
    return models.Molecule(
        molecule_id=record[0],
        molecule_name=record[1]
    )


def record_to_isotopologue_model(record):
    """
    Takes a database record with all isotopologue fields selected and converts it to an isotopologue model object
    :param record: Array containing the values of all isotopologue database fields for a given row
    :return: 
    """
    return models.Isotopologue(
        molecule_id=record[0],
        isotopologue_order_num=record[1],
        afgl_code=record[2],
        abundance=record[3],
        molar_mass=record[4]
    )


def record_to_structure_model(record):
    """
    Takes a database record with all structure fields selected and converts it to an structure model object
    :param record: Array containing the values of all structure database fields for a given row
    :return: 
    """
    return models.Structure(
        structure_id=record[0],
        molecule_id=record[1],
        isotopologue_order_num=record[2],
        global_quanta_lower=models.GlobalQuanta(hitran_raw_quanta=record[3]),
        global_quanta_upper=models.GlobalQuanta(hitran_raw_quanta=record[4]),
        max_line_strength_296K=record[5]
    )


def record_to_structure_model_diatomic(record):
    """
    Takes a database record with all structure fields selected and converts it to an structure model object
    :param record: Array containing the values of all stucture database fields for a given row
    :return:
    """
    return models.Structure(
        structure_id=record[0],
        molecule_id=record[1],
        isotopologue_order_num=record[2],
        global_quanta_lower=models.DiatomicGlobalQuanta(vibrational_quantum=record[6]),
        global_quanta_upper=models.DiatomicGlobalQuanta(vibrational_quantum=record[7]),
        max_line_strength_296K=record[5]
    )


def record_to_structure_model_triatomic_linear_fermi_resonant(record):
    """
    Takes a database record with all structure fields selected and converts it to an structure model object
    :param record: Array containing the values of all stucture database fields for a given row
    :return:
    """
    return models.Structure(
        structure_id=record[0],
        molecule_id=record[1],
        isotopologue_order_num=record[2],
        global_quanta_lower=models.TriatomicLinearFermiResonantGlobalQuanta(
            vibrational_quantum_1=record[8],
            vibrational_quantum_2=record[9],
            vibrational_quantum_3=record[10],
            vibrational_angular_momentum=record[11],
            vibrational_ranking_index=record[12]
        ),
        global_quanta_upper=models.TriatomicLinearFermiResonantGlobalQuanta(
            vibrational_quantum_1=record[13],
            vibrational_quantum_2=record[14],
            vibrational_quantum_3=record[15],
            vibrational_angular_momentum=record[16],
            vibrational_ranking_index=record[17]
        ),
        max_line_strength_296K=record[5]
    )


def record_to_structure_triatomic(record):
    return models.Structure(
        structure_id=record[0],
        molecule_id=record[1],
        isotopologue_order_num=record[2],
        global_quanta_lower=models.TriatomicGlobalQuanta(
            vibrational_quantum_1=record[8],
            vibrational_quantum_2=record[9],
            vibrational_quantum_3=record[10],
        ),
        global_quanta_upper=models.TriatomicGlobalQuanta(
            vibrational_quantum_1=record[11],
            vibrational_quantum_2=record[12],
            vibrational_quantum_3=record[13],
        ),
        max_line_strength_296K=record[5]
    )

def record_to_structure_model_pyramidal_tetratomic(record):
    return models.Structure(
        structure_id=record[0],
        molecule_id=record[1],
        isotopologue_order_num=record[2],
        global_quanta_lower=models.PyramidalTetratomicGlobalQuanta(
            vibrational_quantum_1=record[6],
            vibrational_quantum_2=record[7],
            vibrational_quantum_3=record[8],
            vibrational_quantum_4=record[9],
            vibrational_angular_momentum_3=record[10],
            vibrational_angular_momentum_4=record[11],
            vibrational_angular_momentum=record[12],
            vibrational_symmetry=record[13]
        ),
        global_quanta_upper=models.PyramidalTetratomicGlobalQuanta(
            vibrational_quantum_1=record[14],
            vibrational_quantum_2=record[15],
            vibrational_quantum_3=record[16],
            vibrational_quantum_4=record[17],
            vibrational_angular_momentum_3=record[18],
            vibrational_angular_momentum_4=record[19],
            vibrational_angular_momentum=record[20],
            vibrational_symmetry=record[21]
        ),
        max_line_strength_296K=record[5]
    )


def record_to_line_model(record):
    """
    Takes a database record with all line fields selected and converts it to an line model object
    :param record: Array containing the values of all line database fields for a given row
    :return: 
    """
    return models.Line(
        line_id=record[1],
        structure_id=record[2],
        wavenumber_vacuum=record[3],
        line_strength_296K=record[4],
        einstein_A=record[5],
        broadening_hw_air=record[6],
        broadening_hw_self=record[7],
        energy_lower_state=record[8],
        broadening_temp_coefficient=record[9],
        pressure_shift=record[10],
        local_quanta_lower=models.LocalQuanta(hitran_raw_quanta=record[11]),
        local_quanta_upper=models.LocalQuanta(hitran_raw_quanta=record[12]),
        rotational_statistical_weight_lower=record[13],
        rotational_statistical_weight_upper=record[14],
    )


def record_to_line_model_diatomic_or_linear(record):
    """
    Takes a database record with all line fields selected and converts it to an line model object
    :param record: Array containing the values of all line database fields for a given row
    :return:
    """
    return models.Line(
        line_id=record[1],
        structure_id=record[2],
        wavenumber_vacuum=record[3],
        line_strength_296K=record[4],
        einstein_A=record[5],
        broadening_hw_air=record[6],
        broadening_hw_self=record[7],
        energy_lower_state=record[8],
        broadening_temp_coefficient=record[9],
        pressure_shift=record[10],
        local_quanta_lower=models.DiatomicOrLinearLocalQuanta(
            J=record[15],
            symmetry=record[16],
            F=record[17]
        ),
        local_quanta_upper=models.DiatomicOrLinearLocalQuanta(
            J=record[18],
            symmetry=record[19],
            F=record[20]
        ),
        rotational_statistical_weight_lower=record[13],
        rotational_statistical_weight_upper=record[14],
    )


def record_to_line_model_triatomic(record):
    return models.Line(
        line_id=record[1],
        structure_id=record[2],
        wavenumber_vacuum=record[3],
        line_strength_296K=record[4],
        einstein_A=record[5],
        broadening_hw_air=record[6],
        broadening_hw_self=record[7],
        energy_lower_state=record[8],
        broadening_temp_coefficient=record[9],
        pressure_shift=record[10],
        local_quanta_lower=models.TriatomicNonLinearLocalQuanta(
            J=record[15],
            K=record[16],
        ),
        local_quanta_upper=models.TriatomicNonLinearLocalQuanta(
            J=record[17],
            K=record[18]
        ),
        rotational_statistical_weight_lower=record[13],
        rotational_statistical_weight_upper=record[14],
    )


def record_to_line_model_pyramidal_tetratomic(record):
    """
    Takes a database record with all line fields selected and converts it to an line model object
    :param record: Array containing the values of all line database fields for a given row
    :return:
    """
    return models.Line(
        line_id=record[1],
        structure_id=record[2],
        wavenumber_vacuum=record[3],
        line_strength_296K=record[4],
        einstein_A=record[5],
        broadening_hw_air=record[6],
        broadening_hw_self=record[7],
        energy_lower_state=record[8],
        broadening_temp_coefficient=record[9],
        pressure_shift=record[10],
        local_quanta_lower=models.PyramidalTetratomicLocalQuanta(
            J=record[15],
            K=record[16],
            inversion_symmetry=record[17],
            symmetry_rotational=record[18],
            symmetry_total=record[19]
        ),
        local_quanta_upper=models.PyramidalTetratomicLocalQuanta(
            J=record[20],
            K=record[21],
            inversion_symmetry=record[22],
            symmetry_rotational=record[23],
            symmetry_total=record[24]
        ),
        rotational_statistical_weight_lower=record[13],
        rotational_statistical_weight_upper=record[14],
    )


class ResultSet:
    """
    Object to represent a set of results from the database. This object implements a Python interator to
    step through the records one-by-one to prevent high memory usage for large datasets. Alternatively, the full set of
    results can be loaded using the all() method.
    """
    def __init__(self, cursor: sqlite3.Cursor):
        self._cursor = cursor

    def __iter__(self):
        """
        Return the iterator object
        :return: 
        """
        return self


class MoleculeResultSet(ResultSet):
    """
    Object to represent a set of molecule results from the database. This object implements a Python interator to
    step through the records one-by-one to prevent high memory usage for large datasets. Alternatively, the full set of
    results can be loaded using the all() method.
    """
    def __next__(self) -> models.Molecule:
        """
        Fetch the next record from the database, and convert it to a Molecule object
        :return: 
        """
        record = self._cursor.fetchone()

        if record is not None:
            return record_to_molecule_model(record)
        else:
            raise StopIteration

    def all(self) -> List[models.Molecule]:
        """
        Fetch the full list of molecules from the resultset and convert them to Molecule objects
        :return: 
        """
        return [molecule for molecule in self]


class IsotopologueResultSet(ResultSet):
    """
    Object to represent a set of isotopologue results from the database. This object implements a Python interator to
    step through the records one-by-one to prevent high memory usage for large datasets. Alternatively, the full set of
    results can be loaded using the all() method.
    """
    def __next__(self) -> models.Isotopologue:
        """
        Fetch the next record from the database, and convert it to an isotopologue object
        :return: 
        """
        record = self._cursor.fetchone()

        if record is not None:
            return record_to_isotopologue_model(record)
        else:
            raise StopIteration

    def all(self) -> List[models.Isotopologue]:
        """
        Fetch the full list of isotopologues from the resultset and convert them to isotopologue objects
        :return: 
        """
        return [isotopologue for isotopologue in self]


class StructureResultSet(ResultSet):
    """
    Object to represent a set of structure results from the database. This object implements a Python interator to
    step through the records one-by-one to prevent high memory usage for large datasets. Alternatively, the full set of
    results can be loaded using the all() method.
    """
    def __next__(self) -> models.Structure:
        """
        Fetch the next record from the database, and convert it to a structure object
        :return: 
        """
        record = self._cursor.fetchone()

        if record is not None:
            if models.Structure.get_global_quanta_model_type_by_molecule_id(record[1]) == models.DiatomicGlobalQuanta:
                return record_to_structure_model_diatomic(record)
            elif models.Structure.get_global_quanta_model_type_by_molecule_id(record[1]) == models.TriatomicLinearFermiResonantGlobalQuanta:
                return record_to_structure_model_triatomic_linear_fermi_resonant(record)
            elif models.Structure.get_global_quanta_model_type_by_molecule_id(record[1]) == models.TriatomicGlobalQuanta:
                return record_to_structure_triatomic(record)
            elif models.Structure.get_global_quanta_model_type_by_molecule_id(record[1]) == models.PyramidalTetratomicGlobalQuanta:
                return record_to_structure_model_pyramidal_tetratomic(record)
            return record_to_structure_model(record)
        else:
            raise StopIteration

    def all(self) -> List[models.Structure]:
        """
        Fetch the full list of structures from the resultset and convert them to structure objects
        :return: 
        """
        return [structure for structure in self]


class LineResultSet(ResultSet):
    """
    Object to represent a set of line results from the database. This object implements a Python interator to
    step through the records one-by-one to prevent high memory usage for large datasets. Alternatively, the full set of
    results can be loaded using the all() method.
    """
    def __next__(self) -> models.Line:
        """
        Fetch the next record from the database, and convert it to a line object
        :return: 
        """
        record = self._cursor.fetchone()

        if record is not None:
            if models.Line.get_local_quanta_model_type_by_molecule_id(record[0]) == models.DiatomicOrLinearLocalQuanta:
                return record_to_line_model_diatomic_or_linear(record)
            elif models.Line.get_local_quanta_model_type_by_molecule_id(record[0]) == models.TriatomicNonLinearLocalQuanta:
                return record_to_line_model_triatomic(record)
            elif models.Line.get_local_quanta_model_type_by_molecule_id(record[0]) == models.PyramidalTetratomicLocalQuanta:
                return record_to_line_model_pyramidal_tetratomic(record)
            return record_to_line_model(record)
        else:
            raise StopIteration

    def all(self) -> List[models.Line]:
        """
        Fetch the full list of lines from the resultset and convert them to line objects
        :return: 
        """
        return [line for line in self]


class PartitionSumSet:
    """
    Object to represent the temperature dependent partition sum of a given isotopologue. This partition sum in stored
    as a series of datapoints (two seperate arrays) with a partition sum for a given set of temperatures.
    """
    def __init__(self, temperatures, partition_sums):
        self._temperatures = temperatures
        self._partition_sums = partition_sums

        self._sort_data()

    def _sort_data(self):
        """
        If interpolation is used to determine the partition sum in between datapoints, the data should be sorted on
        temperature incrementally. This method performs the sorting operation.
        :return: 
        """
        sort_indexes = np.argsort(self._temperatures)

        self._temperatures = [self._temperatures[index] for index in sort_indexes]
        self._partition_sums = [self._partition_sums[index] for index in sort_indexes]

    def get_partition_sum_for_temperature(self, temperature):
        """
        Return the partition sum for the given temperature passed as a parameter.
        :param temperature: The temperature for which the partition sum should be returned
        :return: The value of the partition sum at the given temperature
        """
        if temperature >= min(self._temperatures) and temperature <= max(self._temperatures):
            return np.interp(temperature, self._temperatures, self._partition_sums)

        raise ValueError("No partition sum can be calculated for the given temperature, since the temperature is ouside of the temperature range for which partition sums are known")
