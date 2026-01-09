from pypika import Query, Table, Order


class QueryFilterInterface:
    """
    Interface defining the required methods for a filter of a pypika Query object. Implementation of this interface can
    define any mutations on the Query object to filter by column values, sort the output, etc.
    """
    def filter_query(self, query: Query):
        """
        Perform a filter operation on a Query object, typically by calling where() on the Query object. This method
        should return a Query object with the filter applied.
        :param query: A pypika Query object
        :return: A pypika Query object with the filter applied
        """
        pass


class QueryFilterComposite(QueryFilterInterface):
    """
    This object represents a query filter composite, containing multiple query filters, which can be added using the
    add method. All added filters are applied to the given Query object when calling filter_query
    """
    def __init__(self):
        self._query_filters = []

    def add(self, query_filter: QueryFilterInterface):
        """
        Add a QueryFilter object to the composite filter
        :param query_filter: QueryFilter object
        :return: 
        """
        self._query_filters.append(query_filter)

    def filter_query(self, query: Query):
        """
        Apply the composite query to a given pypika Query object
        :param query: A pypika Query object 
        :return: Query object
        """
        for query_filter in self._query_filters:
            query = query_filter.filter_query(query)

        return query


class MoleculeIdMoleculeQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the molecules table by a molecule's id
    """
    def __init__(self, molecule_id):
        self._molecule_id = molecule_id
        self._molecule_table = Table('molecules')

    def filter_query(self, query: Query):
        return query.where(self._molecule_table.molecule_id == self._molecule_id)


class MoleculeIdListMoleculeQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the molecules table by a list of molecule ids. Records are included in
    the results if the molecule id matches any of the ids in the list provided in the constructor.
    """
    def __init__(self, molecule_ids):
        self._molecule_ids = molecule_ids
        self._molecule_table = Table('molecules')

    def filter_query(self, query: Query):
        condition = self._molecule_table.molecule_id == self._molecule_ids[0]

        for molecule_id in self._molecule_ids[1:]:
            condition |= self._molecule_table.molecule_id == molecule_id

        return query.where(condition)


class MoleculeNameMoleculeQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the molecules table by a molecule's name
    """
    def __init__(self, molecule_name):
        self._molecule_name = molecule_name
        self._molecule_table = Table('molecules')

    def filter_query(self, query: Query):
        return query.where(self._molecule_table.molecule_name == self._molecule_name)


class IsotopologueIdIsotopologueQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the isotopologues table by an isotopologue's id
    """
    def __init__(self, molecule_id, isotopologue_order_num):
        self._molecule_id = molecule_id
        self._isotopologue_order_num = isotopologue_order_num
        self._isotopologue_table = Table('isotopologues')

    def filter_query(self, query: Query):
        return query.where((self._isotopologue_table.molecule_id == self._molecule_id) & (self._isotopologue_table.isotopologue_order_num == self._isotopologue_order_num))


class MoleculeIdIsotopologueQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the isotopologues table by a molecule's id
    """
    def __init__(self, molecule_id):
        self._molecule_id = molecule_id
        self._isotopologue_table = Table('isotopologues')

    def filter_query(self, query: Query):
        return query.where(self._isotopologue_table.molecule_id == self._molecule_id)


class MostAbundantIsotopologueQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the isotopologues table by a molecule's id
    """
    def __init__(self, isotopologue_count):
        self._isotopologue_count = isotopologue_count
        self._isotopologue_table = Table('isotopologues')

    def filter_query(self, query: Query):
        return query.orderby(self._isotopologue_table.isotopologue_order_num, order=Order.asc).limit(self._isotopologue_count)


class StructureIdStructureQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the structures table by a structure's id
    """
    def __init__(self, structure_id):
        self._structure_id = structure_id
        self._structure_table = Table('structures')

    def filter_query(self, query: Query):
        return query.where(self._structure_table.structure_id == self._structure_id)


class IsotopologueIdStructureQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the structures table by an isotopologue's id
    """
    def __init__(self, molecule_id, isotopologue_order_num):
        self._molecule_id = molecule_id
        self._isotopologue_order_num = isotopologue_order_num
        self._structure_table = Table('structures')

    def filter_query(self, query: Query):
        return query.where((self._structure_table.molecule_id == self._molecule_id) & (self._structure_table.isotopologue_order_num == self._isotopologue_order_num))


class StrongestStructureQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the structures table by the maximum line strength at 296K
    """
    def __init__(self, structure_count):
        self._structure_count = structure_count
        self._structure_table = Table('structures')

    def filter_query(self, query: Query):
        return query.orderby(self._structure_table.max_line_strength_296K, order=Order.desc).limit(self._structure_count)


class WavenumberRangeStructureQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the lines table by an structures's id
    """
    def __init__(self, wavenumber_min, wavenumber_max):
        self._wavenumber_min = wavenumber_min
        self._wavenumber_max = wavenumber_max
        self._structure_table = Table('structures')

    def filter_query(self, query: Query):
        return query.where((self._structure_table.min_wavenumber_vacuum < self._wavenumber_max) & (self._structure_table.max_wavenumber_vacuum > self._wavenumber_min))


class DiatomicVibrationalQuantumLowerStructureQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the diatomic structures table by the vibrational quantum
    """
    def __init__(self, vibrational_quantum_lower):
        self._vibrational_quantum_lower = vibrational_quantum_lower
        self._structure_table = Table('structures_diatomic')

    def filter_query(self, query: Query):
        return query.where(self._structure_table.vibrational_quantum_lower == self._vibrational_quantum_lower)


class TriatomicLinearFermiResonantVibrationalQuantum3LowerStructureQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the triatomic Fermi resonant structures table by the third vibrational
    quantum (asymmetric stretch)
    """
    def __init__(self, vibrational_quantum_3_lower):
        self._vibrational_quantum_3_lower = vibrational_quantum_3_lower
        self._structure_table = Table('structures_triatomic_linear_fermi_resonant')

    def filter_query(self, query: Query):
        return query.where(self._structure_table.vibrational_quantum_3_lower == self._vibrational_quantum_3_lower)


class PyramidalTetratomicStructureQueryFilter(QueryFilterInterface):
    def __init__(self, vibrational_quantum_1_lower: int):
        self._vibrational_quantum_1_lower = vibrational_quantum_1_lower
        self._structure_table = Table("structures_pyramidal_tetratomic")

    def filter_query(self, query: Query):
        return query.where(self._structure_table._vibrational_quantum_1_lower == self._vibrational_quantum_1_lower)


class LineIdLineQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the lines table by an line's id
    """
    def __init__(self, line_id):
        self._line_id = line_id
        self._line_table = Table('lines')

    def filter_query(self, query: Query):
        return query.where(self._line_table.line_id == self._line_id)


class StructureIdLineQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the lines table by an structures's id
    """
    def __init__(self, structure_id):
        self._structure_id = structure_id
        self._line_table = Table('lines')

    def filter_query(self, query: Query):
        return query.where(self._line_table.structure_id == self._structure_id)


class WavenumberRangeLineQueryFilter(QueryFilterInterface):
    """
    This object filters a pypika Query object on the lines table by an structures's id
    """
    def __init__(self, wavenumber_min, wavenumber_max):
        self._wavenumber_min = wavenumber_min
        self._wavenumber_max = wavenumber_max
        self._line_table = Table('lines')

    def filter_query(self, query: Query):
        return query.where(self._line_table.wavenumber_vacuum[self._wavenumber_min:self._wavenumber_max])
