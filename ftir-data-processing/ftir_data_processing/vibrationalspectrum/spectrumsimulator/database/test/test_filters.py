import unittest
from pypika import Query, Table
from spectrumsimulator.database import filters


class MoleculeFilterTestCase(unittest.TestCase):
    def setUp(self):
        self.table = Table('molecules')
        self.query = Query().from_(self.table).select('molecule_id')

    def test_filter_by_id(self):
        query = filters.MoleculeIdMoleculeQueryFilter(1).filter_query(self.query)
        self.assertEqual(str(query), 'SELECT "molecule_id" FROM "molecules" WHERE "molecule_id"=1')

    def test_filter_by_id_list(self):
        query = filters.MoleculeIdListMoleculeQueryFilter([1, 2]).filter_query(self.query)
        self.assertEqual(str(query), 'SELECT "molecule_id" FROM "molecules" WHERE "molecule_id"=1 OR "molecule_id"=2')

    def test_filter_by_name(self):
        query = filters.MoleculeNameMoleculeQueryFilter("H2O").filter_query(self.query)
        self.assertEqual(str(query), 'SELECT "molecule_id" FROM "molecules" WHERE "molecule_name"=\'H2O\'')

    def test_filter_aggregate(self):
        composite_filter = filters.QueryFilterComposite()
        composite_filter.add(filters.MoleculeIdMoleculeQueryFilter(1))
        composite_filter.add(filters.MoleculeIdListMoleculeQueryFilter([1, 2]))
        composite_filter.add(filters.MoleculeNameMoleculeQueryFilter("H2O"))

        query = composite_filter.filter_query(self.query)
        self.assertEqual(str(query),
                         'SELECT "molecule_id" '
                         'FROM "molecules" '
                         'WHERE "molecule_id"=1 AND ("molecule_id"=1 OR "molecule_id"=2) AND "molecule_name"=\'H2O\''
                         )


class IsotopologueFilterTestCase(unittest.TestCase):
    def setUp(self):
        self.table = Table('isotopologues')
        self.query = Query().from_(self.table).select('molecule_id', 'isotopologue_order_num')

    def test_filter_by_isotopologue_id(self):
        query = filters.IsotopologueIdIsotopologueQueryFilter(1, 2).filter_query(self.query)
        self.assertEqual(str(query), 'SELECT "molecule_id","isotopologue_order_num" FROM "isotopologues" WHERE "molecule_id"=1 AND "isotopologue_order_num"=2')

    def test_filter_by_molecule_id(self):
        query = filters.MoleculeIdIsotopologueQueryFilter(1).filter_query(self.query)
        self.assertEqual(str(query), 'SELECT "molecule_id","isotopologue_order_num" FROM "isotopologues" WHERE "molecule_id"=1')

    def test_filter_most_abundant(self):
        query = filters.MostAbundantIsotopologueQueryFilter(3).filter_query(self.query)
        self.assertEqual(str(query), 'SELECT "molecule_id","isotopologue_order_num" FROM "isotopologues" ORDER BY "isotopologue_order_num" ASC LIMIT 3')


class StructureFilterTestCase(unittest.TestCase):
    def setUp(self):
        self.table = Table('structures')
        self.query = Query().from_(self.table).select('structure_id')

    def test_filter_by_structure_id(self):
        query = filters.StructureIdStructureQueryFilter(1).filter_query(self.query)
        self.assertEqual(str(query), 'SELECT "structure_id" FROM "structures" WHERE "structure_id"=1')

    def test_filter_by_isotopologue_id(self):
        query = filters.IsotopologueIdStructureQueryFilter(1, 2).filter_query(self.query)
        self.assertEqual(str(query), 'SELECT "structure_id" FROM "structures" WHERE "molecule_id"=1 AND "isotopologue_order_num"=2')

    def test_filter_strongest(self):
        query = filters.StrongestStructureQueryFilter(3).filter_query(self.query)
        self.assertEqual(str(query), 'SELECT "structure_id" FROM "structures" ORDER BY "max_line_strength_296K" DESC LIMIT 3')

    def test_filter_wavenumber_range(self):
        query = filters.WavenumberRangeStructureQueryFilter(1800, 2400).filter_query(self.query)
        self.assertEqual(str(query), 'SELECT "structure_id" FROM "structures" WHERE "min_wavenumber_vacuum"<2400 AND "max_wavenumber_vacuum">1800')


class LineFilterTestCase(unittest.TestCase):
    def setUp(self):
        self.table = Table('lines')
        self.query = Query().from_(self.table).select('wavenumber_vacuum')

    def test_filter_by_line_id(self):
        query = filters.LineIdLineQueryFilter(1).filter_query(self.query)
        self.assertEqual(str(query), 'SELECT "wavenumber_vacuum" FROM "lines" WHERE "line_id"=1')

    def test_filter_by_structure_id(self):
        query = filters.StructureIdLineQueryFilter(1).filter_query(self.query)
        self.assertEqual(str(query), 'SELECT "wavenumber_vacuum" FROM "lines" WHERE "structure_id"=1')

    def test_filter_wavenumber_range(self):
        query = filters.WavenumberRangeLineQueryFilter(1800, 2400).filter_query(self.query)
        self.assertEqual(str(query), 'SELECT "wavenumber_vacuum" FROM "lines" WHERE "wavenumber_vacuum" BETWEEN 1800 AND 2400')


if __name__ == '__main__':
    unittest.main()
