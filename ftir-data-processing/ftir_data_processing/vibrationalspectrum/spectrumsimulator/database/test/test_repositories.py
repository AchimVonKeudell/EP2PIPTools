import unittest
from spectrumsimulator.database import hitran, repositories, models, queries


class MoleculeInsertRepositoryTestCase(unittest.TestCase):
    def setUp(self):
        self._database = hitran.HitranDatabase(":memory:")
        self._database.executescript(queries.molecule_table_create)
        self._repository = repositories.MoleculeRepository(self._database)

    def tearDown(self):
        self._database.close()

    def test_save_returns_inserted_id(self):
        molecule_id = 1
        molecule_name = "H2O"
        molecule = models.Molecule(
            molecule_id=molecule_id,
            molecule_name=molecule_name
        )

        inserted_id = self._repository.save(molecule)
        self.assertEqual(molecule_id, inserted_id)


class MoleculeQueryRepositoryTestCase(unittest.TestCase):
    def setUp(self):
        self._database = hitran.HitranDatabase(":memory:")
        self._database.executescript(queries.molecule_table_create)
        self._repository = repositories.MoleculeRepository(self._database)

        # insert two molecules
        self._molecule_h2o = models.Molecule(
            molecule_id=1,
            molecule_name="H2O"
        )
        self._molecule_co2 = models.Molecule(
            molecule_id=2,
            molecule_name="CO2"
        )

        self._repository.save(self._molecule_h2o)
        self._repository.save(self._molecule_co2)

    def tearDown(self):
        self._database.close()

    def test_get_all(self):
        molecules = self._repository.get_all().all()

        self.assertEqual(len(molecules), 2)

    def test_get_properties(self):
        h2o_molecule = self._repository.get_by_id(self._molecule_h2o.molecule_id)

        self.assertEqual(self._molecule_h2o.molecule_id, h2o_molecule.molecule_id)
        self.assertEqual(self._molecule_h2o.molecule_name, h2o_molecule.molecule_name)

    def test_get_non_existing(self):
        with self.assertRaises(IndexError):
            self._repository.get_by_id(6)


class IsotopologueInsertRepositoryTestCase(unittest.TestCase):
    def setUp(self):
        self._database = hitran.HitranDatabase(":memory:")
        self._database.executescript(queries.molecule_table_create)
        self._database.executescript(queries.isotopologue_table_create)
        self._repository = repositories.IsotopologueRepository(self._database)

    def tearDown(self):
        self._database.close()

    def test_save_returns_isotopologue_object(self):
        molecule_id = 1
        isotopologue_order_num = 1
        afgl_code = "161"
        abundance = 0.996
        molar_mass = 16.04

        isotopologue = models.Isotopologue(
            molecule_id=molecule_id,
            isotopologue_order_num=isotopologue_order_num,
            afgl_code=afgl_code,
            abundance=abundance,
            molar_mass=molar_mass
        )

        self._repository.save(isotopologue)
        inserted_isotopologue = self._repository.get_by_isotopologue_order_num(molecule_id, isotopologue_order_num)

        self.assertIsInstance(inserted_isotopologue, models.Isotopologue)


class IsotopologueQueryRepositoryTestCase(unittest.TestCase):
    def setUp(self):
        self._database = hitran.HitranDatabase(":memory:")
        self._database.executescript(queries.molecule_table_create)
        self._database.executescript(queries.isotopologue_table_create)
        self._repository = repositories.IsotopologueRepository(self._database)

        # insert three isotopologues
        self._isotopologue_h2o_1 = models.Isotopologue(
            molecule_id=1,
            isotopologue_order_num=1,
            afgl_code="161",
            abundance=0.99,
            molar_mass=1
        )
        self._isotopologue_h2o_2 = models.Isotopologue(
            molecule_id=1,
            isotopologue_order_num=2,
            afgl_code="171",
            abundance=0.01,
            molar_mass=2
        )
        self._isotopologue_co2_1 = models.Isotopologue(
            molecule_id=2,
            isotopologue_order_num=1,
            afgl_code="282",
            abundance=0.99,
            molar_mass=3
        )

        self._repository.save(self._isotopologue_h2o_1)
        self._repository.save(self._isotopologue_h2o_2)
        self._repository.save(self._isotopologue_co2_1)

    def tearDown(self):
        self._database.close()

    def test_get_properties(self):
        isotopologue_h2o_2 = self._repository.get_by_isotopologue_order_num(
            self._isotopologue_h2o_2.molecule_id,
            self._isotopologue_h2o_2.isotopologue_order_num
        )

        self.assertEqual(self._isotopologue_h2o_2.molecule_id, isotopologue_h2o_2.molecule_id)
        self.assertEqual(self._isotopologue_h2o_2.isotopologue_order_num, isotopologue_h2o_2.isotopologue_order_num)
        self.assertEqual(self._isotopologue_h2o_2.afgl_code, isotopologue_h2o_2.afgl_code)
        self.assertEqual(self._isotopologue_h2o_2.abundance, isotopologue_h2o_2.abundance)
        self.assertEqual(self._isotopologue_h2o_2.molar_mass, isotopologue_h2o_2.molar_mass)

    def test_get_isotopologues_by_molecule_id(self):
        h2o_isotopologues = self._repository.get_by_molecule_id(self._isotopologue_h2o_1.molecule_id).all()

        self.assertEqual(len(h2o_isotopologues), 2)

    def test_get_non_existing(self):
        with self.assertRaises(IndexError):
            self._repository.get_by_isotopologue_order_num(2, 2)


class StructureInsertRepositoryTestCase(unittest.TestCase):
    def setUp(self):
        self._database = hitran.HitranDatabase(":memory:")
        self._database.executescript(queries.molecule_table_create)
        self._database.executescript(queries.isotopologue_table_create)
        self._database.executescript(queries.structure_table_create)
        self._database.executescript(queries.structure_diatomic_table_create)
        self._database.executescript(queries.structure_triatomic_linear_fermi_resonant_table_create)
        self._repository = repositories.StructureRepository(self._database)

    def tearDown(self):
        self._database.close()

    def test_save_returns_structure_object(self):
        molecule_id = 1
        isotopologue_order_num = 1
        global_quanta_lower = models.GlobalQuanta("0 0 0")
        global_quanta_upper = models.GlobalQuanta("0 0 1")
        max_line_strength_296K = 1.2

        structure = models.Structure(
            molecule_id=molecule_id,
            isotopologue_order_num=isotopologue_order_num,
            global_quanta_lower=global_quanta_lower,
            global_quanta_upper=global_quanta_upper,
            max_line_strength_296K=max_line_strength_296K
        )

        inserted_structure_id = self._repository.save(structure)
        inserted_structure = self._repository.get_by_id(inserted_structure_id)

        self.assertIsInstance(inserted_structure, models.Structure)
        self.assertIsInstance(inserted_structure.global_quanta_lower, models.GlobalQuanta)
        self.assertIsInstance(inserted_structure.global_quanta_upper, models.GlobalQuanta)

    def test_save_diatomic_returns_structure_object(self):
        molecule_id = 5
        isotopologue_order_num = 1
        global_quanta_lower = models.DiatomicGlobalQuanta(vibrational_quantum=0)
        global_quanta_upper = models.DiatomicGlobalQuanta(vibrational_quantum=1)
        max_line_strength_296K = 1.2

        structure = models.Structure(
            molecule_id=molecule_id,
            isotopologue_order_num=isotopologue_order_num,
            global_quanta_lower=global_quanta_lower,
            global_quanta_upper=global_quanta_upper,
            max_line_strength_296K=max_line_strength_296K
        )

        inserted_structure_id = self._repository.save(structure)
        inserted_structure = self._repository.get_by_id(inserted_structure_id)

        self.assertIsInstance(inserted_structure, models.Structure)
        self.assertIsInstance(inserted_structure.global_quanta_lower, models.DiatomicGlobalQuanta)
        self.assertIsInstance(inserted_structure.global_quanta_upper, models.DiatomicGlobalQuanta)

    def test_save_triatomic_linear_fermi_resonant_returns_structure_object(self):
        molecule_id = 2
        isotopologue_order_num = 1
        global_quanta_lower = models.TriatomicLinearFermiResonantGlobalQuanta(0, 1, 1, 1, 1)
        global_quanta_upper = models.TriatomicLinearFermiResonantGlobalQuanta(0, 1, 0, 1, 1)
        max_line_strength_296K = 1.2

        structure = models.Structure(
            molecule_id=molecule_id,
            isotopologue_order_num=isotopologue_order_num,
            global_quanta_lower=global_quanta_lower,
            global_quanta_upper=global_quanta_upper,
            max_line_strength_296K=max_line_strength_296K
        )

        inserted_structure_id = self._repository.save(structure)
        inserted_structure = self._repository.get_by_id(inserted_structure_id)

        self.assertIsInstance(inserted_structure, models.Structure)
        self.assertIsInstance(inserted_structure.global_quanta_lower, models.TriatomicLinearFermiResonantGlobalQuanta)
        self.assertIsInstance(inserted_structure.global_quanta_upper, models.TriatomicLinearFermiResonantGlobalQuanta)


class StructureQueryRepositoryTestCase(unittest.TestCase):
    def setUp(self):
        self._database = hitran.HitranDatabase(":memory:")
        self._database.executescript(queries.molecule_table_create)
        self._database.executescript(queries.isotopologue_table_create)
        self._database.executescript(queries.structure_table_create)
        self._database.executescript(queries.structure_diatomic_table_create)
        self._database.executescript(queries.structure_triatomic_linear_fermi_resonant_table_create)
        self._repository = repositories.StructureRepository(self._database)

        # insert three isotopologues
        self._structure_1 = models.Structure(
            molecule_id=1,
            isotopologue_order_num=1,
            global_quanta_lower=models.GlobalQuanta("0 0 0"),
            global_quanta_upper=models.GlobalQuanta("0 0 1"),
            max_line_strength_296K=0.18
        )
        self._structure_2 = models.Structure(
            molecule_id=1,
            isotopologue_order_num=1,
            global_quanta_lower=models.GlobalQuanta("0 0 0"),
            global_quanta_upper=models.GlobalQuanta("0 0 1"),
            max_line_strength_296K=0.18
        )
        self._structure_3 = models.Structure(
            molecule_id=1,
            isotopologue_order_num=2,
            global_quanta_lower=models.GlobalQuanta("0 0 0"),
            global_quanta_upper=models.GlobalQuanta("0 0 1"),
            max_line_strength_296K=0.18
        )
        self._structure_4 = models.Structure(
            molecule_id=5,
            isotopologue_order_num=1,
            global_quanta_lower=models.DiatomicGlobalQuanta(vibrational_quantum=0),
            global_quanta_upper=models.DiatomicGlobalQuanta(vibrational_quantum=1)
        )
        self._structure_5 = models.Structure(
            molecule_id=2,
            isotopologue_order_num=2,
            global_quanta_lower=models.TriatomicLinearFermiResonantGlobalQuanta(
                vibrational_quantum_1=0,
                vibrational_quantum_2=0,
                vibrational_quantum_3=0,
                vibrational_angular_momentum=0,
                vibrational_ranking_index=1
            ),
            global_quanta_upper=models.TriatomicLinearFermiResonantGlobalQuanta(
                vibrational_quantum_1=0,
                vibrational_quantum_2=1,
                vibrational_quantum_3=0,
                vibrational_angular_momentum=1,
                vibrational_ranking_index=0
            )
        )

        self._repository.save(self._structure_1)
        self._repository.save(self._structure_2)
        self._repository.save(self._structure_3)
        self._repository.save(self._structure_4)
        self._repository.save(self._structure_5)

    def tearDown(self):
        self._database.close()

    def test_get_properties(self):
        structure_1 = self._repository.get_by_id(1)

        self.assertEqual(self._structure_1.structure_id, structure_1.structure_id)
        self.assertEqual(self._structure_1.molecule_id, structure_1.molecule_id)
        self.assertEqual(self._structure_1.isotopologue_order_num, structure_1.isotopologue_order_num)
        self.assertEqual(self._structure_1.global_quanta_lower.hitran_raw_quanta, structure_1.global_quanta_lower.hitran_raw_quanta)
        self.assertEqual(self._structure_1.global_quanta_upper.hitran_raw_quanta, structure_1.global_quanta_upper.hitran_raw_quanta)
        self.assertEqual(self._structure_1.max_line_strength_296K, structure_1.max_line_strength_296K)

        self.assertEqual(self._structure_1.global_quanta_lower.hitran_raw_quanta, structure_1.global_quanta_lower.hitran_raw_quanta)
        self.assertEqual(self._structure_1.global_quanta_upper.hitran_raw_quanta, structure_1.global_quanta_upper.hitran_raw_quanta)

    def test_get_properties_diatomic(self):
        structure_4 = self._repository.get_by_id(4)

        self.assertEqual(self._structure_4.global_quanta_lower.hitran_raw_quanta, structure_4.global_quanta_lower.hitran_raw_quanta)
        self.assertEqual(self._structure_4.global_quanta_lower.vibrational_quantum, structure_4.global_quanta_lower.vibrational_quantum)
        self.assertEqual(self._structure_4.global_quanta_upper.hitran_raw_quanta, structure_4.global_quanta_upper.hitran_raw_quanta)
        self.assertEqual(self._structure_4.global_quanta_upper.vibrational_quantum, structure_4.global_quanta_upper.vibrational_quantum)

    def test_get_properties_triatomic_linear_fermi_resonant(self):
        structure_5 = self._repository.get_by_id(5)

        self.assertEqual(self._structure_5.global_quanta_lower.hitran_raw_quanta, structure_5.global_quanta_lower.hitran_raw_quanta)
        self.assertEqual(self._structure_5.global_quanta_lower.vibrational_quantum_1, structure_5.global_quanta_lower.vibrational_quantum_1)
        self.assertEqual(self._structure_5.global_quanta_lower.vibrational_quantum_2, structure_5.global_quanta_lower.vibrational_quantum_2)
        self.assertEqual(self._structure_5.global_quanta_lower.vibrational_quantum_3, structure_5.global_quanta_lower.vibrational_quantum_3)
        self.assertEqual(self._structure_5.global_quanta_lower.vibrational_angular_momentum, structure_5.global_quanta_lower.vibrational_angular_momentum)
        self.assertEqual(self._structure_5.global_quanta_lower.vibrational_ranking_index, structure_5.global_quanta_lower.vibrational_ranking_index)
        self.assertEqual(self._structure_5.global_quanta_upper.hitran_raw_quanta, structure_5.global_quanta_upper.hitran_raw_quanta)
        self.assertEqual(self._structure_5.global_quanta_upper.vibrational_quantum_1, structure_5.global_quanta_upper.vibrational_quantum_1)
        self.assertEqual(self._structure_5.global_quanta_upper.vibrational_quantum_2, structure_5.global_quanta_upper.vibrational_quantum_2)
        self.assertEqual(self._structure_5.global_quanta_upper.vibrational_quantum_3, structure_5.global_quanta_upper.vibrational_quantum_3)
        self.assertEqual(self._structure_5.global_quanta_upper.vibrational_angular_momentum, structure_5.global_quanta_upper.vibrational_angular_momentum)
        self.assertEqual(self._structure_5.global_quanta_upper.vibrational_ranking_index, structure_5.global_quanta_upper.vibrational_ranking_index)

    def test_get_structures_by_isotopologue_id(self):
        structures = self._repository.get_by_isotopologue_order_num(self._structure_1.molecule_id, self._structure_1.isotopologue_order_num).all()

        self.assertEqual(len(structures), 2)

    def test_get_non_existing(self):
        with self.assertRaises(IndexError):
            self._repository.get_by_id(6)


class LineInsertRepositoryTestCase(unittest.TestCase):
    def setUp(self):
        self._database = hitran.HitranDatabase(":memory:")
        self._database.executescript(queries.molecule_table_create)
        self._database.executescript(queries.isotopologue_table_create)
        self._database.executescript(queries.structure_table_create)
        self._database.executescript(queries.structure_diatomic_table_create)
        self._database.executescript(queries.line_table_create)
        self._database.executescript(queries.line_diatomic_or_linear_table_create)
        self._repository_lines = repositories.LineRepository(self._database)
        self._repository_structures = repositories.StructureRepository(self._database)

    def tearDown(self):
        self._database.close()

    def test_save_returns_line_object(self):
        structure = models.Structure(
            molecule_id=1,
            isotopologue_order_num=1,
            global_quanta_lower=models.LocalQuanta(hitran_raw_quanta="0"),
            global_quanta_upper=models.LocalQuanta(hitran_raw_quanta="1")
        )

        self._repository_structures.save(structure)

        line = models.Line(
            structure_id=structure.structure_id,
            wavenumber_vacuum=2000.0,
            line_strength_296K=1e-20,
            einstein_A=2.3,
            broadening_hw_air=1.1,
            broadening_hw_self=1.4,
            energy_lower_state=100.0,
            broadening_temp_coefficient= 1.4,
            pressure_shift=-0.1,
            local_quanta_lower=models.LocalQuanta(hitran_raw_quanta="12"),
            local_quanta_upper=models.LocalQuanta(hitran_raw_quanta="13"),
            rotational_statistical_weight_lower=6.0,
            rotational_statistical_weight_upper=3.0
        )

        inserted_line_id = self._repository_lines.save(line, structure.molecule_id)
        inserted_line = self._repository_lines.get_by_id(inserted_line_id)

        self.assertIsInstance(inserted_line, models.Line)
        self.assertIsInstance(inserted_line.local_quanta_lower, models.LocalQuanta)
        self.assertIsInstance(inserted_line.local_quanta_upper, models.LocalQuanta)


class LineQueryRepositoryTestCase(unittest.TestCase):
    def _get_diatomic_linear_local_quanta(self, J, symmetry, F):
        return models.DiatomicOrLinearLocalQuanta(J, symmetry, F)

    def setUp(self):
        self._database = hitran.HitranDatabase(":memory:")
        self._database.executescript(queries.molecule_table_create)
        self._database.executescript(queries.isotopologue_table_create)
        self._database.executescript(queries.structure_table_create)
        self._database.executescript(queries.structure_diatomic_table_create)
        self._database.executescript(queries.line_table_create)
        self._database.executescript(queries.line_diatomic_or_linear_table_create)
        self._repository_lines = repositories.LineRepository(self._database)
        self._repository_structures = repositories.StructureRepository(self._database)

        # insert two structures to link the lines to
        self._structure_1 = models.Structure(
            molecule_id=1,
            isotopologue_order_num=1,
            global_quanta_lower=models.GlobalQuanta("0 0 0"),
            global_quanta_upper=models.GlobalQuanta("0 0 1"),
            max_line_strength_296K=0.18
        )
        self._structure_2 = models.Structure(
            molecule_id=5,
            isotopologue_order_num=1,
            global_quanta_lower=models.DiatomicGlobalQuanta(vibrational_quantum=0),
            global_quanta_upper=models.DiatomicGlobalQuanta(vibrational_quantum=1),
            max_line_strength_296K=0.18
        )

        self._repository_structures.save(self._structure_1)
        self._repository_structures.save(self._structure_2)

        # insert three lines
        self._line_1 = models.Line(
            structure_id=self._structure_1.structure_id,
            wavenumber_vacuum=2000.0,
            line_strength_296K=1e-20,
            einstein_A=2.3,
            broadening_hw_air=1.1,
            broadening_hw_self=1.4,
            energy_lower_state=100.0,
            broadening_temp_coefficient=1.5,
            pressure_shift=-0.1,
            local_quanta_lower=models.LocalQuanta(hitran_raw_quanta="12"),
            local_quanta_upper=models.LocalQuanta(hitran_raw_quanta="13"),
            rotational_statistical_weight_lower=6.0,
            rotational_statistical_weight_upper=3.0
        )
        self._line_2 = models.Line(
            structure_id=self._structure_1.structure_id,
            wavenumber_vacuum=2100.0,
            line_strength_296K=1e-21,
            einstein_A=2.4,
            broadening_hw_air=1.2,
            broadening_hw_self=1.5,
            energy_lower_state=110.0,
            broadening_temp_coefficient=1.6,
            pressure_shift=-0.2,
            local_quanta_lower=models.LocalQuanta(hitran_raw_quanta="14"),
            local_quanta_upper=models.LocalQuanta(hitran_raw_quanta="15"),
            rotational_statistical_weight_lower=7.0,
            rotational_statistical_weight_upper=4.0
        )
        self._line_3 = models.Line(
            structure_id=self._structure_1.structure_id,
            wavenumber_vacuum=2200.0,
            line_strength_296K=1e-22,
            einstein_A=2.5,
            broadening_hw_air=1.3,
            broadening_hw_self=1.6,
            energy_lower_state=120.0,
            broadening_temp_coefficient=1.7,
            pressure_shift=-0.3,
            local_quanta_lower=models.LocalQuanta(hitran_raw_quanta="17"),
            local_quanta_upper=models.LocalQuanta(hitran_raw_quanta="18"),
            rotational_statistical_weight_lower=8.0,
            rotational_statistical_weight_upper=5.0
        )
        self._line_4 = models.Line(
            structure_id=self._structure_2.structure_id,
            wavenumber_vacuum=2200.0,
            line_strength_296K=1e-22,
            einstein_A=2.5,
            broadening_hw_air=1.3,
            broadening_hw_self=1.6,
            energy_lower_state=120.0,
            broadening_temp_coefficient=1.7,
            pressure_shift=-0.3,
            local_quanta_lower=self._get_diatomic_linear_local_quanta(J=17, F="a", symmetry="+"),
            local_quanta_upper=self._get_diatomic_linear_local_quanta(J=18, F="b", symmetry="-"),
            rotational_statistical_weight_lower=8.0,
            rotational_statistical_weight_upper=5.0
        )

        self._repository_lines.save(self._line_1, self._structure_1.molecule_id)
        self._repository_lines.save(self._line_2, self._structure_1.molecule_id)
        self._repository_lines.save(self._line_3, self._structure_1.molecule_id)
        self._repository_lines.save(self._line_4, self._structure_2.molecule_id)

    def tearDown(self):
        self._database.close()

    def test_get_properties(self):
        line_1 = self._repository_lines.get_by_id(1)

        self.assertEqual(self._line_1.line_id, line_1.structure_id)
        self.assertEqual(self._line_1.structure_id, line_1.structure_id)
        self.assertEqual(self._line_1.wavenumber_vacuum, line_1.wavenumber_vacuum)
        self.assertEqual(self._line_1.line_strength_296K, line_1.line_strength_296K)
        self.assertEqual(self._line_1.einstein_A, line_1.einstein_A)
        self.assertEqual(self._line_1.broadening_hw_air, line_1.broadening_hw_air)
        self.assertEqual(self._line_1.broadening_hw_self, line_1.broadening_hw_self)
        self.assertEqual(self._line_1.energy_lower_state, line_1.energy_lower_state)
        self.assertEqual(self._line_1.broadening_temp_coefficient, line_1.broadening_temp_coefficient)
        self.assertEqual(self._line_1.pressure_shift, line_1.pressure_shift)
        self.assertEqual(self._line_1.local_quanta_lower.hitran_raw_quanta, line_1.local_quanta_lower.hitran_raw_quanta)
        self.assertEqual(self._line_1.local_quanta_upper.hitran_raw_quanta, line_1.local_quanta_upper.hitran_raw_quanta)
        self.assertEqual(self._line_1.rotational_statistical_weight_lower, line_1.rotational_statistical_weight_lower)
        self.assertEqual(self._line_1.rotational_statistical_weight_upper, line_1.rotational_statistical_weight_upper)

    def test_get_properties_diatomic_or_linear(self):
        line_4 = self._repository_lines.get_by_id(4)

        self.assertIsInstance(line_4.local_quanta_lower, models.DiatomicOrLinearLocalQuanta)
        self.assertIsInstance(line_4.local_quanta_upper, models.DiatomicOrLinearLocalQuanta)

        self.assertEqual(self._line_4.local_quanta_lower.hitran_raw_quanta, line_4.local_quanta_lower.hitran_raw_quanta)
        self.assertEqual(self._line_4.local_quanta_lower.J, line_4.local_quanta_lower.J)
        self.assertEqual(self._line_4.local_quanta_lower.symmetry, line_4.local_quanta_lower.symmetry)
        self.assertEqual(self._line_4.local_quanta_lower.F, line_4.local_quanta_lower.F)
        self.assertEqual(self._line_4.local_quanta_upper.hitran_raw_quanta, line_4.local_quanta_upper.hitran_raw_quanta)
        self.assertEqual(self._line_4.local_quanta_upper.J, line_4.local_quanta_upper.J)
        self.assertEqual(self._line_4.local_quanta_upper.symmetry, line_4.local_quanta_upper.symmetry)
        self.assertEqual(self._line_4.local_quanta_upper.F, line_4.local_quanta_upper.F)

    def test_get_structures_by_isotopologue_id(self):
        lines = self._repository_lines.get_by_structure_id(self._line_1.structure_id).all()

        self.assertEqual(len(lines), 3)

    def test_get_non_existing(self):
        with self.assertRaises(IndexError):
            self._repository_lines.get_by_id(5)


if __name__ == '__main__':
    unittest.main()
