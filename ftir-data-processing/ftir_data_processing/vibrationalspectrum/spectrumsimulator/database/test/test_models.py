import unittest
from spectrumsimulator.database import models


class HitranRawQuantaTestCase(unittest.TestCase):
    def test_diatomic_global_quanta(self):
        quanta = models.DiatomicGlobalQuanta(vibrational_quantum=1)

        self.assertEqual(len(quanta.hitran_raw_quanta), 15)
        self.assertEqual(quanta.hitran_raw_quanta, "              1")

    def test_triatomic_global_quanta(self):
        quanta = models.TriatomicGlobalQuanta(
            vibrational_quantum_1=1,
            vibrational_quantum_2=2,
            vibrational_quantum_3=3
        )

        self.assertEqual(len(quanta.hitran_raw_quanta), 15)
        self.assertEqual(quanta.hitran_raw_quanta, "          1 2 3")

    def test_triatomic_linear_global_quanta(self):
        quanta = models.TriatomicLinearGlobalQuanta(
            vibrational_quantum_1=1,
            vibrational_quantum_2=2,
            vibrational_quantum_3=3,
            vibrational_angular_momentum=4
        )

        self.assertEqual(len(quanta.hitran_raw_quanta), 15)
        self.assertEqual(quanta.hitran_raw_quanta, "        1 2 4 3")

    def test_triatomic_linear_fermi_resonant_global_quanta(self):
        quanta = models.TriatomicLinearFermiResonantGlobalQuanta(
            vibrational_quantum_1=1,
            vibrational_quantum_2=2,
            vibrational_quantum_3=3,
            vibrational_angular_momentum=4,
            vibrational_ranking_index=5
        )

        self.assertEqual(len(quanta.hitran_raw_quanta), 15)
        self.assertEqual(quanta.hitran_raw_quanta, "       1 2 4 35")

    def test_diatomic_linear_local_quanta(self):
        quanta = models.DiatomicOrLinearLocalQuanta(
            J=15,
            symmetry="e",
            F="14.5"
        )

        self.assertEqual(len(quanta.hitran_raw_quanta), 15)
        self.assertEqual(quanta.hitran_raw_quanta, "       15e 14.5")


if __name__ == '__main__':
    unittest.main()
