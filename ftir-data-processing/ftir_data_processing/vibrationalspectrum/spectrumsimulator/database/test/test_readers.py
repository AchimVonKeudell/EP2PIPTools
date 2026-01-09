import os
import unittest
from spectrumsimulator.database import readers


class HITRAN160ReaderTestCase(unittest.TestCase):
    def setUp(self):
        self._test_file = os.path.join(os.path.dirname(__file__), 'data/hitran_160.par')
        self._broken_test_file = os.path.join(os.path.dirname(__file__), 'data/hitran_160_broken.par')
        self._nonexisting_file = os.path.join(os.path.dirname(__file__), 'data/nonexisting.par')

    def test_init(self):
        reader = readers.HITRANReader(self._test_file)

    def test_nonexisting_file(self):
        with self.assertRaises(FileNotFoundError):
            reader = readers.HITRANReader(self._nonexisting_file)

    def test_read_all_lines(self):
        reader = readers.HITRANReader(self._test_file)
        line_count = 0

        for line_data in reader:
            line_count += 1

        self.assertEquals(line_count, 50)

    def test_incomplete_line_exception(self):
        reader = readers.HITRANReader(self._broken_test_file)

        with self.assertRaises(ValueError):
            for line_data in reader:
                pass


class HITRAN160LineDataTestCase(unittest.TestCase):
    def setUp(self):
        self._mapper = readers.HITRAN160ParameterMapper()

        with open(os.path.join(os.path.dirname(__file__), 'data/hitran_160.par'), 'r') as valid_file:
            self._valid_data_line = valid_file.readline().rstrip()

        with open(os.path.join(os.path.dirname(__file__), 'data/hitran_160_broken.par'), 'r') as broken_file:
            broken_file.readline()
            self._invalid_data_line = broken_file.readline().rstrip()

    def test_from_invalid_line(self):
        with self.assertRaises(ValueError):
            readers.HITRANLineData.from_line_string(self._invalid_data_line, self._mapper)

    def test_from_line_return_value(self):
        line_data = readers.HITRANLineData.from_line_string(self._valid_data_line)

        self.assertTrue(isinstance(line_data, readers.HITRANLineData))

    def test_data_correctness(self):
        line_data = readers.HITRANLineData.from_line_string(self._valid_data_line)

        self.assertEqual(line_data.molecule_id, 5)
        self.assertEqual(line_data.isotopologue_order_num, 5)
        self.assertEqual(line_data.wavenumber_vacuum, 3.462498)
        self.assertEqual(line_data.line_strength_296K, 1.599e-33)
        self.assertEqual(line_data.einstein_A, 3.155e-8)
        self.assertEqual(line_data.broadening_hw_air, 0.0797)
        self.assertEqual(line_data.broadening_hw_self, 0.086)
        self.assertEqual(line_data.energy_lower_state, 2043.6929)
        self.assertEqual(line_data.broadening_temp_coefficient, 0.76)
        self.assertEqual(line_data.pressure_shift, 0.000000)
        self.assertEqual(line_data.global_quanta_upper, "              1")
        self.assertEqual(line_data.global_quanta_lower, "              1")
        self.assertEqual(line_data.local_quanta_upper, "               ")
        self.assertEqual(line_data.local_quanta_lower, "     R  0      ")
        self.assertEqual(line_data.uncertainty_indices, [5, 8, 7, 6, 6, 0])
        self.assertEqual(line_data.reference_indices, [1, 1, 2, 2, 1, 0])
        self.assertEqual(line_data.flag, " ")
        self.assertEqual(line_data.rotational_statistical_weight_upper, 6.0)
        self.assertEqual(line_data.rotational_statistical_weight_lower, 2.0)
