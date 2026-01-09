import unittest
from spectrumsimulator.database import results


class PartitionSumSetTestCase(unittest.TestCase):
    def setUp(self):
        temperatures = [100, 200, 300]
        partition_sums = [1, 2, 4]

        self._partition_sum_set = results.PartitionSumSet(temperatures, partition_sums)

    def tearDown(self):
        self._partition_sum_set = None

    def test_linear_interpolation(self):
        self.assertEqual(self._partition_sum_set.get_partition_sum_for_temperature(150), 1.5)
        self.assertEqual(self._partition_sum_set.get_partition_sum_for_temperature(250), 3.0)

    def test_upper_bound_error(self):
        with self.assertRaises(ValueError):
            self._partition_sum_set.get_partition_sum_for_temperature(350)

    def test_lower_bound_error(self):
        with self.assertRaises(ValueError):
            self._partition_sum_set.get_partition_sum_for_temperature(-50)
