import unittest
import numpy as np
from initialization_functions import create_grid, create_rates_list


class TestInitializationMethods(unittest.TestCase):

    def setUp(self):
        self.grid = create_grid(2)
        self.rates_list =create_rates_list(self.grid, 1)

    def test_initialize_grid_size(self):
        self.assertEqual(self.grid.shape, (2, 2))


    def test_initial_rates_are_correcct(self):
        initial_rates = np.zeros(2*2*10)
        initial_rates[0] = 1
        initial_rates[10] = 1
        initial_rates[20] = 1
        initial_rates[30] = 1
        self.assertListEqual(self.rates_list.tolist(), initial_rates.tolist())

if __name__ == "__main__":
    unittest.main()
