import unittest
import numpy as np
from initialization_functions import create_grid, create_rates_list
from events import adsorption, desorption, diffusion, update_rates_list


class TestEvents(unittest.TestCase):


    def setUp(self):
        self.grid = np.ones((3, 3))
        self.p_adsorption = 1
        self.rates_list = create_rates_list(self.grid, self.p_adsorption)


    def test_adsorption(self):
        self.grid = adsorption(self.grid, 0, 0)
        self.assertEqual(self.grid[0][0], 2)


    def test_desorption(self):
        self.grid = desorption(self.grid, 0, 0)
        self.assertEqual(self.grid[0][0], 0)


    def test_diffusion(self):
        self.grid = diffusion(self.grid, 0, 0, 1, 1)
        self.assertEqual(self.grid[0][0], 0)
        self.assertEqual(self.grid[1][1], 2)


    def test_desorption_rates_are_correct(self):
        self.grid = adsorption(self.grid,0,0)
        self.rates_list = update_rates_list(self.grid, 0, 0, 0, 0, self.rates_list, self.p_adsorption, 0.1, 0)
        self.grid = adsorption(self.grid,0,1)
        self.rates_list = update_rates_list(self.grid, 0, 1, 0, 1,  self.rates_list, self.p_adsorption, 0.1, 0)
        self.assertEqual(self.rates_list[0], 1)
        self.assertEqual(self.rates_list[1], 0.1)
        self.assertEqual(self.rates_list[10],1)
        self.assertEqual(self.rates_list[11], 0.1)


if __name__ == "__main__":
    unittest.main()
