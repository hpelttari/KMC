import unittest
import numpy as np
from initialization_functions import create_grid
from events import adsorption, desorption, diffusion


class TestEvents(unittest.TestCase):

    def setUp(self):
        self.grid = np.ones((3, 3))

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

        
if __name__ == "__main__":
    unittest.main()
