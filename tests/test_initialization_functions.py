import unittest
from initialization_functions import create_grid


class TestInitializationMethods(unittest.TestCase):

    def test_initialize_grid_size(self):
        self.assertEqual(create_grid(2).shape, (2, 2))


if __name__ == "__main__":
    unittest.main()
