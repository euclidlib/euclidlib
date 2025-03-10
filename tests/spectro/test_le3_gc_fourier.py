import unittest
from numpy import ndarray
from euclidlib import spectro

class TestPkReadingRoutine(unittest.TestCase):

    def test_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            pk = spectro.power_spectrum("random_filename.nothing")

    def test_invalid_filename_type(self):
        with self.assertRaises(TypeError):
            pk = spectro.power_spectrum(1)
            pk = spectro.power_spectrum([])
            pk = spectro.power_spectrum({})

    def test_not_fits(self):
        with self.assertRaises(ValueError):
            pk = spectro.power_spectrum("data/dummy_ascii.txt")
            pk = spectro.power_spectrum(
                ["data/dummy_pk.fits", "data/dummy_ascii.txt"]
            )

    def test_bad_fits(self):
        with self.assertRaises(ValueError):
            pk = spectro.power_spectrum("data/dummy_bad_pk.fits")
            pk = spectro.power_spectrum(
                ["data/dummy_pk.fits", "data/dummy_bad_pk.fits"]
            )

    def test_read_single_pk(self):
        pk = spectro.power_spectrum("data/dummy_pk.fits")
        self.assertTrue(isinstance(pk, dict))
        self.assertEqual(list(pk.keys()), [("POS", "POS", 1, 1),])
        self.assertTrue(isinstance(pk["POS", "POS", 1, 1], ndarray))

    def test_read_multiple_pk(self):
        bin_number = 4
        pk = spectro.power_spectrum(["data/dummy_pk.fits"]*bin_number)
        self.assertTrue(isinstance(pk, dict))
        self.assertEqual(list(pk.keys()), [("POS", "POS", (n+1), (n+1)) for n in range(bin_number)])


if __name__ == "__main__":
    unittest.main()
