import unittest
from numpy import ndarray
from euclidlib import spectro
from euclidlib.spectro._datamodel import Result

class TestPkReadingRoutine(unittest.TestCase):

    def test_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            pk = spectro.power_spectrum("random_filename.nowhere")

    def test_invalid_filename_type(self):
        for bad_input in (0, [], {}):
            with self.assertRaises(TypeError):
                pk = spectro.power_spectrum(bad_input)

    def test_not_fits(self):
        with self.assertRaises(ValueError):
            pk = spectro.power_spectrum("data/dummy_ascii.txt")

    def test_bad_fits(self):
        bad_fits_input_lists_pkbk = [
            [
                "data/dummy_noextname_{}.fits".format(stat),
                "data/dummy_badextname_{}.fits".format(stat),
                "data/dummy_nopk_{}.fits".format(stat)
            ] for stat in ("pk", "bk")
        ]
        for bad_fits_input_list in bad_fits_input_lists_pkbk:
            for bad_fits_input in bad_fits_input_list:
                with self.assertRaises(ValueError):
                    pk = spectro.power_spectrum(bad_fits_input)

    def test_read_single_pk(self):
        pk = spectro.power_spectrum("data/dummy_pk.fits")
        self.assertTrue(isinstance(pk, Result))

if __name__ == "__main__":
    unittest.main()
