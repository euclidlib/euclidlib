import unittest
from numpy import ndarray
from euclidlib import spectro

class TestPkReadingRoutine(unittest.TestCase):

    def test_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            pk = spectro.power_spectrum("random_filename.nowhere")

    def test_invalid_filename_type(self):
        for bad_input in (0, [], {}):
            with self.assertRaises(TypeError):
                pk = spectro.power_spectrum(bad_input)

    def test_not_fits(self):
        test_file_names = ["data/dummy_pk.fits", "data/dummy_ascii.txt"]
        for not_fits_input in (test_file_names[1], test_file_names):
            with self.assertRaises(ValueError):
                pk = spectro.power_spectrum(not_fits_input)

    def test_bad_fits(self):
        test_file_names_list_pkbk = [
            [
                ["data/dummy_{}.fits".format(stat), "data/dummy_noextname_{}.fits".format(stat)],
                ["data/dummy_{}.fits".format(stat), "data/dummy_nopk_{}.fits".format(stat)]
            ] for stat in ("pk", "bk")
        ]
        for test_file_names_list in test_file_names_list_pkbk:
            for test_file_names in test_file_names_list:
                for bad_fits_input in (test_file_names[1], test_file_names):
                    with self.assertRaises(ValueError):
                        pk = spectro.power_spectrum(bad_fits_input)

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
        self.assertTrue(all(isinstance(pk["POS", "POS", (n+1), (n+1)], ndarray) for n in range(bin_number)))


if __name__ == "__main__":
    unittest.main()
