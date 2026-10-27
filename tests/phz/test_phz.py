import numpy as np  # type: ignore
import fitsio

import euclidlib as el


def test_redshift_distributions_q1(data_path):
    z, nz = el.phz.redshift_distributions(data_path / "nz_q1.fits", version="q1")

    np.testing.assert_array_equal(z, np.linspace(0.0, 6.0, 3001)[:-1])

    assert list(nz.keys()) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]

    mean = np.array(
        [0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2.0]
    )
    var = (1 + mean) ** 2 * 0.01

    got = np.stack(list(nz.values()))
    expected = np.exp(-0.5 * (z - mean[:, None]) ** 2 / var[:, None])

    np.testing.assert_allclose(got, expected, atol=1e-10, rtol=1e-5)


def test_redshift_distributions_dr1(data_path):
    z, nz = el.phz.redshift_distributions(data_path / "nz_dr1.fits")

    np.testing.assert_array_equal(z, 0.01 * np.arange(600))

    assert list(nz.keys()) == [1, 2, 3, 4, 5, 6]

    mean = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    var = (1 + mean) ** 2 * 0.01

    got = np.stack(list(nz.values()))
    expected = np.exp(-0.5 * (z - mean[:, None]) ** 2 / var[:, None])

    np.testing.assert_allclose(got, expected, atol=1e-10, rtol=1e-5)


def test_redshift_distributions_write_dr1(tmp_path):
    z_step = 0.001
    z = z_step * np.arange(5000)

    bins = np.arange(1, 6)
    mean = np.array([0.5, 1.0, 1.5, 2.0, 2.5])
    var = (1 + mean) ** 2 * 0.001

    nz = np.exp(-0.5 * (z - mean[:, None]) ** 2 / var[:, None])

    nz_dict = {bin_id: dist for bin_id, dist in zip(bins, nz)}

    el.phz.redshift_distributions.write(tmp_path / "nz.fits", z, nz_dict)

    with fitsio.FITS(tmp_path / "nz.fits", "r") as fits:
        assert len(fits) == 2
        hdu = fits[-1]

        header = hdu.read_header()
        assert header["EXTNAME"] == "BIN_INFO"
        assert header["WEIGHT_METHOD"] == "NO_WEIGHT"
        assert header["BIN_TYPE"] == "TOM_BIN"
        assert header["NBIN"] == nz.shape[0]
        assert header["Z_STEP"] == z_step
        assert header["Z_STEP_NUMBER"] == z.size

        data = hdu.read()
        assert data.size == nz.shape[0]
        np.testing.assert_allclose(data["MEAN_REDSHIFT"], mean)
        np.testing.assert_allclose(data["VAR_REDSHIFT"], var)
        np.testing.assert_allclose(data["N_Z"], nz)
