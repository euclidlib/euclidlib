import numpy as np  # type: ignore
import pytest  # type: ignore

import euclidlib as el


@pytest.mark.parametrize("read_hist", [True, False])
@pytest.mark.parametrize("write_hist", [True, False])
def test_redshift_distributions(tmp_path, write_hist, read_hist):
    def fn(z, hist):
        """produce a mock n(z), optionally binned into a histogram"""
        nz = np.exp(-((z - 1.0) ** 2) / (0.5) ** 2 / 2)
        if hist:
            nz = (nz[:-1] + nz[1:]) / 2 * np.diff(z)
        return nz

    # write test data from higher redshift resolution
    z = np.linspace(0.0, 6.0, 10001)
    nz = fn(z, write_hist)

    path = tmp_path / "nz.fits"

    el.phz.redshift_distributions.write(path, z, nz, hist=write_hist)

    # read test data

    z, nz = el.phz.redshift_distributions(path, hist=read_hist)

    np.testing.assert_array_equal(z, np.linspace(0.0, 6.0, 3001))

    assert len(nz) == 1
    assert list(nz.keys()) == [1]

    nz_ = fn(z, read_hist)
    np.testing.assert_allclose(nz[1], nz_, rtol=1e-2, atol=1e-4)


def test_redshift_distributions_ident(tmp_path):
    """
    Test that writing a histogram in the correct format returns the data
    unchanged.
    """
    z_ = np.linspace(0.0, 6.0, 3001)
    zmid = (z_[:-1] + z_[1:]) / 2
    nz_ = np.exp(-((zmid - 1.0) ** 2) / (0.5) ** 2 / 2)

    path = tmp_path / "nz.fits"

    el.phz.redshift_distributions.write(path, z_, nz_, hist=True)
    z, nz = el.phz.redshift_distributions(path, hist=True)

    np.testing.assert_array_equal(z, z_)
    np.testing.assert_array_equal(nz[1], nz_.astype(nz[1].dtype))
