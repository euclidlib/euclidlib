import numpy as np  # type: ignore
import pytest  # type: ignore
import euclidlib as el


@pytest.mark.parametrize("read_hist", [True, False])
@pytest.mark.parametrize("write_hist", [True, False])
def test_redshift_distributions(tmp_path, write_hist, read_hist):
    def fn(z, hist):
        nz = np.exp(-((z - 1.0) ** 2) / (0.5) ** 2 / 2)
        if hist:
            zmid = (z[:-1] + z[1:]) / 2
            nz = np.exp(-((zmid - 1.0) ** 2) / (0.5) ** 2 / 2) * np.diff(z)
            return {1: nz.astype(np.float32)}
        return {1: nz.astype(np.float32)}

    z = np.linspace(0.0, 6.0, 10001)
    nz = fn(z, write_hist)
    path = tmp_path / "nz.fits"

    el.photo.redshift_distributions.write(path, z, nz, hist=write_hist)
    z_read, nz_read = el.photo.redshift_distributions(path, hist=read_hist)

    np.testing.assert_array_equal(z_read, np.linspace(0.0, 6.0, 3001))
    assert list(nz_read.keys()) == [1]

    z_ref = np.linspace(0.0, 6.0, 3001)
    nz_ref = fn(z_ref, read_hist)
    np.testing.assert_allclose(nz_read[1], nz_ref[1], rtol=1e-1, atol=5e-2)


def test_redshift_distributions_ident(tmp_path):
    z_ = np.linspace(0.0, 6.0, 3001)
    zmid = (z_[:-1] + z_[1:]) / 2
    nz_ = np.exp(-((zmid - 1.0) ** 2) / (0.5) ** 2 / 2) * np.diff(z_)
    dict_nz = {1: nz_.astype(np.float32)}

    path = tmp_path / "nz.fits"
    el.photo.redshift_distributions.write(path, z_, dict_nz, hist=True)
    z_read, nz_read = el.photo.redshift_distributions(path, hist=True)

    np.testing.assert_array_equal(z_read, z_)
    np.testing.assert_allclose(nz_read[1], dict_nz[1], rtol=1e-1, atol=3e-4)
