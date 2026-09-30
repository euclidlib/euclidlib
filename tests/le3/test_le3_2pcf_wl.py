import euclidlib as el


def test_correlation_functions(data_path):
    xis = el.le3.twopcf_wl.correlation_functions(
        data_path / "two_point_correlation_function_Xi.fits"
    )

    expected_keys = [
        ("SHE", "SHE", 2, 2),
        ("SHE", "SHE", 2, 5),
        ("SHE", "SHE", 5, 5),
    ]

    assert list(xis.keys()) == expected_keys


def _make_cosebis(nbins=3, nmodes=5):
    import numpy as np
    from cosmolib.data import COSEBI

    rng = np.random.default_rng(42)
    mode = np.arange(1, nmodes + 1)
    cb = {}
    for i in range(1, nbins + 1):
        for j in range(i, nbins + 1):
            ee, bb, eb = rng.normal(size=(3, nmodes))
            arr = np.array([[ee, eb], [eb, bb]])
            cb[("SHE", "SHE", i, j)] = COSEBI(
                array=arr, axis=(2,), mode=mode, nmodes=nmodes, thmin=1.0, thmax=400.0
            )
    return cb


def test_cosebis_roundtrip(tmp_path):
    import fitsio
    import numpy as np

    cb = _make_cosebis()
    path = tmp_path / "cosebis.fits"
    el.le3.twopcf_wl.cosebis.write(path, cb)

    # columns on disk hold the diagonal / off-diagonal components
    raw = fitsio.read(path, ext="SHEARSHEAR2D_COSEBI_1_2")
    np.testing.assert_array_equal(raw["EE"], cb[("SHE", "SHE", 1, 2)].array[0, 0])
    np.testing.assert_array_equal(raw["BB"], cb[("SHE", "SHE", 1, 2)].array[1, 1])
    np.testing.assert_array_equal(raw["EB"], cb[("SHE", "SHE", 1, 2)].array[0, 1])

    back = el.le3.twopcf_wl.cosebis(path)
    assert list(back.keys()) == list(cb.keys())
    for key, c in back.items():
        assert c.array.shape == (2, 2, 5)
        assert c.axis == (2,)
        np.testing.assert_array_equal(c.array, cb[key].array)
        np.testing.assert_array_equal(c.mode, cb[key].mode)
        assert (c.thmin, c.thmax, c.nmodes) == (1.0, 400.0, 5)


def test_cosebis_write_rejects_wrong_shape(tmp_path):
    import numpy as np
    import pytest
    from cosmolib.data import COSEBI

    bad = {
        ("SHE", "SHE", 1, 1): COSEBI(
            array=np.zeros((3, 5)),
            mode=np.arange(1, 6),
            nmodes=5,
            thmin=1.0,
            thmax=400.0,
        )
    }
    with pytest.raises(ValueError, match=r"\(2, 2, n_modes\)"):
        el.le3.twopcf_wl.cosebis.write(tmp_path / "bad.fits", bad)
