import fitsio  # type: ignore
import numpy as np  # type: ignore
import pytest  # type: ignore

import euclidlib as el


@pytest.fixture
def bias_file(tmp_path):
    """mock she.biasParams product with the layout of the SGS files"""
    nbin = 3
    rng = np.random.default_rng(42)
    data = np.zeros(
        nbin + 1,
        dtype=[
            ("TOM_BIN_ID", ">i8"),
            ("NUM_OBJECTS", ">i8"),
            ("PARAMS", ">f8", (6,)),
            ("PARAM_COVARS", ">f8", (36,)),
            ("EMPIRICAL_PARAMS", ">f8", (4,)),
            ("EMPIRICAL_PARAM_COVARS", ">f8", (16,)),
        ],
    )
    data["TOM_BIN_ID"] = np.arange(nbin + 1)
    data["NUM_OBJECTS"] = [-1, *rng.integers(1000, 2000, nbin)]
    data["PARAMS"] = np.nan
    data["PARAMS"][:, :4] = rng.normal(0, 0.01, (nbin + 1, 4))
    data["PARAM_COVARS"] = np.nan
    for i in range(4):
        data["PARAM_COVARS"][:, 7 * i] = rng.uniform(1e-6, 1e-3, nbin + 1)
    data["EMPIRICAL_PARAMS"] = rng.normal(0, 0.01, (nbin + 1, 4))
    data["EMPIRICAL_PARAM_COVARS"] = rng.uniform(1e-6, 1e-3, (nbin + 1, 16))
    header = {
        "FITS_DEF": "she.biasParams",
        "METHOD": "LensMC",
        "MODEL": "NOT_COMPUTED",
        "EMPIRIC": "constant_shear_binned_weighted_linear_fit",
    }
    path = tmp_path / "bias.fits"
    fitsio.write(path, data, extname="SHEBIASPARAMS", header=header)
    return path, data


@pytest.mark.parametrize("empirical", [True, False])
def test_multiplicative_shear_bias(bias_file, empirical):
    path, data = bias_file
    bias = el.shecal.multiplicative_shear_bias(path, empirical=empirical)

    if empirical:
        par, cov = data["EMPIRICAL_PARAMS"], data["EMPIRICAL_PARAM_COVARS"]
        cov = cov.reshape(-1, 4, 4)
    else:
        par, cov = data["PARAMS"], data["PARAM_COVARS"]
        cov = cov.reshape(-1, 6, 6)[:, :4, :4]

    assert list(bias.keys()) == list(data["TOM_BIN_ID"])
    for i, (key, b) in enumerate(bias.items()):
        assert set(b) == {"m1", "m2", "c1", "c2", "cov", "num_objects"}
        assert b["num_objects"] == data["NUM_OBJECTS"][i]
        for j, name in enumerate(["m1", "m2", "c1", "c2"]):
            assert b[name] == par[i, j]
        assert b["cov"].shape == (4, 4)
        np.testing.assert_array_equal(b["cov"], cov[i])

    # the model parameters only have variances
    if not empirical:
        np.testing.assert_array_equal(np.diag(bias[1]["cov"]), np.diag(cov[1]))
        assert np.isnan(bias[1]["cov"][0, 1])


def test_multiplicative_shear_bias_bad_format(tmp_path):
    path = tmp_path / "bad.fits"
    data = np.zeros(2, dtype=[("TOM_BIN_ID", ">i8"), ("NUM_OBJECTS", ">i8")])
    fitsio.write(path, data, extname="SHEBIASPARAMS")
    with pytest.raises(ValueError, match="EMPIRICAL_PARAMS"):
        el.shecal.multiplicative_shear_bias(path)
