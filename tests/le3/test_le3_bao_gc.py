import euclidlib as el


def test_bao_gc_with_dummy_data(data_path):
    # Load the data
    data = el.le3.bao_gc.BAO_alphas(data_path / "EUC_LE3_BAO_ALPHAS_Z1.0.fits")

    assert data is not None, "Data should not be None"
    assert isinstance(data, dict), "Data should be a dictionary"
    assert ("SPE", "SPE", 0, 0) in data, "Data should contain ('SPE', 'SPE', 0, 0)"
    assert data[("SPE", "SPE", 0, 0)] is not None, "BAO alphas should not be None"
