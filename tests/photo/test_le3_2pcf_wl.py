import euclidlib as el


def test_xi_tpcf(data_path):
    xis = el.photo.xi_tpcf(data_path / "two_point_correlation_function_Xi.fits")

    expected_keys = ["+", "-"]

    assert list(xis.keys()) == expected_keys
