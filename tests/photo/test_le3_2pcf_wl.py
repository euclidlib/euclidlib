import euclidlib as el


def test_two_point_correlation(data_path):
    xis = el.photo.two_point_correlation(
        data_path / "two_point_correlation_function_Xi.fits")

    expected_keys = ["THETA", "2-2", "2-5", "5-5", "THMIN", "THMAX"]

    assert list(xis.keys()) == expected_keys
