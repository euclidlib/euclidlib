import euclidlib as el


def test_correlation_functions(data_path):
    xis = el.photo.correlation_functions(
        data_path / "two_point_correlation_function_Xi.fits"
    )

    expected_keys = ["THETA", "2-2", "2-5", "5-5", "THMIN", "THMAX"]

    assert list(xis.keys()) == expected_keys
