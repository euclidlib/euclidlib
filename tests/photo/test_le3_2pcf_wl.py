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
