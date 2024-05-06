import euclidlib as el


def test_angular_power_spectra(data_path):
    cls = el.photo.angular_power_spectra(data_path / "angular_power_spectra.fits")

    expected_keys = [
        ("P", "P", 0, 0),
        ("P", "G_E", 0, 0),
        ("P", "G_B", 0, 0),
        ("G_E", "G_E", 0, 0),
        ("G_B", "G_B", 0, 0),
        ("G_E", "G_B", 0, 0),
        ("P", "P", 0, 1),
        ("P", "G_E", 0, 1),
        ("P", "G_B", 0, 1),
        ("G_E", "G_E", 0, 1),
        ("G_B", "G_B", 0, 1),
        ("G_E", "G_B", 0, 1),
        ("P", "G_E", 1, 0),
        ("P", "G_B", 1, 0),
        ("G_E", "G_B", 1, 0),
        ("P", "P", 1, 1),
        ("P", "G_E", 1, 1),
        ("P", "G_B", 1, 1),
        ("G_E", "G_E", 1, 1),
        ("G_B", "G_B", 1, 1),
        ("G_E", "G_B", 1, 1),
    ]

    assert list(cls.keys()) == expected_keys


def test_mixing_matrices(data_path):
    mms = el.photo.mixing_matrices(data_path / "mixing_matrices.fits")

    expected_keys = [
        ("P", "P", 0, 1),
        ("P", "G_E", 1, 2),
        ("G_E", "G_E", 2, 3),
    ]

    assert list(mms.keys()) == expected_keys

def test_xi_tpcf(data_path):
    xis = el.photo.xi_tpcf(data_path / "two_point_correlation_function_Xi.fits")

    expected_keys = ["+",
                     "-"]

    assert list(xis.keys()) == expected_keys