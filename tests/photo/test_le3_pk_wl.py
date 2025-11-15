import euclidlib as el


def test_angular_power_spectra(data_path):
    cls = el.le3.pk_wl.angular_power_spectra(data_path / "angular_power_spectra.fits")

    expected_keys = {
        ("POS", "POS", 2, 2),
        ("POS", "POS", 2, 5),
        ("POS", "POS", 5, 5),
        ("SHE", "SHE", 2, 2),
        ("SHE", "SHE", 2, 5),
        ("SHE", "SHE", 5, 5),
        ("POS", "SHE", 2, 2),
        ("POS", "SHE", 2, 5),
        ("POS", "SHE", 5, 2),
        ("POS", "SHE", 5, 5),
    }

    assert set(cls.keys()) == expected_keys


def test_mixing_matrices(data_path):
    mms = el.le3.pk_wl.mixing_matrices(data_path / "mixing_matrices.fits")

    expected_keys = {
        ("POS", "POS", 2, 2),
        ("POS", "POS", 2, 5),
        ("POS", "POS", 5, 5),
        ("SHE", "SHE", 2, 2),
        ("SHE", "SHE", 2, 5),
        ("SHE", "SHE", 5, 5),
        ("POS", "SHE", 2, 2),
        ("POS", "SHE", 2, 5),
        ("POS", "SHE", 5, 2),
        ("POS", "SHE", 5, 5),
    }

    assert set(mms.keys()) == expected_keys
