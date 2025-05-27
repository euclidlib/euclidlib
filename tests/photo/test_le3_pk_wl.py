import euclidlib as el


def test_angular_power_spectra(data_path):
    cls = el.photo.harmonic_space.angular_power_spectra(
        data_path / "angular_power_spectra.fits"
    )

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
    mms = el.photo.harmonic_space.mixing_matrices(data_path / "mixing_matrices.fits")

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

def test_covariance_matrices(data_path):
    mms = el.photo.harmonic_space.covariance_matrices(data_path / "cov.fits")

    expected_keys = {
            ('POS', 'POS', 'POS', 'POS', 0, 0, 0, 0),
            ('POS', 'POS', 'POS', 'SHE', 0, 0, 0, 0),
            ('POS', 'POS', 'SHE', 'SHE', 0, 0, 0, 0),
            ('POS', 'POS', 'POS', 'POS', 0, 0, 0, 1),
            ('POS', 'POS', 'POS', 'SHE', 0, 0, 0, 1),
            ('POS', 'POS', 'SHE', 'SHE', 0, 0, 0, 1),
            ('POS', 'POS', 'POS', 'SHE', 0, 0, 1, 0),
            ('POS', 'POS', 'POS', 'POS', 0, 0, 1, 1),
            ('POS', 'POS', 'POS', 'SHE', 0, 0, 1, 1),
            ('POS', 'POS', 'SHE', 'SHE', 0, 0, 1, 1),
            ('POS', 'SHE', 'POS', 'SHE', 0, 0, 0, 0),
            ('POS', 'SHE', 'SHE', 'SHE', 0, 0, 0, 0),
            ('POS', 'SHE', 'POS', 'POS', 0, 0, 0, 1),
            ('POS', 'SHE', 'POS', 'SHE', 0, 0, 0, 1),
            ('POS', 'SHE', 'SHE', 'SHE', 0, 0, 0, 1),
            ('POS', 'SHE', 'POS', 'SHE', 0, 0, 1, 0),
            ('POS', 'SHE', 'POS', 'POS', 0, 0, 1, 1),
            ('POS', 'SHE', 'POS', 'SHE', 0, 0, 1, 1),
            ('POS', 'SHE', 'SHE', 'SHE', 0, 0, 1, 1),
            ('SHE', 'SHE', 'SHE', 'SHE', 0, 0, 0, 0),
            ('SHE', 'SHE', 'POS', 'POS', 0, 0, 0, 1),
            ('SHE', 'SHE', 'POS', 'SHE', 0, 0, 0, 1),
            ('SHE', 'SHE', 'SHE', 'SHE', 0, 0, 0, 1),
            ('SHE', 'SHE', 'POS', 'SHE', 0, 0, 1, 0),
            ('SHE', 'SHE', 'POS', 'POS', 0, 0, 1, 1),
            ('SHE', 'SHE', 'POS', 'SHE', 0, 0, 1, 1),
            ('SHE', 'SHE', 'SHE', 'SHE', 0, 0, 1, 1),
            ('POS', 'POS', 'POS', 'POS', 0, 1, 0, 1),
            ('POS', 'POS', 'POS', 'SHE', 0, 1, 0, 1),
            ('POS', 'POS', 'SHE', 'SHE', 0, 1, 0, 1),
            ('POS', 'POS', 'POS', 'SHE', 0, 1, 1, 0),
            ('POS', 'POS', 'POS', 'POS', 0, 1, 1, 1),
            ('POS', 'POS', 'POS', 'SHE', 0, 1, 1, 1),
            ('POS', 'POS', 'SHE', 'SHE', 0, 1, 1, 1),
            ('POS', 'SHE', 'POS', 'SHE', 0, 1, 0, 1),
            ('POS', 'SHE', 'SHE', 'SHE', 0, 1, 0, 1),
            ('POS', 'SHE', 'POS', 'SHE', 0, 1, 1, 0),
            ('POS', 'SHE', 'POS', 'POS', 0, 1, 1, 1),
            ('POS', 'SHE', 'POS', 'SHE', 0, 1, 1, 1),
            ('POS', 'SHE', 'SHE', 'SHE', 0, 1, 1, 1),
            ('SHE', 'SHE', 'SHE', 'SHE', 0, 1, 0, 1),
            ('SHE', 'SHE', 'POS', 'SHE', 0, 1, 1, 0),
            ('SHE', 'SHE', 'POS', 'POS', 0, 1, 1, 1),
            ('SHE', 'SHE', 'POS', 'SHE', 0, 1, 1, 1),
            ('SHE', 'SHE', 'SHE', 'SHE', 0, 1, 1, 1),
            ('POS', 'SHE', 'POS', 'SHE', 1, 0, 1, 0),
            ('POS', 'SHE', 'POS', 'POS', 1, 0, 1, 1),
            ('POS', 'SHE', 'POS', 'SHE', 1, 0, 1, 1),
            ('POS', 'SHE', 'SHE', 'SHE', 1, 0, 1, 1),
            ('POS', 'POS', 'POS', 'POS', 1, 1, 1, 1),
            ('POS', 'POS', 'POS', 'SHE', 1, 1, 1, 1),
            ('POS', 'POS', 'SHE', 'SHE', 1, 1, 1, 1),
            ('POS', 'SHE', 'POS', 'SHE', 1, 1, 1, 1),
            ('POS', 'SHE', 'SHE', 'SHE', 1, 1, 1, 1),
            ('SHE', 'SHE', 'SHE', 'SHE', 1, 1, 1, 1),
    }

    assert set(mms.keys()) == expected_keys
