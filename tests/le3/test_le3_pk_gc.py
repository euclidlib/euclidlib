import euclidlib as el  # Replace with the actual import


def test_pk_gc_with_dummy_data(data_path):
    # Load the dummy data from the FITS file
    data = el.le3.pk_gc.power_spectrum(
        data_path / "dummy_pk.fits"
    )  # Adjust the function as necessary

    # Example assertions based on expected output
    assert data is not None, "Data should not be None"
    assert isinstance(data, dict), "Data should be a dictionary"
    assert ("SPE", "SPE", 0, 0) in data, "Data should contain ('SPE', 'SPE', 0, 0)"
    assert data[("SPE", "SPE", 0, 0)] is not None, (
        "Power spectrum data should not be None"
    )
