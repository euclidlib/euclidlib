import pytest  # type: ignore


@pytest.fixture
def data_path(request):
    """
    Return data path relative to test being run.
    """
    print(request.path.parent / "data")
    return request.path.parent / "data"
