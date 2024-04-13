import pytest


@pytest.fixture
def data_path(request):
    """
    Return data path relative to test being run.
    """
    return request.path.parent / "data"
