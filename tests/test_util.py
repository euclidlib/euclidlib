import euclidlib as el


def test_writable():
    def func():
        pass

    func.__module__ = "mymodule"

    @el._util.writer(func)
    def writer():
        pass

    assert func.write is writer
    assert func.write.__module__ == "mymodule"
    assert func.write.__name__ == "write"
    assert func.write.__qualname__ == f"{func.__qualname__}.write"
