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


def test_writer_decorator_sets_history(monkeypatch):
    called = {}

    def fake_writer(func):
        def wrapper(*args, **kwargs):
            called["called"] = True
            return func(*args, **kwargs)

        wrapper._history = "test history"
        return wrapper

    monkeypatch.setattr(el.le3.pk_wl.angular_power_spectra, "write", fake_writer)

    @el.le3.pk_wl.angular_power_spectra.write
    def dummy_writer(path, results):
        return "ok"

    assert dummy_writer("foo", {}) == "ok"
    assert called["called"]
    assert getattr(dummy_writer, "_history", None) == "test history"
