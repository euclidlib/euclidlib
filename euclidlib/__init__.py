__all__ = [
    "__version__",
    "__version_tuple__",
    "le3",
    "phz",
]

# generated version information
try:
    from ._version import __version__, __version_tuple__  # type: ignore [import-not-found, unused-ignore]
except ModuleNotFoundError:
    pass

# import submodules here so everything is available with package import
from . import le3, phz
