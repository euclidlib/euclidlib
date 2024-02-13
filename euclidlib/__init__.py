__all__ = [
    "__version__",
    "__version_tuple__",
    "photo",
    "spectro",
]

# generated version information
try:
    from ._version import __version__, __version_tuple__
except ModuleNotFoundError:
    pass

# import submodules here so everything is available with package import
from . import photo, spectro
