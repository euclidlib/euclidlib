from __future__ import annotations

import re
from os import PathLike
from typing import TYPE_CHECKING, Tuple, Union

import fitsio  # type: ignore [import-not-found]
import numpy as np # type: ignore
from numpy.typing import NDArray, DTypeLike # type: ignore

from dataclasses import dataclass

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence
    from typing import Any, TypeAlias # type: ignore


# type alias
_DictKey: "TypeAlias" = Union[str, int, Tuple["_DictKey", ...]]


def normalize_result_axis(axis, result, ell):
    """Return an axis tuple for a result."""
    try:
        from numpy.lib.array_utils import normalize_axis_tuple # type: ignore
    except ModuleNotFoundError:
        from numpy.lib.stride_tricks import normalize_axis_tuple # type: ignore

    if axis is None:
        if result.ndim == 0:
            axis = ()
        elif isinstance(ell, tuple):
            axis = tuple(range(-len(ell), 0))
        else:
            axis = -1
    return normalize_axis_tuple(axis, result.ndim, "axis")

@dataclass(frozen=True, repr=False)
class Result:
    """
    Container for results.
    """
    array: NDArray[Any]
    ell: NDArray[Any] | tuple[NDArray[Any], ...] | None = None
    axis: int | tuple[int, ...] | None = None
    lower: NDArray[Any] | tuple[NDArray[Any], ...] | None = None
    upper: NDArray[Any] | tuple[NDArray[Any], ...] | None = None
    weight: NDArray[Any] | tuple[NDArray[Any], ...] | None = None

    def __post_init__(self) -> None:
        # Ensure array is of float dtype
        float_array = np.asarray(self.array, dtype=float)
        object.__setattr__(self, "array", float_array)

        # Normalize the axis after setting the array
        axis = normalize_result_axis(self.axis, self.array, self.ell)
        object.__setattr__(self, "axis", axis)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(axis={self.axis!r})"

    def __array__(self, dtype=None, *, copy=None) -> NDArray[Any]:
        if copy is not None:
            return self.array.__array__(dtype, copy=copy)
        return self.array.__array__(dtype)

    def __getitem__(self, key):
        return self.array[key]

    @property
    def ndim(self) -> int:
        return self.array.ndim

    @property
    def shape(self) -> tuple[int, ...]:
        return self.array.shape

    @property
    def dtype(self) -> DTypeLike:
        return self.array.dtype

def _key_from_string(s: str) -> _DictKey:
    """
    Return key for a given string representation.
    """
    parts = re.split(r"(?<!\\)-", s.replace("\\\\", "\0"))
    if len(parts) > 1:
        return tuple(map(_key_from_string, parts))
    key = parts[0]
    key = key.replace("\\-", "-")
    key = key.replace("\0", "\\")
    return int(key) if key.removeprefix("-").isdigit() else key

def _read_metadata(hdu):
    """read array metadata from FITS HDU"""
    h = hdu.read_header()
    md = {}
    for key in h:
        if key.startswith("META "):
            md[key[5:].lower()] = h[key]
    return md

def _read_result(hdu):
    """
    Read a result array from FITS.
    """

    from ast import literal_eval

    # read columnar data from extension
    data = hdu.read()
    h = hdu.read_header()

    # the angular axis
    axis = literal_eval(h["ELLAXIS"])

    # get data array and move axis back to right position
    arr = np.moveaxis(data["ARRAY"], tuple(range(len(axis))), axis)

    # sort ell axes into natural order
    order = np.argsort(axis)

    # get ells
    _ell = data["ELL"]
    if _ell.ndim == 1:
        ell = _ell
    else:
        ell = tuple(_ell[: arr.shape[axis[i]], i] for i in order)

    # get lower bounds
    _lower = data["LOWER"]
    if _lower.ndim == 1:
        lower = _lower
    else:
        lower = tuple(_lower[: arr.shape[axis[i]], i] for i in order)

    # get upper bounds
    _upper = data["UPPER"]
    if _upper.ndim == 1:
        upper = _upper
    else:
        upper = tuple(_upper[: arr.shape[axis[i]], i] for i in order)

    # get weights
    _weight = data["WEIGHT"]
    if _weight.ndim == 1:
        weight = _weight
    else:
        weight = tuple(_weight[: arr.shape[axis[i]], i] for i in order)

    # construct result array with ancillary arrays and metadata
    return Result(
        arr.view(np.dtype(arr.dtype, metadata=_read_metadata(hdu))),
        axis=tuple(axis[i] for i in order),
        ell=ell,
        lower=lower,
        upper=upper,
        weight=weight,
    )

def read(path) -> dict[_DictKey, NDArray[Any]]:
    """
    Read a set of results from a FITS file.
    """

    # the returned set of cls
    results = {}

    # read all HDUs in file into dict keys
    with fitsio.FITS(path) as fits:
        for hdu in fits:
            if not hdu.has_data():
                continue
            ext = hdu.get_extname()
            if not ext:
                continue
            key = _key_from_string(ext)
            if not key:
                continue
            results[key] = _read_result(hdu)

    return results

def angular_power_spectra(path) -> dict[_DictKey, NDArray[Any]]:
    """
    Read a set of results from a FITS file.
    """

    # the returned set of cls
    results = read(path)
    return results

def mixing_matrices(path) -> dict[_DictKey, NDArray[Any]]:
    """
    Read a set of results from a FITS file.
    """

    # the returned set of mixing matrices
    results = read(path)
    results
    return results
