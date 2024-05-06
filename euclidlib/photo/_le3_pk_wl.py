from __future__ import annotations

import re
from os import PathLike, path
from typing import TYPE_CHECKING, Tuple, Union

import fitsio  # type: ignore [import-not-found]
import numpy as np
from numpy.typing import NDArray

from euclidlib.photo import photo_data

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence
    from typing import Any, TypeAlias


# type alias
_DictKey: "TypeAlias" = Union[str, int, Tuple["_DictKey", ...]]


def toc_match(
    key: _DictKey,
    include: Sequence[tuple[Any, ...]] | None = None,
    exclude: Sequence[tuple[Any, ...]] | None = None,
) -> bool:
    """return whether a dict key matches include/exclude criteria"""
    if not isinstance(key, tuple):
        key = (key,)
    if include is not None:
        for pattern in include:
            if all(p is Ellipsis or p == k for p, k in zip(pattern, key)):
                break
        else:
            return False
    if exclude is not None:
        for pattern in exclude:
            if all(p is Ellipsis or p == k for p, k in zip(pattern, key)):
                return False
    return True


def _key_from_string(s: str) -> _DictKey:
    """return key from string representation"""
    keys = s.split(";")
    if len(keys) > 1:
        return tuple(map(_key_from_string, keys))
    keys = keys[0].split(",")
    if len(keys) > 1:
        return tuple(map(_key_from_string, keys))
    key = keys[0]
    return int(key) if key.isdigit() else key


def _read_metadata(hdu: fitsio.TableHDU) -> dict[str, Any]:
    """read array metadata from FITS HDU"""
    h = hdu.read_header()
    md = {}
    for key in h:
        if key.startswith("META "):
            md[key[5:].lower()] = h[key]
    return md


def _read_twopoint(hdu: fitsio.TableHDU) -> NDArray[Any]:
    """read two-point data from FITS HDU"""
    # read data from extension
    arr: NDArray[Any] = hdu.read()
    # attach metadata to dtype
    dt = np.dtype(arr.dtype, metadata=_read_metadata(hdu))
    # return a copy (unfortunately) with changed dtype
    return arr.astype(dt, casting="no", copy=True)


def _iterfits(
    path: str | PathLike[str],
    tag: str,
    include: Sequence[tuple[Any, ...]] | None = None,
    exclude: Sequence[tuple[Any, ...]] | None = None,
) -> Iterator[fitsio.TableHDU]:
    """
    Iterate over HDUs that correspond to *tag* and have valid keys.
    """
    with fitsio.FITS(path) as fits:
        for hdu in fits:
            if not re.match(f"^{tag}\\d+$", hdu.get_extname()):
                continue
            h = hdu.read_header()
            s = h.get("DICTKEY")
            if s is None:
                continue
            key = _key_from_string(s)
            if not toc_match(key, include=include, exclude=exclude):
                continue
            yield key, hdu


def angular_power_spectra(
    path: str | PathLike[str],
    *,
    include: Sequence[tuple[Any, ...]] | None = None,
    exclude: Sequence[tuple[Any, ...]] | None = None,
) -> dict[_DictKey, NDArray[Any]]:
    """
    Read angular power spectra from a Euclid data product.
    """

    cls = {}
    for key, hdu in _iterfits(path, "CL", include, exclude):
        cls[key] = _read_twopoint(hdu)
    return cls

def xi_tpcf(
    path: str | PathLike[str],
    *,
    include: Sequence[tuple[Any, ...]] | None = None,
    exclude: Sequence[tuple[Any, ...]] | None = None,
) -> dict[_DictKey, NDArray[Any]]:
    """
    Read the Weak Lensing two-point correlation functions for LE3.
    """

    corr_str = "ShearShear"

    # Load LE3 output file: 2PCF
    tpcf_le3 = photo_data.TpcfDataReader.from_fits(
        path,
        corr_str=corr_str,
    )

    z_combinations = list(tpcf_le3._data[corr_str].keys())

    # Get LE3 data
    xis = {}

    for component in ("+", "-"):
        xis[component] = {}
        for i, z_comb in enumerate(z_combinations):
            xis[component][z_comb]= tpcf_le3.get_xi_component(component,
                                                             z_comb)
    return xis

def mixing_matrices(
    path: str | PathLike[str],
    *,
    include: Sequence[tuple[Any, ...]] | None = None,
    exclude: Sequence[tuple[Any, ...]] | None = None,
) -> dict[_DictKey, NDArray[Any]]:
    """
    Read mixing matrices from a Euclid data product.
    """

    mms = {}
    for key, hdu in _iterfits(path, "MM", include, exclude):
        mms[key] = _read_twopoint(hdu)
    return mms
