from __future__ import annotations

from collections.abc import Sequence
from os import PathLike
from typing import Any, Tuple, TypeAlias

import fitsio  # type: ignore [import-not-found]
import numpy as np
from numpy.typing import NDArray


# type alias
TwoPointKey: TypeAlias = Tuple[str, str, int, int]


def toc_match(
    key: tuple[Any, ...],
    include: Sequence[tuple[Any, ...]] | None = None,
    exclude: Sequence[tuple[Any, ...]] | None = None,
) -> bool:
    """return whether a dict key matches include/exclude criteria"""
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


def _read_metadata(hdu: fitsio.TableHDU) -> dict[str, Any]:
    """read array metadata from FITS HDU"""
    h = hdu.read_header()
    md = {}
    for key in h:
        if key.startswith("META "):
            md[key[5:].lower()] = h[key]
    return md


def _read_twopoint(fits: fitsio.FITS, ext: str) -> NDArray[Any]:
    """read two-point data from FITS"""
    # read data from extension
    arr: NDArray[Any] = fits[ext].read()
    # attach metadata to dtype
    dt = np.dtype(arr.dtype, metadata=_read_metadata(fits[ext]))
    # return a copy (unfortunately) with changed dtype
    return arr.astype(dt, casting="no", copy=True)


def angular_power_spectra(
    path: str | PathLike[str],
    *,
    include: Sequence[tuple[Any, ...]] | None = None,
    exclude: Sequence[tuple[Any, ...]] | None = None,
) -> dict[TwoPointKey, NDArray[Any]]:
    """
    Read angular power spectra from a Euclid data product.
    """

    # the returned set of cls
    cls = {}

    # open the FITS file for reading
    with fitsio.FITS(path) as fits:
        # get the TOC from the FITS file
        fits_toc = fits["CLTOC"].read()

        # read every entry in the TOC, add it to the list, then read the cls
        for entry in fits_toc:
            ext, k1, k2, i1, i2 = entry[["EXT", "NAME1", "NAME2", "BIN1", "BIN2"]]

            # skip if not selected
            if not toc_match((k1, k2, i1, i2), include=include, exclude=exclude):
                continue

            # read the cl from the extension and store in set of cls
            cls[k1, k2, i1, i2] = _read_twopoint(fits, ext)

    # return the dictionary of cls
    return cls


def mixing_matrices(
    path: str | PathLike[str],
    *,
    include: Sequence[tuple[Any, ...]] | None = None,
    exclude: Sequence[tuple[Any, ...]] | None = None,
) -> dict[TwoPointKey, NDArray[Any]]:
    """
    Read mixing matrices from a Euclid data product.
    """

    # the returned set of mms
    mms = {}

    # open the FITS file for reading
    with fitsio.FITS(path) as fits:
        # get the TOC from the FITS file
        fits_toc = fits["MMTOC"].read()

        # read every entry in the TOC, add it to the list, then read the mms
        for entry in fits_toc:
            ext, k1, k2, i1, i2 = entry[["EXT", "NAME1", "NAME2", "BIN1", "BIN2"]]

            # skip if not selected
            if not toc_match((k1, k2, i1, i2), include=include, exclude=exclude):
                continue

            # read the mixing matrix from the extension and store in set of mms
            mms[k1, k2, i1, i2] = _read_twopoint(fits, ext)

    # return the dictionary of mms
    return mms
