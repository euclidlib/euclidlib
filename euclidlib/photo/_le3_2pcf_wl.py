from __future__ import annotations

import re
from os import PathLike
from typing import TYPE_CHECKING, Tuple, Union

import fitsio  # type: ignore [import-not-found]
import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence
    from typing import Any, TypeAlias


# type alias
_DictKey: "TypeAlias" = Union[str, int, Tuple["_DictKey", ...]]

def _key_from_string(s: str) -> _DictKey:
    """return key from string representation"""
    tuple_key_dict = None
    if 'POSPOS' in s:
        name_split = s.split('_')
        tuple_key_dict = ('POS', int(name_split[-2]), int(name_split[-1]))
    if 'POSSHEAR' in s:
        name_split = s.split('_')
        tuple_key_dict = ('POSSHEAR', name_split[-2], name_split[-1])
    else:
        pass
    return tuple_key_dict


def correlation_functions(
    path: str | PathLike[str],
    *,
    include: Sequence[tuple[Any, ...]] | None = None,
    exclude: Sequence[tuple[Any, ...]] | None = None,
) -> dict[_DictKey, NDArray[Any]]:
    """
    Read correlation functions from a Euclid data product.
    """

    xi = {}
    with fitsio.FITS(path) as fits:
        for hdu in fits:
            if not '2D' in hdu.get_extname():
                continue
            tuple_key = _key_from_string(hdu.get_extname())
            data = hdu.read()
            xi[tuple_key] = data
    return xi

def bandpowers(
    path: str | PathLike[str],
    *,
    include: Sequence[tuple[Any, ...]] | None = None,
    exclude: Sequence[tuple[Any, ...]] | None = None,
) -> dict[_DictKey, NDArray[Any]]:
    """
    Read bandpowers from a Euclid data product.
    """

    bandp = {}
    with fitsio.FITS(path) as fits:
        for hdu in fits:
            if not '2D' in hdu.get_extname():
                continue
            tuple_key = _key_from_string(hdu.get_extname())
            data = hdu.read()
            bandp[tuple_key] = data
    return bandp
