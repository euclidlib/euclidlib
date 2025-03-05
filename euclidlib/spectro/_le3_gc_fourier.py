from __future__ import annotations
from typing import Any, Dict, Tuple, Union
from os import PathLike
import numpy as np
import fitsio
from numpy.typing import NDArray
from ._utils import _get_hdu_header, _get_hdu_data

def _read_powerspectrum(path: Union[str, PathLike[str]]) -> Tuple[Dict[str, Any], NDArray[Any]]:
    """
    Reads power spectrum from Euclid LE3-GC fits file

    Parameters
    ----------
    filename : str | PathLike
        path to the fits file

    Returns
    -------
    header, data : Dict, NDArray
    """
    with fitsio.FITS(path) as fits_input:
        header = _get_powerspectrum_header(fits_input[1])
        data = _get_powerspectrum_data(fits_input[1])
    return header, data
