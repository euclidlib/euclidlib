from __future__ import annotations
from typing import Any, Dict, Tuple, Union
from os import PathLike
import numpy as np
import fitsio
from numpy.typing import NDArray
from ._utils import _verify_input_file, _get_hdu_header, _get_hdu_data

def _read_power_spectrum(path: Union[str, PathLike[str]]) -> Tuple[Dict[str, Any], NDArray[Any]]:
    """
    Reads power spectrum from Euclid LE3-GC fits file, returning both header and data

    Parameters
    ----------
    filename : str | PathLike
        path to the fits file

    Returns
    -------
    header, data : Dict, NDArray
    """
    # File quality check
    _verify_input_file(path)
    # Read file and verify it contains information on the pk
    with fitsio.FITS(path) as fits_input:
        header = _get_hdu_header(fits_input[1])
        if ("EXTNAME" not in header) or ("SPECTRUM" not in header["EXTNAME"]):
            raise ValueError(
                "Invalid fits file provided, cannot locate Power Spectrum data."
            )
        elif header["EXTNAME"] == "SPECTRUM":
            data = _get_hdu_data(fits_input[1])
        # If Bispectrum file provided, retrieve pk from next hdu (yet to be tested!)
        elif header["EXTNAME"] == "BISPECTRUM":
            header = _get_hdu_header(fits_input[2])
            data = _get_hdu_data(fits_input[2])
    return header, data
