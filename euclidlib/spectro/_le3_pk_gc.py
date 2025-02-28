from __future__ import annotations
from typing import Any, Dict, Tuple, Union
from os import PathLike
import numpy as np
import fitsio
from numpy.typing import NDArray

def _get_powerspectrum_data(hdu: fitsio.TableHDU) -> NDArray[Any]:
    """
    Reads power spectrum data from fits file HDU

    Parameters
    ----------
    hdu : fitsio.TableHDU
        HDU from fits table (for any typical LE3-GC-PK file, the second in the list obtained when reading the fits file)

    Returns
    -------
    pk_arr : NDArray
        record array containing information on the power spectrum

    Notes
    -----
    The columns in the pk record array are: "K", "K_EFF" (for now, it will just display the same values of "K"), "PKL" where "L" goes from "0" to "4" and "NUM_MOD".
    """
    pk_arr: NDArray[Any] = hdu.read()
    return pk_arr

def _get_powerspectrum_header(hdu: fitsio.TableHDU) -> Dict[str, Any]:
    """
    Reads power spectrum header from fits file HDU

    Parameters
    ----------
    hdu : fitsio.TableHDU
        HDU from fits table (for any typical LE3-GC-PK file, the second in the list obtained when reading the fits file)

    Returns
    -------
    pk_head : Dict
        dictionary containing all header entries
    """
    pk_head: Dict[str, Any] = dict(hdu.read_header())
    return pk_head

def read_powerspectrum(path: Union[str, PathLike[str]]) -> Tuple[Dict[str, Any], NDArray[Any]]:
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
