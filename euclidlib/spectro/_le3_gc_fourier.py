from __future__ import annotations
from warnings import warn
from typing import Any, Dict, Tuple, Iterable, Union
from os import PathLike
import fitsio
from numpy.typing import NDArray
from ._utils import _verify_input_file, _get_hdu_header, _get_hdu_data
from ._datamodel import Result

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
    # Multiply-used error message
    pk_not_found_message = "Invalid fits file provided, cannot locate power spectrum data."
    # File quality check
    _verify_input_file(path)
    # Read file and verify it contains information on the pk
    with fitsio.FITS(path) as fits_input:
        header = _get_hdu_header(fits_input[1])
        if "SPECTRUM" not in header["EXTNAME"]:
           raise ValueError(pk_not_found_message)
        elif header["EXTNAME"] == "SPECTRUM":
            data = _get_hdu_data(fits_input[1])
        # If a bispectrum file was provided, retrieve pk from next hdu and warn the user
        elif header["EXTNAME"] == "BISPECTRUM":
            warn(
                "\n\033[0;31m[!]\033[0m This looks to be a bispectrum file. "
                "Falling back to next HDU looking for data on the power spectrum."
            )
            _verify_input_file(path, check_extra_hdu=True)
            header = _get_hdu_header(fits_input[2])
            if header["EXTNAME"] != "SPECTRUM":
                raise ValueError(pk_not_found_message)
            data = _get_hdu_data(fits_input[2])
        else:
            raise ValueError(pk_not_found_message)
    return header, data

def power_spectrum(
        path: Union[str, PathLike[str]]
    ) -> Result:
    """
    Returns power spectrum data in the cloe-compatible euclidlib data format

    Parameters
    ----------
    path : str | PathLike
        path to the pk file

    Returns
    -------
    pk_result : Result
        class containing the pk data

    Examples
    --------
    >>> filename = "pk_z0.9-1.1.fits"
    >>> pk = power_spectrum(filename)
    >>> pk.k # returns bin centers
    >>> pk.k_eff # returns effective bin values
    >>> pk.mode_number # returns number of modes in each bin
    >>> pk.p[l] # returns multipole l (l = 0, ..., 4)
    >>> pk[l] # equivalent to above, returns multipole l
    >>> len(pk) # returns length of k array
    """
    header, data = _read_power_spectrum(path)
    pk_result = Result(
        data["K"],
        data["K_EFF"],
        data["NUM_MOD"],
        {l: data["PK{}".format(str(l))] for l in range(5)},
        header
    )
    return pk_result
    

