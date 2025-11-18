from __future__ import annotations
from warnings import warn
from typing import Any, Dict, Tuple, Iterable, Union
from os import PathLike
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
        path: Union[str, PathLike[str], Iterable[str], Iterable[PathLike[str]]]
    ) -> Dict[Tuple[str, str, int, int], NDArray]:
    """
    Returns power spectrum data in the cloe-compatible euclidlib data format

    Parameters
    ----------
    path : str | PathLike | Iterable[str] | Iterable[PathLike]
        path(s) to the pk file(s) (see notes for more info on providing more than one file)

    Returns
    -------
    pk_dict : dict[NDArray]
        dictionaty containing the pk data

    Notes
    -----
    Multiple-file input is intended to load multiple redshift bins, not any collection of measurements.
    When passing multiple redshift bins, provide them as a sequence of filenames ordered by growing value of z.

    Examples
    --------
    >>> filename = "pk_z0.9-1.1.fits"
    >>> pk = power_spectrum(filename)
    >>> pk.keys() # returns dict_keys([('POS', 'POS', 1, 1)])
    >>> filenames = ["pk_z0.9-1.1.fits", "pk_1.1-1.3.fits"]
    >>> pk = power_spectrum(filenames)
    >>> pk.keys() # returns dict_keys([('POS', 'POS', 1, 1), ('POS', POS', 2, 2)])
    """
    if isinstance(path, Iterable) and (len(path) == 0):
        raise TypeError("Invalid input provided, it cannot be an empty sequence.")
    if (not isinstance(path, Iterable)) or isinstance(path, str):
        path = (path,)
    pk_dict = dict()
    for n, filename in enumerate(path):
        # TODO
        # For now, header data is ignored.
        # It could be necessary to read the effective z for each bin...
        _, pk = _read_power_spectrum(filename)
        pk_dict["POS", "POS", (n+1), (n+1)] = pk
    return pk_dict
    

