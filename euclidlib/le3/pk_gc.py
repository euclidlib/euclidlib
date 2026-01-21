from __future__ import annotations
from warnings import warn
from typing import Any, Dict, Tuple, Union
from os import PathLike
import fitsio
import os
from numpy.typing import NDArray

from cosmolib.data import PowerSpectrum  # type: ignore [import-not-found]


def _verify_input_file(
    path: Union[str, PathLike[str]], check_extra_hdu: bool = False
) -> None:
    """
    Verifies that input file is a compatible LE3-GC fits file

    Parameters
    ----------
    path : str | PathLike
        path to the (hopefully) fits file
    check_extra_hdu : bool
        verify file contains information on two different hdu (e.g. bispectrum and power spectrum)

    Raises
    ------
    TypeError
        if `path` is not a valid string or PathLike object
    FileNotFoundError
        if file does not exist
    ValueError
        if file is not a fits file or is it a fits file with an incompatible structure
    """
    hdul_target_length = 3 if check_extra_hdu else 2
    if not any((isinstance(path, str), isinstance(path, PathLike))):
        raise TypeError(
            "Provided fits file name must be a string or a PathLike object."
        )
    if not os.path.isfile(path):
        raise FileNotFoundError("Could not find file '{}'.".format(str(path)))
    try:
        with fitsio.FITS(path) as fits_input:
            assert len(fits_input) > hdul_target_length - 1
            assert "EXTNAME" in _get_hdu_header(fits_input[hdul_target_length - 1])
    except OSError:
        raise ValueError("Provided file is not a valid fits file.")
    except AssertionError:
        raise ValueError(
            "Provided fits file does not match the structure of a valid LE3-GC product."
        )
    except Exception:
        raise RuntimeError("Invalid file provided.")


def _get_hdu_header(hdu: fitsio.TableHDU) -> Dict[str, Any]:
    """
    Reads header from fits file HDU

    Parameters
    ----------
    hdu : fitsio.TableHDU
        HDU from fits table

    Returns
    -------
    head : Dict
        dictionary containing all header entries
    """
    head: Dict[str, Any] = dict(hdu.read_header())
    return head


def _get_hdu_data(hdu: fitsio.TableHDU) -> NDArray[Any]:
    """
    Reads data from fits file HDU

    Parameters
    ----------
    hdu : fitsio.TableHDU
        HDU from fits table

    Returns
    -------
    arr : NDArray
        record array containing information stored in the hdu
    """
    arr: NDArray[Any] = hdu.read()
    return arr


def _read_power_spectrum(
    path: Union[str, PathLike[str]],
) -> Tuple[Dict[str, Any], NDArray[Any]]:
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
    pk_not_found_message = (
        "Invalid fits file provided, cannot locate power spectrum data."
    )
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


def power_spectrum(path: Union[str, PathLike[str]]) -> PowerSpectrum:
    """
    Returns power spectrum data in the cloe-compatible euclidlib data format

    Parameters
    ----------
    path : str | PathLike
        path to the pk file

    Returns
    -------
    pk_result : PowerSpectrum
        class containing the pk data
    """
    header, data = _read_power_spectrum(path)
    pk_result = PowerSpectrum(
        data["K"],
        data["K_EFF"],
        data["NUM_MOD"],
        {l: data["PK{}".format(str(l))] for l in range(5)},
        0,  # Fiducial cosmology is not stored in the fits file
        0,  # redshift_eff is not stored in the fits file
        1.0 / header["SN_VALUE"],
        header["SN_VALUE"],
    )
    return pk_result
