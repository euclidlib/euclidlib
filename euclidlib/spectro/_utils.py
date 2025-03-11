from __future__ import annotations
from typing import Any, Dict, Union
import os
from os import PathLike
import fitsio
from numpy.typing import NDArray

def _verify_input_file(path: Union[str, PathLike[str]]) -> None:
    """
    Verifies that input file is a compatible LE3-GC fits file

    Parameters
    ----------
    path : str | PathLike 
        path to the (hopefully) fits file

    Raises
    ------
    TypeError
        if `path` is not a valid string or PathLike object
    FileNotFoundError
        if file does not exist
    ValueError
        if file is not a fits file or is it a fits file with an incompatible structure
    """
    if not any((isinstance(path, str), isinstance(path, PathLike))):
        raise TypeError("Provided fits file name must be a string or a PathLike object.")
    if not os.path.isfile(path):
        raise FileNotFoundError("Could not find file '{}'.".format(str(path)))
    try:
        with fitsio.FITS(path) as fits_input:
            assert len(fits_input) > 1
            assert "EXTNAME" in _get_hdu_header(fits_input[1])
    except OSError:
        raise ValueError(
            "Provided file is not a valid fits file."
        )
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
