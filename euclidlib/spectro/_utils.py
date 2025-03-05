from __future__ import annotations
from typing import Any, Dict
import fitsio
from numpy.typing import NDArray

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
