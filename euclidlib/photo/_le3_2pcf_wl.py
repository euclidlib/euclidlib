from __future__ import annotations

from os import PathLike
from typing import TYPE_CHECKING, Tuple, Union, Any

import fitsio  # type: ignore [import-not-found]
from numpy.typing import NDArray # type: ignore


# type alias
_DictKey: "TypeAlias" = Union[str, int, Tuple["_DictKey", ...]]

def _key_from_string(s: str) -> _DictKey:
    """return key from string representation"""
    tuple_key_dict = None
    if 'POSPOS' in s:
        name_split = s.split('_')
        tuple_key_dict = ('POS', 'POS', int(name_split[-2]), int(name_split[-1]))
    elif 'POSSHEAR' in s:
        name_split = s.split('_')
        tuple_key_dict = ('POS', 'SHE', int(name_split[-2]), int(name_split[-1]))
    elif 'SHEAR' in s:
        name_split = s.split('_')
        tuple_key_dict = ('SHE', 'SHE', int(name_split[-2]), int(name_split[-1]))
    else:
        pass
    return tuple_key_dict


def correlation_functions(
    path: str | PathLike[str],
) -> dict[_DictKey, NDArray[Any]]:
    """
    Reads 2D correlation functions from a Euclid data product.

    Parameters
    ----------
    path : str | PathLike[str]
        The path to the FITS file containing the Euclid data product.

    Returns
    -------
    dict[_DictKey, NDArray[Any]]
        A dictionary where each key is a tuple representing the HDU name (generated by 
        `_key_from_string`), and each value is a NumPy array containing the data of 
        the corresponding HDU.

    Notes
    -----
    - Only HDUs whose extension name contains '2D' are considered.
    - The `fitsio.FITS` context manager is used to read the FITS file efficiently.
    - The function depends on `_key_from_string()` to convert the HDU names into keys 
      for the output dictionary.
    """

    xi = {}
    with fitsio.FITS(path) as fits:
        for hdu in fits:
            if '2D' not in hdu.get_extname():
                continue
            tuple_key = _key_from_string(hdu.get_extname())
            data = hdu.read()
            xi[tuple_key] = data
    return xi

def bandpowers(
    path: str | PathLike[str],
) -> dict[_DictKey, NDArray[Any]]:
    """
    Reads 2D bandpowers from a Euclid data product.

    Parameters
    ----------
    path : str | PathLike[str]
        The path to the FITS file containing the Euclid bandpower data product.

    Returns
    -------
    dict[_DictKey, NDArray[Any]]
        A dictionary where each key is a tuple representing the HDU name (generated by 
        `_key_from_string`), and each value is a NumPy array containing the data of 
        the corresponding HDU.

    Notes
    -----
    - Only HDUs whose extension name contains '2D' are considered.
    - The `fitsio.FITS` context manager is used to read the FITS file efficiently.
    - The function depends on `_key_from_string()` to convert the HDU names into keys 
      for the output dictionary.
    """
    bandp = {}
    with fitsio.FITS(path) as fits:
        for hdu in fits:
            if '2D' not in hdu.get_extname():
                continue
            if 'POSPOS' in hdu.get_extname():
                tuple_key = _key_from_string(hdu.get_extname())
                data = hdu.read()
                bandp[tuple_key] = data
                dtype_common = ([('L', '>f8'), ('CL', '>f8'), 
                                 ('LMIN', '>f8'), ('LMAX', '>f8')])
                bandp[tuple_key].dtype = dtype_common
            elif 'POSSHEAR' in hdu.get_extname():
                tuple_key = _key_from_string(hdu.get_extname())
                data = hdu.read()
                bandp[tuple_key] = data
                dtype_common = ([('L', '>f8'), ('CL_E', '>f8'), ('CL_B', '>f8'), 
                                 ('LMIN', '>f8'), ('LMAX', '>f8')])
                bandp[tuple_key].dtype = dtype_common
            elif 'SHEARSHEAR' in hdu.get_extname():
                tuple_key = _key_from_string(hdu.get_extname())
                data = hdu.read()
                bandp[tuple_key] = data
                dtype_common = ([('L', '>f8'), ('CL_E', '>f8'), ('CL_B', '>f8'), 
                                 ('LMIN', '>f8'), ('LMAX', '>f8')])
                bandp[tuple_key].dtype = dtype_common

    return bandp

def cosebis(
    path: str | PathLike[str],
) -> dict[_DictKey, NDArray[Any]]:
    """
    Reads 2D bandpowers from a Euclid data product.

    Parameters
    ----------
    path : str | PathLike[str]
        The path to the FITS file containing the Euclid bandpower data product.

    Returns
    -------
    dict[_DictKey, NDArray[Any]]
        A dictionary where each key is a tuple representing the HDU name (generated by 
        `_key_from_string`), and each value is a NumPy array containing the data of 
        the corresponding HDU.

    Notes
    -----
    - Only HDUs whose extension name contains '2D' are considered.
    - The `fitsio.FITS` context manager is used to read the FITS file efficiently.
    - The function depends on `_key_from_string()` to convert the HDU names into keys 
      for the output dictionary.
    """
    
    cb = {}
    with fitsio.FITS(path) as fits:
        for hdu in fits:
            if '2D' not in hdu.get_extname():
                continue
            tuple_key = _key_from_string(hdu.get_extname())
            data = hdu.read()
            cb[tuple_key] = data
    return cb
