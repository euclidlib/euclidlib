from __future__ import annotations

import os
from os import PathLike
import numpy as np
import fitsio  # type: ignore [import-not-found]
from numpy.typing import NDArray
from .._util import writer

from cosmolib.data import TwoPointCorrelationFunction, COSEBI  # type: ignore [import-not-found]


TYPE_CHECKING = True
if TYPE_CHECKING:
    from typing import Any

    try:
        from typing import TypeAlias
    except ImportError:
        from typing_extensions import TypeAlias

    _DictKey: TypeAlias = tuple[str, str, int, int]


def _key_from_string(s: str) -> _DictKey | None:
    """
    Decode a Euclid-style EXTNAME into an internal dictionary key.

    Parameters
    ----------
    s : str
        The EXTNAME string from a FITS HDU, e.g. 'SHEARSHEAR2D_1_2'.

    Returns
    -------
    tuple[str, str, int, int] or None
        A tuple of the form (TYPE1, TYPE2, BIN1, BIN2), or ``None`` if the
        EXTNAME does not match any known correlation type.

    Notes
    -----
    - This function is shared by the correlation function, bandpower, and
      COSEBI readers.
    - Recognized EXTNAME prefixes: ``POSPOS``, ``POSSHEAR``, ``SHEAR``.
    """
    tuple_key_dict = None
    if "POSPOS" in s:
        name_split = s.split("_")
        tuple_key_dict = ("POS", "POS", int(name_split[-2]), int(name_split[-1]))
    elif "POSSHEAR" in s:
        name_split = s.split("_")
        tuple_key_dict = ("POS", "SHE", int(name_split[-2]), int(name_split[-1]))
    elif "SHEAR" in s:
        name_split = s.split("_")
        tuple_key_dict = ("SHE", "SHE", int(name_split[-2]), int(name_split[-1]))
    return tuple_key_dict


def correlation_functions(
    path: str | PathLike[str],
) -> dict[_DictKey, TwoPointCorrelationFunction]:
    """
    Read Euclid 2D two-point correlation functions from a FITS file.

    Parameters
    ----------
    path : str or PathLike
        Path to the FITS file containing the correlation-function data.

    Returns
    -------
    dict[_DictKey, TwoPointCorrelationFunction]
        A dictionary mapping keys of the form
        ``(TYPE1, TYPE2, BIN1, BIN2)`` to ``TwoPointCorrelationFunction`` objects.

    Notes
    -----
    - Only HDUs whose EXTNAME contains ``'2D'`` are read.
    - The mapping of FITS columns to TPCF components follows the official
      Euclid naming conventions:
        * SHE–SHE → (XI_P, XI_M, XI_X)
        * POS–SHE → (GAMMA_T, GAMMA_X)
        * POS–POS → (WTHETA)
    """

    xi: dict[_DictKey, TwoPointCorrelationFunction] = {}
    with fitsio.FITS(path) as fits:
        for hdu in fits[1:]:
            extname = hdu.get_extname()
            key = _key_from_string(extname)
            if key is None:
                raise ValueError("Encountered None for key")
            data = hdu.read()
            THETA = data["THETA"]
            WEIGHT = data["WEIGHT"]
            if "SHEARSHEAR" in extname:
                array = np.array(
                    [[data["XI_P"], data["XI_X"]], [data["XI_X"], data["XI_M"]]]
                )
                axis = (2,)
            elif "POSSHEAR" in extname:
                array = np.array([data["GAMMA_T"], data["GAMMA_X"]])
                axis = (1,)
            elif "POSPOS" in extname:
                array = np.array(data["WTHETA"])
                axis = (0,)
            else:
                raise ValueError(f"Unknown file type: {key}")

            xi[key] = TwoPointCorrelationFunction(
                array, theta=THETA, axis=axis, weight=WEIGHT
            )
    return xi


@writer(correlation_functions)
def _(
    path: str | PathLike[str], results: dict[_DictKey, TwoPointCorrelationFunction]
) -> None:
    """
    Write Euclid-style 2D two-point correlation functions to a FITS file.

    Parameters
    ----------
    path : str or PathLike
        Destination path for the output FITS file.
    results : dict
        Dictionary mapping keys of the form
        ``(TYPE1, TYPE2, BIN1, BIN2)`` to ``TwoPointCorrelationFunction`` objects.
    """

    if os.path.exists(path):
        os.remove(path)

    with fitsio.FITS(path, "rw") as fits:
        # ---- PRIMARY HDU ----
        hdr = {
            "THMIN": 1.0,
            "THMAX": 400.0,
            "THUNIT": "arcmin",
            "BINSCALE": "LOG",
            "COORDSYS": "SPHERICAL",
        }
        fits.write(None, header=hdr)

        # ---- EXTENSIONS ----
        for key, tpcf in results.items():
            # Key format is: (TYPE1, TYPE2, BIN1, BIN2)
            type1, type2, bin1, bin2 = key
            corr_type = (type1.upper(), type2.upper())

            # Build Euclid-style EXTNAME, e.g. "SHEARSHEAR2D_1_2"
            if corr_type == ("SHE", "SHE"):
                extname = f"SHEARSHEAR2D_{bin1}_{bin2}"
            elif corr_type == ("POS", "SHE"):
                extname = f"POSSHEAR2D_{bin1}_{bin2}"
            elif corr_type == ("POS", "POS"):
                extname = f"POSPOS2D_{bin1}_{bin2}"
            else:
                raise ValueError(f"Unknown correlation type from key: {key}")

            theta = np.asarray(tpcf.theta)
            weight = np.asarray(tpcf.weight)
            arr = np.asarray(tpcf.array)

            # ---- SHEAR–SHEAR ----
            if corr_type == ("SHE", "SHE"):
                # arr shape: (2,2,N)
                xi_p = arr[0, 0]
                xi_x = arr[0, 1]
                xi_m = arr[1, 1]

                dtype = [
                    ("THETA", "f8"),
                    ("XI_P", "f8"),
                    ("XI_M", "f8"),
                    ("XI_X", "f8"),
                    ("WEIGHT", "f8"),
                ]

                data = np.zeros(len(theta), dtype=dtype)
                data["THETA"] = theta
                data["XI_P"] = xi_p
                data["XI_M"] = xi_m
                data["XI_X"] = xi_x
                data["WEIGHT"] = weight

            # ---- GALAXY–SHEAR ----
            elif corr_type == ("POS", "SHE"):
                # arr shape: (2, N)
                gamma_t = arr[0]
                gamma_x = arr[1]

                dtype = [
                    ("THETA", "f8"),
                    ("GAMMA_T", "f8"),
                    ("GAMMA_X", "f8"),
                    ("WEIGHT", "f8"),
                ]

                data = np.zeros(len(theta), dtype=dtype)
                data["THETA"] = theta
                data["GAMMA_T"] = gamma_t
                data["GAMMA_X"] = gamma_x
                data["WEIGHT"] = weight

            # ---- GALAXY–GALAXY ----
            elif corr_type == ("POS", "POS"):
                # arr shape: (N,)
                wtheta = arr

                dtype = [
                    ("THETA", "f8"),
                    ("WTHETA", "f8"),
                    ("WEIGHT", "f8"),
                ]

                data = np.zeros(len(theta), dtype=dtype)
                data["THETA"] = theta
                data["WTHETA"] = wtheta
                data["WEIGHT"] = weight

            # ---- Extension header ----
            hdr = {
                "BIN_ID1": int(bin1),
                "BIN_ID2": int(bin2),
            }
            # Write the header to the FITS file
            history = getattr(_, "_history", None)
            if history is not None:
                hdr["HISTORY"] = history
            fits.write(data, extname=extname, header=hdr)


def bandpowers(path: str | PathLike[str]) -> dict[_DictKey, NDArray[Any]]:
    """
    [SOON TO BE DEPRECATED] Read Euclid 2D bandpowers from a FITS product.

    Parameters
    ----------
    path : str or PathLike
        Path to the FITS file containing the bandpower data.

    Returns
    -------
    dict[_DictKey, numpy.ndarray]
        Mapping from decoded EXTNAME keys to NumPy structured arrays containing
        the bandpower data.

    Notes
    -----
    - Only HDUs whose EXTNAME contains ``'2D'`` are read.
    - Column names are normalized to Euclid conventions:
        * POS–POS → (L, CL)
        * POS–SHE → (L, CL_E, CL_B)
        * SHE–SHE → (L, CL_E, CL_B)
    """
    bandp: dict[_DictKey, NDArray[Any]] = {}
    with fitsio.FITS(path) as fits:
        for hdu in fits:
            if "2D" not in hdu.get_extname():
                continue
            key = _key_from_string(hdu.get_extname())
            if key is None:
                continue
            data = hdu.read()
            if key[:2] == ("POS", "POS"):
                data.dtype.names = ["L", "CL", "LMIN", "LMAX"]
            elif key[:2] == ("POS", "SHE"):
                data.dtype.names = ["L", "CL_E", "CL_B", "LMIN", "LMAX"]
            elif key[:2] == ("SHE", "SHE"):
                data.dtype.names = ["L", "CL_E", "CL_B", "LMIN", "LMAX"]
            bandp[key] = data
    return bandp


def cosebis(path: str | PathLike[str]) -> dict[_DictKey, COSEBI]:
    """
    Read Euclid COSEBI bandpowers from a FITS product.

    Parameters
    ----------
    path : str or PathLike
        The path to the FITS file containing the COSEBI data product.

    Returns
    -------
    dict[_DictKey, COSEBI]
        Mapping from decoded EXTNAME keys to ``COSEBI`` dataclass instances.

    Notes
    -----
    - The PRIMARY HDU is expected to contain:
        ``THMIN``, ``THMAX``, ``NMODES``.
    - Each extension HDU must have EXTNAME formatted as:
        ``SHEARSHEAR2D_COSEBI_j_k``.
    - Columns read: MODE, EE, EB, BB (Euclid ordering).
    """

    raw = {}
    cb: dict[_DictKey, COSEBI] = {}

    with fitsio.FITS(path) as fits:
        # ---- Extract metadata from PRIMARY ----
        primary = fits[0].read_header()
        thmin = primary.get("THMIN")
        thmax = primary.get("THMAX")
        nmodes = primary.get("NMODES")

        # ---- Loop over extension HDUs ----
        for hdu in fits:
            if "2D" not in hdu.get_extname():
                continue

            key = _key_from_string(hdu.get_extname())
            if key is None:
                continue

            raw[key] = hdu.read()

            # Build float array (EE, EB, BB)
            arr = np.array(
                [
                    raw[key]["EE"],
                    raw[key]["EB"],
                    raw[key]["BB"],
                ]
            )

            # Create COSEBI dataclass
            cb[key] = COSEBI(
                array=arr,
                mode=raw[key]["MODE"].astype(int),
                thmin=thmin,
                thmax=thmax,
                nmodes=nmodes,
            )

    return cb


@writer(cosebis)
def _(path: str | PathLike[str], results: dict[_DictKey, COSEBI]) -> None:
    """
    Write Euclid COSEBI bandpowers to a FITS file.

    Parameters
    ----------
    path : str or PathLike
        Destination path for the output FITS file.
    results : dict
        Mapping from keys of the form ``(TYPE1, TYPE2, BIN1, BIN2)``
        to COSEBI dataclass objects.

    Notes
    -----
    - The PRIMARY HDU contains angular metadata:
        ``THMIN``, ``THMAX``, ``NMODES`` (taken from any COSEBI instance).
    - Only SHE–SHE COSEBIs are supported.
    - EXTNAME is set to ``SHEARSHEAR2D_COSEBI_j_k``.
    - Columns written per HDU:
        ``MODE``, ``EE``, ``BB``, ``EB``.
    """

    # Overwrite existing file
    if os.path.exists(path):
        os.remove(path)

    with fitsio.FITS(path, "rw") as fits:
        # ---------------------------------------------------------
        # PRIMARY HDU: Euclid COSEBI metadata
        # ---------------------------------------------------------
        # These values come from *any* COSEBI object because they
        # are identical across bins in official Euclid products.
        some_cosebi = next(iter(results.values()))

        hdr = {
            "THMIN": some_cosebi.thmin,
            "THMAX": some_cosebi.thmax,
            "NMODES": some_cosebi.nmodes,
            "THUNIT": "arcmin",
            "BINSCALE": "LOG",
            "COORDSYS": "SPHERICAL",
        }

        fits.write(None, header=hdr)

        # ---------------------------------------------------------
        # EXTENSIONS: one per tomographic bin pair
        # ---------------------------------------------------------
        for key, cosebi in results.items():
            type1, type2, bin1, bin2 = key
            corr_type = (type1.upper(), type2.upper())

            if corr_type != ("SHE", "SHE"):
                raise ValueError("COSEBI writer supports only SHE-SHE components.")

            # EXTNAME like SHEARSHEAR2D_COSEBI_1_1
            extname = f"SHEARSHEAR2D_COSEBI_{bin1}_{bin2}"

            # Extract arrays
            mode = np.asarray(cosebi.mode, dtype="i8")
            EE = np.asarray(cosebi.array[0], dtype="f8")
            EB = (
                np.asarray(cosebi.array[1], dtype="f8")
                if cosebi.array.shape[0] > 1
                else np.zeros_like(mode)
            )
            BB = (
                np.asarray(cosebi.array[2], dtype="f8")
                if cosebi.array.shape[0] > 2
                else np.zeros_like(mode)
            )

            # Structured dtype
            dtype = [
                ("MODE", "i8"),
                ("EE", "f8"),
                ("BB", "f8"),
                ("EB", "f8"),
            ]

            # Build table
            data = np.zeros(len(mode), dtype=dtype)
            data["MODE"] = mode
            data["EE"] = EE
            data["BB"] = BB
            data["EB"] = EB

            # Extension header
            h = {
                "BIN_ID1": int(bin1),
                "BIN_ID2": int(bin2),
            }

            history = getattr(_, "_history", None)
            if history is not None:
                h["HISTORY"] = history

            # Write the extension
            fits.write(data, extname=extname, header=h)
