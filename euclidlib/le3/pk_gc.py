from __future__ import annotations

from warnings import warn
from os import PathLike
import os
import fitsio  # type: ignore [import-not-found]

from typing import Optional, cast

from numpy.typing import NDArray
from cosmolib.data import PowerSpectrum  # type: ignore [import-not-found]

TYPE_CHECKING = True
if TYPE_CHECKING:
    from typing import Any, Dict, Tuple, Union

    try:
        from typing import TypeAlias
    except ImportError:
        from typing_extensions import TypeAlias

    _DictKey: TypeAlias = Union[str, int, tuple["_DictKey", ...]]


def _verify_input_file(
    path: Union[str, PathLike[str]], check_extra_hdu: bool = False
) -> None:
    """
    Verifies that input file is a compatible LE3-GC fits file
    """
    hdul_target_length = 3 if check_extra_hdu else 2

    if not any((isinstance(path, str), isinstance(path, PathLike))):
        raise TypeError(
            "Provided fits file name must be a string or a PathLike object."
        )

    if not os.path.isfile(path):
        raise FileNotFoundError(f"Could not find file '{path}'.")

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
    """
    head: Dict[str, Any] = dict(hdu.read_header())
    return head


def _get_hdu_data(hdu: fitsio.TableHDU) -> NDArray[Any]:
    """
    Reads data from fits file HDU
    """
    arr: NDArray[Any] = hdu.read()
    return arr


def _read_power_spectrum(
    path: Union[str, PathLike[str]],
) -> Tuple[Dict[str, Any], NDArray[Any]]:
    """
    Reads power spectrum from Euclid LE3-GC fits file
    """
    pk_not_found_message = (
        "Invalid fits file provided, cannot locate power spectrum data."
    )

    _verify_input_file(path)

    with fitsio.FITS(path) as fits_input:
        header = _get_hdu_header(fits_input[1])

        if "SPECTRUM" not in header["EXTNAME"]:
            raise ValueError(pk_not_found_message)

        if header["EXTNAME"] == "SPECTRUM":
            data = _get_hdu_data(fits_input[1])

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
    path: Union[str, PathLike[str]],
) -> dict[_DictKey, PowerSpectrum]:
    """
    Returns power spectrum data in the cloe-compatible euclidlib data format
    """
    header, data = _read_power_spectrum(path)

    fidual_cosmology = {
        "OMEGA_M": header["OMEGA_M"],
        "OMEGA_R": header["OMEGA_R"],
        "OMEGA_B": header["OMEGA_B"],
        "OMEGA_V": header["OMEGA_V"],
        "OMEGA_K": header["OMEGA_K"],
        "HUBBLE": header["HUBBLE"],
        "INDEX_N": header["INDEX_N"],
        "SIGMA_8": header["SIGMA_8"],
        "W_STATE": header["W_STATE"],
        "N_EFF": header["N_EFF"],
        "T_CMB": header["T_CMB"],
    }

    # Explicit annotation required for mypy
    results: dict[_DictKey, Optional[PowerSpectrum]] = {}

    # Create all keys (auto + cross)
    for i in range(5):
        for j in range(5):
            results[("SPE", "SPE", i, j)] = None

    # Fill only auto-correlations
    for i in range(5):
        results[("SPE", "SPE", i, i)] = PowerSpectrum(
            data["K"],
            data["K_EFF"],
            data["NUM_MOD"],
            data[f"PK{i}"],
            fidual_cosmology,
            0,  # redshift_eff is not stored
            1.0 / header["SN_VALUE"],
            header["SN_VALUE"],
        )

    return cast(dict[_DictKey, PowerSpectrum], results)
