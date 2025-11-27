from __future__ import annotations
import re
import os
import fitsio  # type: ignore [import-not-found]
import numpy as np
from .._util import writer

from cosmolib.data import AngularPowerSpectrum  # type: ignore [import-not-found]

TYPE_CHECKING = True
if TYPE_CHECKING:
    from os import PathLike
    from typing import Any, Union

    try:
        from typing import TypeAlias
    except ImportError:
        from typing_extensions import TypeAlias
    from numpy.typing import NDArray

    _DictKey: TypeAlias = Union[str, int, tuple["_DictKey", ...]]


def normalize_result_axis(
    axis: tuple[int, ...] | int | None,
    result: NDArray[Any],
    ell: tuple[NDArray[Any], ...] | NDArray[Any] | None,
) -> tuple[int, ...]:
    """
    Normalize the axis used to store results.

    Parameters
    ----------
    axis : tuple of int, int, or None
        Axis or axes over which the result is computed.
    result : NDArray
        The result array.
    ell : tuple of NDArray, NDArray, or None
        Associated ell values.

    Returns
    -------
    axis : tuple of int
        Normalized tuple of axes.
    """
    try:
        from numpy.lib.array_utils import normalize_axis_tuple
    except ModuleNotFoundError:
        from numpy.lib.stride_tricks import normalize_axis_tuple  # type: ignore

    if axis is None:
        if result.ndim == 0:
            axis = ()
        elif isinstance(ell, tuple):
            axis = tuple(range(-len(ell), 0))
        else:
            axis = -1
    return normalize_axis_tuple(axis, result.ndim, "axis")


def _key_from_string(s: str) -> _DictKey:
    """
    Decode a string key from a FITS extension name.

    Parameters
    ----------
    s : str
        Encoded key as a string.

    Returns
    -------
    key : str, int, or tuple
        Decoded dictionary key.
    """
    parts = re.split(r"(?<!\\)-", s.replace("\\\\", "\0"))
    if len(parts) > 1:
        return tuple(map(_key_from_string, parts))
    key = parts[0]
    key = key.replace("\\-", "-")
    key = key.replace("\0", "\\")
    return int(key) if key.removeprefix("-").isdigit() else key


def _read_metadata(hdu: Any) -> dict[str, Any]:
    """
    Read array metadata from a FITS HDU header.

    Parameters
    ----------
    hdu : fitsio.FITS_rec
        FITS HDU object.

    Returns
    -------
    metadata : dict
        Metadata dictionary from header entries prefixed with 'META '.
    """
    h = hdu.read_header()
    md = {}
    for key in h:
        if key.startswith("META "):
            md[key[5:].lower()] = h[key]
    return md


def _read_result(hdu: Any) -> AngularPowerSpectrum:
    """
    Read a result object from a FITS HDU.

    Parameters
    ----------
    hdu : fitsio.FITS_rec
        FITS HDU object.

    Returns
    -------
    result : Result
        Parsed result object.
    """
    from ast import literal_eval

    data = hdu.read()
    h = hdu.read_header()

    axis = literal_eval(h["ELLAXIS"])
    arr = np.moveaxis(data["ARRAY"], tuple(range(len(axis))), axis)
    order = np.argsort(axis)

    _ell = data["ELL"]
    ell = (
        _ell if _ell.ndim == 1 else tuple(_ell[: arr.shape[axis[i]], i] for i in order)
    )

    _lower = data["LOWER"]
    lower = (
        _lower
        if _lower.ndim == 1
        else tuple(_lower[: arr.shape[axis[i]], i] for i in order)
    )

    _upper = data["UPPER"]
    upper = (
        _upper
        if _upper.ndim == 1
        else tuple(_upper[: arr.shape[axis[i]], i] for i in order)
    )

    _weight = data["WEIGHT"]
    weight = (
        _weight
        if _weight.ndim == 1
        else tuple(_weight[: arr.shape[axis[i]], i] for i in order)
    )

    return AngularPowerSpectrum(
        arr.view(np.dtype(arr.dtype, metadata=_read_metadata(hdu))),
        axis=tuple(axis[i] for i in order),
        ell=ell,
        lower=lower,
        upper=upper,
        weight=weight,
    )


def read(path: str | PathLike[str]) -> dict[_DictKey, AngularPowerSpectrum]:
    """
    Read a dictionary of results from a FITS file.

    Parameters
    ----------
    path : str or PathLike
        Path to the FITS file.

    Returns
    -------
    results : dict
        Dictionary mapping decoded keys to dataclass objects.
    """
    results = {}
    with fitsio.FITS(path) as fits:
        for hdu in fits:
            if not hdu.has_data():
                continue
            ext = hdu.get_extname()
            if not ext:
                continue
            key = _key_from_string(ext)
            if not key:
                continue
            results[key] = _read_result(hdu)
    return results


def angular_power_spectra(
    path: str | PathLike[str],
) -> dict[_DictKey, AngularPowerSpectrum]:
    """
    Read angular power spectra results from a FITS file.

    Parameters
    ----------
    path : str or PathLike
        Path to the FITS file.

    Returns
    -------
    spectra : dict
        Dictionary of angular power spectra dataclass objects.
    """
    return read(path)


def mixing_matrices(path: str | PathLike[str]) -> dict[_DictKey, AngularPowerSpectrum]:
    """
    Read mixing matrices from a FITS file.

    Parameters
    ----------
    path : str or PathLike
        Path to the FITS file.

    Returns
    -------
    matrices : dict
        Dictionary of mixing matrix dataclass objects.
    """
    return read(path)


@writer(angular_power_spectra)
def _(path: str | PathLike[str], results: dict[_DictKey, AngularPowerSpectrum]) -> None:
    """
    Write angular power spectra results to a FITS file.
    Parameters
    ----------
    path : str or PathLike
        Path to the output FITS file.
    results : dict
        Dictionary mapping keys to dataclass objects containing the data to write.
    """

    def _key_to_string(key: _DictKey) -> str:
        return "-".join(str(k) for k in (key if isinstance(key, tuple) else (key,)))

    if os.path.exists(path):
        os.remove(path)

    with fitsio.FITS(path, "rw") as fits:
        for key, result in results.items():
            name = _key_to_string(key)
            arr = np.asarray(result.array)
            axis = result.axis
            if isinstance(axis, int):
                axis = (axis,)
            elif axis is None:
                axis = ()
            if arr.ndim == 1:
                nrows = arr.shape[0]
                reshaped_arr = arr
                tdim = None
            elif arr.ndim in (2, 3):
                dim_shape = arr.shape[:-1]
                nrows = arr.shape[-1]
                reshaped_arr = arr.reshape(-1, nrows).T
                tdim = f"({','.join(map(str, dim_shape[::-1]))})"
            else:
                raise ValueError(f"Unsupported array shape: {arr.shape}")

            def get_tuple_or_default(
                attr: str, default_dtype: np.dtype[Any]
            ) -> np.ndarray[Any, Any]:
                val = getattr(result, attr, None)
                if val is None:
                    return np.zeros(nrows, dtype=default_dtype)
                val = np.asarray(val)
                if isinstance(val, (tuple, list)):
                    val = np.stack([np.asarray(v) for v in val], axis=-1)
                return val.reshape(nrows).astype(default_dtype)

            ell = get_tuple_or_default("ell", np.dtype(np.int64))
            lower = get_tuple_or_default("lower", np.dtype(np.int64))
            upper = get_tuple_or_default("upper", np.dtype(np.int64))
            weight = get_tuple_or_default("weight", np.dtype(np.float64))

            array_shape = (
                tuple(map(int, tdim.strip("()").split(",")))
                if tdim
                else reshaped_arr.shape[1:]
            )

            dtype = [
                ("ARRAY", "f8", array_shape),
                ("ELL", "i8"),
                ("LOWER", "i8"),
                ("UPPER", "i8"),
                ("WEIGHT", "f8"),
            ]

            structured_array = np.empty(nrows, dtype=dtype)
            structured_array["ARRAY"] = reshaped_arr.reshape((nrows,) + array_shape)
            structured_array["ELL"] = ell
            structured_array["LOWER"] = lower
            structured_array["UPPER"] = upper
            structured_array["WEIGHT"] = weight

            header = {
                "ELLAXIS": repr(axis),
            }
            if tdim is not None:
                header["TDIM1"] = tdim
            if meta := arr.dtype.metadata:
                for k, v in meta.items():
                    header[f"META {k.upper()}"] = str(v)

            # Write the header to the FITS file
            history = getattr(_, "_history", None)
            if history is not None:
                header["HISTORY"] = history
            fits.write(structured_array, extname=name, header=header)
