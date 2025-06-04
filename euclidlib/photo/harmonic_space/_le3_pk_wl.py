from __future__ import annotations

import re
import os
import fitsio  # type: ignore [import-not-found]
import numpy as np
from datetime import datetime
from ..._util import writer

from dataclasses import dataclass

TYPE_CHECKING = False
if TYPE_CHECKING:
    from os import PathLike
    from typing import Any, TypeAlias
    from numpy.typing import NDArray

    _DictKey: TypeAlias = str | int | tuple["_DictKey", ...]


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


@dataclass(frozen=True, repr=False)
class Result:
    """
    A container for an array and associated data such as errors and ell values.

    Parameters
    ----------
    array : NDArray
        Main data array.
    ell : NDArray or tuple of NDArray, optional
        Angular scale(s) or ell values associated with the data.
    axis : int or tuple of int, optional
        Axis/axes corresponding to ell in the array.
    lower : NDArray or tuple of NDArray, optional
        Lower error bounds.
    upper : NDArray or tuple of NDArray, optional
        Upper error bounds.
    weight : NDArray or tuple of NDArray, optional
        Weights or inverse variances.
    """

    array: NDArray[Any]
    ell: NDArray[Any] | tuple[NDArray[Any], ...] | None = None
    axis: int | tuple[int, ...] | None = None
    lower: NDArray[Any] | tuple[NDArray[Any], ...] | None = None
    upper: NDArray[Any] | tuple[NDArray[Any], ...] | None = None
    weight: NDArray[Any] | tuple[NDArray[Any], ...] | None = None

    def __post_init__(self) -> None:
        float_array = np.asarray(self.array, dtype=float)
        object.__setattr__(self, "array", float_array)
        axis = normalize_result_axis(self.axis, self.array, self.ell)
        object.__setattr__(self, "axis", axis)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(axis={self.axis!r})"

    def __array__(
        self,
        dtype: np.dtype[Any] | None = None,
        *,
        copy: np.bool[bool] | None = None,
    ) -> NDArray[Any]:
        """Return NumPy array representation of the result."""
        if copy is not None:
            return self.array.__array__(dtype, copy=copy)
        return self.array.__array__(dtype)

    def __getitem__(self, key: Any) -> Any:
        """Return a slice or element of the internal array."""
        return self.array[key]

    @property
    def ndim(self) -> int:
        """Return the number of dimensions of the result array."""
        return self.array.ndim

    @property
    def shape(self) -> tuple[int, ...]:
        """Return the shape of the result array."""
        return self.array.shape

    @property
    def dtype(self) -> np.dtype[Any]:
        """Return the data type of the result array."""
        return self.array.dtype


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


def _key_to_string(key: _DictKey) -> str:
    if isinstance(key, tuple):
        return "-".join(_key_to_string(k).replace("-", "\\-") for k in key)
    if isinstance(key, int):
        return str(key)
    return key.replace("-", "\\-").replace("\\", "\\\\")


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


def _read_result(hdu: Any) -> Result:
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

    return Result(
        arr.view(np.dtype(arr.dtype, metadata=_read_metadata(hdu))),
        axis=tuple(axis[i] for i in order),
        ell=ell,
        lower=lower,
        upper=upper,
        weight=weight,
    )


def read(path: str | PathLike[str]) -> dict[_DictKey, Result]:
    """
    Read a dictionary of results from a FITS file.

    Parameters
    ----------
    path : str or PathLike
        Path to the FITS file.

    Returns
    -------
    results : dict
        Dictionary mapping decoded keys to Result objects.
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


def angular_power_spectra(path: str | PathLike[str]) -> dict[_DictKey, Result]:
    """
    Read angular power spectra results from a FITS file.

    Parameters
    ----------
    path : str or PathLike
        Path to the FITS file.

    Returns
    -------
    spectra : dict
        Dictionary of angular power spectra Result objects.
    """
    return read(path)


def mixing_matrices(path: str | PathLike[str]) -> dict[_DictKey, Result]:
    """
    Read mixing matrices from a FITS file.

    Parameters
    ----------
    path : str or PathLike
        Path to the FITS file.

    Returns
    -------
    matrices : dict
        Dictionary of mixing matrix Result objects.
    """
    return read(path)


def covariance_matrices(path: str | PathLike[str]) -> dict[_DictKey, Result]:
    """
    Read covariance matrices from a FITS file.

    Parameters
    ----------
    path : str or PathLike
        Path to the FITS file.

    Returns
    -------
    matrices : dict
        Dictionary of covariance matrix Result objects.
    """
    return read(path)


@writer(angular_power_spectra)
def _(path: str | PathLike[str], results: dict[_DictKey, Result]) -> None:
    """
    Write angular power spectra results to a FITS file.
    Parameters
    ----------
    path : str or PathLike
        Path to the output FITS file.
    results : dict
        Dictionary mapping keys to Result objects containing the data to write.
    """

    def _key_to_string(key):
        return "-".join(map(str, key))

    if os.path.exists(path):
        os.remove(path)

    timestamp = datetime.now().isoformat(timespec="seconds")

    with fitsio.FITS(path, "rw") as fits:
        for key, result in results.items():
            name = _key_to_string(key)
            arr = np.asarray(result.array)
            axis = result.axis
            if not isinstance(axis, tuple):
                axis = (axis,)
            order = np.argsort(axis)

            # Determine the shape of one row and reshape arr
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

            def get_tuple_or_default(attr, default_dtype):
                val = getattr(result, attr, None)
                if val is None:
                    return np.zeros(nrows, dtype=default_dtype)
                val = np.asarray(val)
                if val.ndim > 1:
                    val = np.stack(val, axis=-1)
                return val.reshape(nrows).astype(default_dtype)

            ell = get_tuple_or_default("ell", np.int64)
            lower = get_tuple_or_default("lower", np.int64)
            upper = get_tuple_or_default("upper", np.int64)
            weight = get_tuple_or_default("weight", np.float64)

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
            if tdim:
                structured_array["ARRAY"] = reshaped_arr.reshape((nrows,) + array_shape)
            else:
                structured_array["ARRAY"] = reshaped_arr

            structured_array["ELL"] = ell
            structured_array["LOWER"] = lower
            structured_array["UPPER"] = upper
            structured_array["WEIGHT"] = weight

            # Build header
            header = {
                "ELLAXIS": repr(axis),
                "HISTORY": f"Generated by euclidlib on {timestamp}",
            }
            if tdim is not None:
                header["TDIM1"] = tdim
            if meta := arr.dtype.metadata:
                for k, v in meta.items():
                    header[f"META {k.upper()}"] = str(v)

            fits.write(structured_array, extname=name, header=header)
