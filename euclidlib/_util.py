"""
Module for internal utility functions.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.lib import NumpyVersion

if TYPE_CHECKING:
    from typing import Any, Callable, TypeVar

    AnyT = TypeVar("AnyT", bound=Any)


def writer(func: Any) -> Callable[[AnyT], AnyT]:
    """
    Decorator for writer functions.
    """

    def decorator(writer: AnyT) -> AnyT:
        func.write = writer
        writer.__module__ = func.__module__
        writer.__name__ = "write"
        writer.__qualname__ = f"{func.__qualname__}.write"
        return writer

    return decorator

def trapezoidal_integration(y, x, axis=-1):
    """
    Compute the integral of `y` values using the trapezoidal rule.

    This function uses `np.trapezoid` if available (NumPy version 2.0.0 or later),
    otherwise it falls back to `np.trapz` for older versions of NumPy.

    This function will be eventually deprecated as soon as NumPy version 2.0.0 
    becomes popular enough.

    Parameters
    ----------
    y : array_like
        Input array to integrate.
    x : array_like, optional
        The sample points corresponding to the `y` values. If `x` is None, the
        sample points are assumed to be evenly spaced.
    axis : int, optional
        The axis along which to integrate. Default is -1 (last axis).

    Returns
    -------
    float or ndarray
        Definite integral as approximated by the trapezoidal rule.
    """
    if NumpyVersion(np.__version__) >= NumpyVersion('2.0.0'):
        # NumPy version is 2.0.0' or later
        return np.trapezoid(y, x, axis=axis)
    else:
        # NumPy version is older than 2.0.0'
        return np.trapz(y, x, axis=axis)
