"""
Module for internal utility functions.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.lib import NumpyVersion
from typing import Union, Optional

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

def trapezoidal_integration(y: np.ndarray, x: Optional[np.ndarray] = None, axis: int = -1) -> Union[float, np.ndarray]:
    """
    Compute the integral of `y` values using the trapezoidal rule.
    
    Parameters
    ----------
    y : np.ndarray
        Input array to integrate.
    x : Optional[np.ndarray], optional
        The sample points corresponding to the `y` values. If `x` is None, the
        sample points are assumed to be evenly spaced.
    axis : int, optional
        The axis along which to integrate. Default is -1 (last axis).
    
    Returns
    -------
    float or np.ndarray
        Definite integral as approximated by the trapezoidal rule.
    """
    if NumpyVersion(np.__version__) >= NumpyVersion('2.0.0'):
        # NumPy version is 2.0.0 or later
        if hasattr(np, 'trapezoid'):
            return np.trapezoid(y, x, axis=axis)  # Ensure that `np.trapezoid` exists
        else:
            raise AttributeError("NumPy version is >= 2.0.0 but 'trapezoid' function is not available.")
    else:
        # NumPy version is older than 2.0.0
        if hasattr(np, 'trapz'):
            return np.trapz(y, x, axis=axis)
        else:
            raise AttributeError("NumPy version is < 2.0.0 and 'trapz' function is not available.")
        
