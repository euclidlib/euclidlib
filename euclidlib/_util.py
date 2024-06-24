"""
Module for internal utility functions.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

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
