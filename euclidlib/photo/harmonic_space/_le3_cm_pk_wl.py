from __future__ import annotations

import re

import fitsio  # type: ignore [import-not-found]
import numpy as np

from dataclasses import dataclass

TYPE_CHECKING = False
if TYPE_CHECKING:
    from os import PathLike
    from typing import Any, TypeAlias
    from numpy.typing import NDArray

    _DictKey: TypeAlias = str | int | tuple["_DictKey", ...]