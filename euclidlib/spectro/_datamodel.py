from __future__ import annotations
from dataclass import dataclass
from numpy.typing import NDArray
from typing import Any, Dict

@dataclass(frozen=True)
class Result:
    """
    Output datamodel for read LE3 GC PK measurements, containing information on the header and the power spectrum multipoles.

    Attributes
    ----------
    k : ndarray
        centers of k bins
    k_eff : ndarray
        effective values of k bins
    mode_number : ndarray
        number of modes in each k bin
    p : dict[ndarray]
        multipoles of the power spectrum (keys are integer l's from 0 to 4)
    header : dict
        fits header of the original measurement file
    """
    k: NDArray[Any]
    k_eff: NDArray[Any]
    mode_number: NDArray[Any]
    p: Dict[int, NDArray[Any]]
    header: Dict[str, Any]

    def __getitem__(self, l: int) -> NDArray[Any]:
        return self.p[l]

    @property
    def shot_noise(self) -> float:
        return self.header["SN_VALUE"]
