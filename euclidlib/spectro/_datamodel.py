from __future__ import annotations
from dataclasses import dataclass
from numpy.typing import NDArray
from typing import Any, Dict

@dataclass(frozen=True)
class Result():
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

    def __post_init__(self) -> None:
        # Sanity check on the attributes
        if any(
            [
                len(self.k) != len(self.k_eff),
                any([len(self.k) != len(p_el) for p_el in self.p.values()]),
                len(self.k) != len(self.mode_number)
            ]
        ):
            raise ValueError("Inconsistent class attributes, all arrays must have the same length.")
        for l in range(5):
            if l not in self.p:
                raise ValueError("Power spectrum attribute must contain all multipoles from 0 to 4.")

    def __getitem__(self, l: int) -> NDArray[Any]:
        return self.p[l]
        
    def __len__(self) -> int:
        return len(self.k)

    @property
    def shot_noise(self) -> float:
        return self.header["SN_VALUE"]
