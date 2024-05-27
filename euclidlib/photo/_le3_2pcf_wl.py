from __future__ import annotations

from os import PathLike
from typing import TYPE_CHECKING, Tuple, Union
from numpy.typing import NDArray
import fitsio

from euclidlib.photo import photo_data

if TYPE_CHECKING:
    from typing import Any, TypeAlias


# type alias
_DictKey: "TypeAlias" = Union[str, int, Tuple["_DictKey", ...]]


def xi_tpcf(path: str | PathLike[str]) -> dict[_DictKey, NDArray[Any]]:
    """
    Read the Weak Lensing two-point correlation functions for LE3.
    """

    corr_str = "ShearShear"

    # Load LE3 output file: 2PCF
    tpcf_le3 = photo_data.TpcfDataReader.from_fits(
        path,
        corr_str=corr_str,
    )

    z_combinations = list(tpcf_le3._data[corr_str].keys())

    # Get LE3 data
    tpcf = {}
    tpcf["THETA"] = tpcf_le3._scales["theta_uniq"]

    for i, z_comb in enumerate(z_combinations):
        tpcf[z_comb] = {}
        for component in ("+", "-"):
            tpcf[z_comb][component] = tpcf_le3.get_xi_component(component, z_comb)
        tpcf[z_comb]["WEIGHT"] = fitsio.read(path, columns="WEIGHT")

    header = fitsio.read_header(path)
    tpcf["THMIN"] = header["THMIN"]
    tpcf["THMAX"] = header["THMAX"]
    print(tpcf)
    return tpcf
