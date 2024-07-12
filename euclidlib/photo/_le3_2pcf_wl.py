from __future__ import annotations

from os import PathLike
import fitsio
from euclidlib.photo import photo_data


def correlation_functions(path: str | PathLike[str]) -> dict:
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

    return tpcf
