from __future__ import annotations

import os
import numpy as np
import fitsio # type: ignore[import-untyped]
from os import PathLike

from typing import Optional
from numpy.typing import NDArray

from cosmolib.data import (
    TwoPointCorrelationCartesian,
    TwoPointCorrelationPolar,
    TwoPointCorrelationMultipoles,
    TwoPointCorrelationMultipolesCovariance,
)

from ._common import (
    check_input,
    get_cosmology_from_header,
    read_data_vectors,
    read_and_reshape_covariance_matrix,
    build_2d_correlation,
)

TYPE_CHECKING = True
if TYPE_CHECKING:
    from typing import Union

    try:
        from typing import TypeAlias
    except ImportError:
        from typing_extensions import TypeAlias

    _DictKey: TypeAlias = Union[str, int, tuple["_DictKey", ...]]


def get_TwoPointCorrelationCartesian(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, TwoPointCorrelationCartesian]:
    """
    Returns 2PCF data in the cloe-compatible euclidlib data format
    """
    redshifts, nz = check_input(redshifts)
    result: dict[_DictKey, Optional[TwoPointCorrelationCartesian]] = {}

    for i in range(nz):
        for j in range(nz):
            result[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        header, data = read_data_vectors(str(path).format(zlab), "CORRELATION")
        zeff, fiducial_cosmology = get_cosmology_from_header(header)

        scale_1d = data["SCALE_1D"]
        scale_2d = data["SCALE_2D"]
        xi = data["XI"]

        s_perp, s_para, correlation = build_2d_correlation(scale_1d, scale_2d, xi)

        result[("SPE", "SPE", i, i)] = TwoPointCorrelationCartesian(
            s_perp=s_perp,
            s_para=s_para,
            correlation=correlation,
            fiducial_cosmology=fiducial_cosmology,
            zeff=zeff,
        )

    return result


def get_TwoPointCorrelationPolar(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, TwoPointCorrelationPolar]:
    """
    Returns 2PCF data in the cloe-compatible euclidlib data format
    """
    redshifts, nz = check_input(redshifts)
    result: dict[_DictKey, Optional[TwoPointCorrelationPolar]] = {}

    for i in range(nz):
        for j in range(nz):
            result[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        header, data = read_data_vectors(str(path).format(zlab), "CORRELATION")
        zeff, fiducial_cosmology = get_cosmology_from_header(header)

        s_1d = data["SCALE_1D"]
        mu_1d = data["SCALE_2D"]
        correlation_1d = data["XI"]

        s, mu, correlation = build_2d_correlation(s_1d, mu_1d, correlation_1d)

        result[("SPE", "SPE", i, i)] = TwoPointCorrelationPolar(
            s=s,
            mu=mu,
            correlation=correlation,
            fiducial_cosmology=fiducial_cosmology,
            zeff=zeff,
        )

    return result


def get_TwoPointCorrelationMultipoles(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, TwoPointCorrelationMultipoles]:
    """
    Returns 2PCF data in the cloe-compatible euclidlib data format
    """
    redshifts, nz = check_input(redshifts)
    result: dict[_DictKey, Optional[TwoPointCorrelationMultipoles]] = {}

    for i in range(nz):
        for j in range(nz):
            result[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        header, data = read_data_vectors(str(path).format(zlab), "CORRELATION")
        zeff, fiducial_cosmology = get_cosmology_from_header(header)
        multipoles = np.array([data[f"XI{ell}"] for ell in range(5)])

        result[("SPE", "SPE", i, i)] = TwoPointCorrelationMultipoles(
            s=data["SCALE"],
            multipoles=multipoles,
            fiducial_cosmology=fiducial_cosmology,
            zeff=zeff,
        )

    return result


def get_TwoPointCorrelationMultipolesCovariance(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, TwoPointCorrelationMultipolesCovariance]:
    """
    Returns a single Cov_TPCF_ell object containing the full,
    combined covariance matrix for even multipoles (0, 2, 4), the s-axis,
    and the effective redshift.
    """
    redshifts, nz = check_input(redshifts)
    result: dict[_DictKey, Optional[TwoPointCorrelationMultipolesCovariance]] = {}

    for i in range(nz):
        for j in range(nz):
            result[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        s_values, covariance_blocks, zeff = read_and_reshape_covariance_matrix(
            path=str(path).format(zlab),
            type="CORRELATION",
        )

        result[("SPE", "SPE", i, i)] = TwoPointCorrelationMultipolesCovariance(
            s=s_values,
            covariance=covariance_blocks,
            zeff=zeff,
        )

    return result


def write_TwoPointCorrelationMultipoles(
    path: Union[str, PathLike[str]],
    s: NDArray[np.float64],
    xi: NDArray[np.float64],
    zeff: float,
    cosmology: dict[str, float]
) -> None:
    """
    Writes two-point correlation function multipoles to a FITS file with a structure
    consistent with Euclid LE3 GC products.
    """
    nrows = len(s)
    dtype = [
        ("SCALE", "f8"),
        ("XI0", "f8"),
        ("XI1", "f8"),
        ("XI2", "f8"),
        ("XI3", "f8"),
        ("XI4", "f8"),
    ]
    data = np.zeros(nrows, dtype=dtype)
    data["SCALE"] = s
    for ell in range(min(5, xi.shape[0])):
        data[f"XI{ell}"] = xi[ell]

    header = {
        "TELESCOP": "EUCLID  ",
        "INSTRUME": "LE3GC   ",
        "RUNTYPE": "AUTO    ",
        "Z_EFF": zeff,
        "STAT": "MULTIPOLE",
        "BIN1TYPE": "LIN     ",
        "BIN1NUM": nrows,
        "BIN1MIN": np.min(s) if nrows > 0 else 0.0,
        "BIN1MAX": np.max(s) if nrows > 0 else 0.0,
        "TUNIT1": "Mpc/h   ",
        "TUNIT2": "        ",
        "TUNIT3": "        ",
        "TUNIT4": "        ",
        "TUNIT5": "        ",
        "TUNIT6": "        ",
    }

    header.update(cosmology)

    if os.path.exists(path):
        os.remove(path)

    with fitsio.FITS(path, "rw") as f:
        f.write(None, header={"EXTNAME": "PRIMARY"})
        f.write(data, header=header, extname="CORRELATION")
