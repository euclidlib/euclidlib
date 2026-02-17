from __future__ import annotations

import numpy as np
from os import PathLike

from typing import Optional, Any
from numpy.typing import NDArray

from cosmolib.data import (
    TwoPointCorrelationCartesian,
    TwoPointCorrelationPolar,
    TwoPointCorrelationMultipoles,
    TwoPointCorrelationMultipolesCovariance,
)

from ._common import (
    read_data_vectors,
    read_covariance_data,
    build_2d_correlation,
    get_cosmology_from_header,
    check_input,
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
    even_multipoles = [0, 2, 4]

    redshifts, nz = check_input(redshifts)
    result: dict[_DictKey, Optional[TwoPointCorrelationMultipoles]] = {}

    for i in range(nz):
        for j in range(nz):
            result[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        header, data = read_covariance_data(str(path).format(zlab))
        zeff, _ = get_cosmology_from_header(header, get_fiducial=False)
        mask = np.isin(data["MULTIPOLE-I"], even_multipoles) & np.isin(
            data["MULTIPOLE-J"], even_multipoles
        )
        data = data[mask]

        scale_prefix = "s"
        scale_i_col = f"{scale_prefix.upper()}I"

        k_values = np.unique(data[scale_i_col])
        nk = len(k_values)

        covariance_blocks: dict[str, NDArray[Any]] = {}
        for ell_i in even_multipoles:
            for ell_j in even_multipoles:
                block_mask = (data["MULTIPOLE-I"] == ell_i) & (
                    data["MULTIPOLE-J"] == ell_j
                )
                block_data = data[block_mask]

                if len(block_data) == 0:
                    continue

                current_block = np.zeros((nk, nk))

                k_to_idx = {k_val: i for i, k_val in enumerate(k_values)}

                for row in block_data:
                    k_idx_i = k_to_idx[row[f"{scale_prefix.upper()}I"]]
                    k_idx_j = k_to_idx[row[f"{scale_prefix.upper()}J"]]
                    current_block[k_idx_i, k_idx_j] = row["COVARIANCE"]

                covariance_blocks[f"ELL_{ell_i}-{ell_j}"] = current_block

        result[("SPE", "SPE", i, i)] = TwoPointCorrelationMultipolesCovariance(
            k=k_values,
            covariance=covariance_blocks,
            zeff=zeff,
        )

    return result
