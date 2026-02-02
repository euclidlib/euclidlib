from __future__ import annotations

from os import PathLike
from typing import Optional, cast

from cosmolib.data import (
    TPCF_2Dcart,
    TPCF_2Dpol,
    TPCF_ell,
    Cov_TPCF_ell,
)

from ._common import (
    _read_le3_data,
    _read_covariance_data,
    build_covariance_matrix,
    build_2d_correlation,
    get_cosmology_from_header,
)

TYPE_CHECKING = True
if TYPE_CHECKING:
    from typing import Union

    try:
        from typing import TypeAlias
    except ImportError:
        from typing_extensions import TypeAlias

    _DictKey: TypeAlias = Union[str, int, tuple["_DictKey", ...]]


def get_TPCF_2Dcart(
    path: Union[str, PathLike[str]],
) -> dict[_DictKey, TPCF_2Dpol]:
    """
    Returns 2PCF data in the cloe-compatible euclidlib data format
    """
    header, data = _read_le3_data(path, "CORRELATION")

    z_eff, fiducial_cosmology = get_cosmology_from_header(header)

    scale_1d = data["SCALE_1D"]
    scale_2d = data["SCALE_2D"]
    xi = data["XI"]

    s_perp, s_para, correlation = build_2d_correlation(scale_1d, scale_2d, xi)

    result = TPCF_2Dcart(s_perp, s_para, correlation, fiducial_cosmology, z_eff)

    return result


def get_TPCF_2Dpol(
    path: Union[str, PathLike[str]],
) -> dict[_DictKey, TPCF_2Dpol]:
    """
    Returns 2PCF data in the cloe-compatible euclidlib data format
    """
    header, data = _read_le3_data(path, "CORRELATION")

    z_eff, fiducial_cosmology = get_cosmology_from_header(header)

    s_1d = data["SCALE_1D"]
    mu_1d = data["SCALE_2D"]
    correlation_1d = data["XI"]

    s, mu, correlation = build_2d_correlation(s_1d, mu_1d, correlation_1d)

    result = TPCF_2Dpol(s, mu, correlation, fiducial_cosmology, z_eff)

    return result


def get_TPCF_ell(
    path: Union[str, PathLike[str]],
) -> dict[_DictKey, TPCF_ell]:
    """
    Returns 2PCF data in the cloe-compatible euclidlib data format
    """
    header, data = _read_le3_data(path, "CORRELATION")

    z_eff, fiducial_cosmology = get_cosmology_from_header(header)

    results: dict[_DictKey, Optional[TPCF_ell]] = {}

    # Assuming multipoles 0, 2, 4 are present as XI0, XI2, XI4
    for i in range(5):
        for j in range(5):
            results[("SPE", "SPE", i, j)] = None

    for i in [0, 2, 4]:
        results[("SPE", "SPE", i, i)] = TPCF_ell(
            data["SCALE"],
            data[f"XI{i}"],
            fiducial_cosmology,
            z_eff,
        )

    return cast(dict[_DictKey, TPCF_ell], results)


def get_Cov_TPCF_ell(
    path: Union[str, PathLike[str]],
) -> Cov_TPCF_ell:
    """
    Returns a single Cov_TPCF_ell object containing the full,
    combined covariance matrix for even multipoles (0, 2, 4), the s-axis,
    and the effective redshift.
    """
    header, data = _read_covariance_data(path)

    z_eff = get_cosmology_from_header(header, get_fiducial=False)

    s_values, full_cov_matrix = build_covariance_matrix(data, "TPCF")

    # Create and return the single covariance object
    result = Cov_TPCF_ell(
        s_values,
        full_cov_matrix,
        z_eff,
    )

    return result
