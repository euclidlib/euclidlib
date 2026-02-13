from __future__ import annotations

from warnings import warn
from os import PathLike

import numpy as np
import fitsio  # type: ignore [import-not-found]
from cosmolib.data import (
    PowerSpectrumMultipoles,
    PowerSpectrumMultipolesCovariance,
    PowerSpectrumMultipolesMixingMatrix,
)

from ._common import (
    _read_le3_data,
    _read_covariance_data,
    build_covariance_matrix,
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


def get_PS_ell(
    path: Union[str, PathLike[str]],
) -> dict[_DictKey, PowerSpectrumMultipoles]:
    """
    Returns power spectrum data in the cloe-compatible euclidlib data format
    """
    header, data = _read_le3_data(path, "SPECTRUM", check_extra_hdu=True)

    z_eff, fiducial_cosmology = get_cosmology_from_header(header)

    multipoles: dict[int, np.ndarray] = {}

    for ell in range(5):
        multipoles[f"ELL{ell}"] = data[f"PK{ell}"]

    result = PowerSpectrumMultipoles(
        data["K"],
        data["K_EFF"],
        data["NUM_MOD"],
        multipoles,
        fiducial_cosmology,
        z_eff,
        1.0 / header["SN_VALUE"],
        header["SN_VALUE"],
    )

    return result


def get_Cov_PS_ell(
    path: Union[str, PathLike[str]],
) -> PowerSpectrumMultipolesCovariance:
    """
    Returns a single Cov_PS_ell object containing the full,
    combined covariance matrix for even multipoles (0, 2, 4), the k-axis,
    and the effective redshift.
    """
    header, data = _read_covariance_data(path)

    z_eff = get_cosmology_from_header(header, get_fiducial=False)

    k_values, full_cov_matrix = build_covariance_matrix(data, "PS")

    # Create and return the single PowerSpectrumMultipolesCovariance object
    result = PowerSpectrumMultipolesCovariance(
        k_values,
        full_cov_matrix,
        z_eff,
    )

    return result


def get_MixMat_PS_ell(
    path: Union[str, PathLike[str]],
) -> PowerSpectrumMultipolesMixingMatrix:
    """
    Reads the mixing matrix for power spectrum multipoles from a FITS file.
    """
    with fitsio.FITS(path) as fits:
        kout_data = fits["BINS_OUTPUT"].read()
        kout = kout_data["k"]

        kin_data = fits["BINS_INPUT"].read()
        kin = {
            0: kin_data["kp0"],
            2: kin_data["kp2"],
            4: kin_data["kp4"],
        }

        mixing_data = fits["MIXING_MATRIX"].read()
        header = fits["MIXING_MATRIX"].read_header()

        W00 = mixing_data["W00"]
        W02 = mixing_data["W02"]
        W04 = mixing_data["W04"]
        W20 = mixing_data["W20"]
        W22 = mixing_data["W22"]
        W24 = mixing_data["W24"]
        W40 = mixing_data["W40"]
        W42 = mixing_data["W42"]
        W44 = mixing_data["W44"]

        mixing_matrix = np.block([[W00, W02, W04], [W20, W22, W24], [W40, W42, W44]])

        try:
            z_eff = header["Z_EFF"]
        except KeyError:
            warn(
                "Effective redshift not specified in FITS file header. Setting it to 0."
            )
            z_eff = 0.0

    result = PowerSpectrumMultipolesMixingMatrix(
        kout,
        kin,
        mixing_matrix,
        z_eff,
    )

    return result
