from __future__ import annotations

from warnings import warn
from os import PathLike

from typing import Optional

import numpy as np
from cosmolib.data import (
    PowerSpectrumMultipoles,
    PowerSpectrumMultipolesCovariance,
    PowerSpectrumMultipolesMixingMatrix,
)

from ._common import (
    check_input,
    get_cosmology_from_header,
    read_data_vectors,
    read_and_reshape_covariance_matrix,
    read_mixing_matrix,
)

TYPE_CHECKING = True
if TYPE_CHECKING:
    from typing import Union

    try:
        from typing import TypeAlias
    except ImportError:
        from typing_extensions import TypeAlias

    _DictKey: TypeAlias = Union[str, int, tuple["_DictKey", ...]]


def get_PowerSpectrumMultipoles(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, PowerSpectrumMultipoles]:
    """
    Returns power spectrum data in the cloe-compatible euclidlib data format
    """
    redshifts, nz = check_input(redshifts)
    result: dict[_DictKey, Optional[PowerSpectrumMultipoles]] = {}

    for i in range(nz):
        for j in range(nz):
            result[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        header, data = read_data_vectors(str(path).format(zlab), "SPECTRUM")
        zeff, fiducial_cosmology = get_cosmology_from_header(header)
        nbar = 1.0 / header["SN_VALUE"]
        multipoles = np.array([data[f"PK{ell}"] for ell in range(5)])

        result[("SPE", "SPE", i, i)] = PowerSpectrumMultipoles(
            k=data["K"],
            keff=data["K_EFF"],
            Nmodes=data["NUM_MOD"],
            multipoles=multipoles,
            fiducial_cosmology=fiducial_cosmology,
            zeff=zeff,
            nbar=nbar,
            Psn=header["SN_VALUE"],
        )

    return result


def get_PowerSpectrumMultipolesCovariance(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, PowerSpectrumMultipolesCovariance]:
    """
    Returns a single Cov_PS_ell object containing the full,
    combined covariance matrix for even multipoles (0, 2, 4), the k-axis,
    and the effective redshift.
    """
    redshifts, nz = check_input(redshifts)
    result: dict[_DictKey, Optional[PowerSpectrumMultipolesCovariance]] = {}

    for i in range(nz):
        for j in range(nz):
            result[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        k_values, covariance_blocks, zeff = read_and_reshape_covariance_matrix(
            path=str(path).format(zlab),
            type="SPECTRUM",
        )

        result[("SPE", "SPE", i, i)] = PowerSpectrumMultipolesCovariance(
            k=k_values,
            covariance=covariance_blocks,
            zeff=zeff,
        )

    return result


def get_PowerSpectrumMultipolesMixingMatrix(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, PowerSpectrumMultipolesMixingMatrix]:
    """
    Reads the mixing matrix for power spectrum multipoles from a FITS file.
    """
    even_multipoles = [0, 2, 4]

    redshifts, nz = check_input(redshifts)
    result: dict[_DictKey, Optional[PowerSpectrumMultipolesMixingMatrix]] = {}

    for i in range(nz):
        for j in range(nz):
            result[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        header, data = read_mixing_matrix(str(path).format(zlab))
        kout = data["BINS_OUTPUT"]["k"]
        kin = {ell: data["BINS_INPUT"]["kp{}".format(ell)] for ell in even_multipoles}
        mixing_matrix_blocks = {
            "ELL_{}-{}".format(ell1, ell2): data["MIXING_MATRIX"][
                "W{}{}".format(ell1, ell2)
            ]
            for ell2 in even_multipoles
            for ell1 in even_multipoles
        }

        try:
            zeff = header["MIXING_MATRIX"]["Z_EFF"]
        except KeyError:
            warn(
                "Effective redshift not specified in FITS file header. Setting it to 0."
            )
            zeff = 0.0

        result[("SPE", "SPE", i, i)] = PowerSpectrumMultipolesMixingMatrix(
            kout=kout,
            kin=kin,
            mixing=mixing_matrix_blocks,
            zeff=zeff,
        )

    return result
