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
    read_data_vectors,
    read_covariance_data,
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


def get_PowerSpectrumMultipoles(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, PowerSpectrumMultipoles]:
    """
    Returns power spectrum data in the cloe-compatible euclidlib data format
    """
    nz = len(redshifts)
    result: dict[_DictKey, Optional[PowerSpectrumMultipoles]] = {}

    for i in range(nz):
        for j in range(nz):
            result[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        header, data = read_data_vectors(path.format(zlab), "SPECTRUM")
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
    nz = len(redshifts)
    result: dict[_DictKey, Optional[PowerSpectrumMultipolesCovariance]] = {}

    for i in range(nz):
        for j in range(nz):
            result[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        header, data = read_covariance_data(path.format(zlab))
        zeff = get_cosmology_from_header(header, get_fiducial=False)
        mask = np.isin(data["MULTIPOLE-I"], range(5)) & np.isin(
            data["MULTIPOLE-J"], range(5)
        )
        data = data[mask]

        scale_prefix = "k"
        scale_i_col = f"{scale_prefix.upper()}I"

        k_values = np.unique(data[scale_i_col])
        nk = len(k_values)

        covariance_blocks: dict[str, NDArray[Any]] = {}
        for ell_i in range(5):
            for ell_j in range(5):
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

    nz = len(redshifts)
    result: dict[_DictKey, Optional[PowerSpectrumMultipolesMixingMatrix]] = {}

    for i in range(nz):
        for j in range(nz):
            result[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        with fitsio.FITS(path.format(zlab)) as fits:
            kout_data = fits["BINS_OUTPUT"].read()
            kout = kout_data["k"]

            kin_data = fits["BINS_INPUT"].read()
            kin = {ell: kin_data["kp{}".format(ell)] for ell in even_multipoles}

            mixing_data = fits["MIXING_MATRIX"].read()
            header = fits["MIXING_MATRIX"].read_header()

            mixing_matrix_blocks = {
                "ELL_{}-{}".format(ell1, ell2): mixing_data["W{}{}".format(ell1, ell2)]
                for ell2 in even_multipoles
                for ell1 in even_multipoles
            }

            try:
                zeff = header["Z_EFF"]
            except KeyError:
                warn(
                    "Effective redshift not specified in FITS file header. "
                    "Setting it to 0."
                )
                zeff = 0.0

            result[("SPE", "SPE", i, i)] = PowerSpectrumMultipolesMixingMatrix(
                kout=kout,
                kin=kin,
                mixing=mixing_matrix_blocks,
                zeff=zeff,
            )

    return result
