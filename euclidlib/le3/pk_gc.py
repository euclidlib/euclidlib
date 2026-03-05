from __future__ import annotations

import os
from warnings import warn
from os import PathLike

from typing import Optional

import numpy as np
import fitsio
from cosmolib.data import (
    PowerSpectrumMultipoles,
    PowerSpectrumMultipolesCovariance,
    PowerSpectrumMultipolesMixingMatrix,
)
from numpy.typing import NDArray

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


def write_PowerSpectrumMultipoleShifts(
    path: Union[str, PathLike[str]],
    k: NDArray[np.float64],
    pk: dict[_DictKey, NDArray[np.float64]],
    zeff: float,
    cosmology: dict[str, float],
) -> None:
    """
    Writes power spectrum multipoles to a FITS file with a structure
    consistent with Euclid LE3 GC products.
    """
    nrows = len(k)
    dtype = [
        ("K", "f8"),
        ("K_EFF", "f8"),
        ("PK0", "f8"),
        ("PK1", "f8"),
        ("PK2", "f8"),
        ("PK3", "f8"),
        ("PK4", "f8"),
        ("NUM_MOD", "f8"),
    ]
    data = np.zeros(nrows, dtype=dtype)
    data["K"] = data["K_EFF"] = k
    for ell in range(min(5, pk.shape[0])):
        data[f"PK{ell}"] = pk[ell]

    header = {
        "TELESCOP": "EUCLID  ",
        "INSTRUME": "LE3_GC_PK",
        "RUNTYPE": "AUTO    ",
        "Z_EFF": zeff,
        "STAT": None,
        "USE_NBAR": None,
        "MAS": None,
        "MAS_CORR": None,
        "FKP_CORR": None,
        "P_EST": None,
        "INTERLAC": None,
        "SN_CORR": None,
        "SN_VALUE": None,
        "ALPHA": None,
        "SCALE": "LIN     ",
    }
    header.update(cosmology)

    header.update(
        {
            "TUNIT1": "Bin 1d k scale",
            "TUNIT2": "Effective k bin 1D scale",
            "TUNIT3": "Blinding shift for multipole ell=0",
            "TUNIT4": "Blinding shift for multipole ell=1",
            "TUNIT5": "Blinding shift for multipole ell=2",
            "TUNIT6": "Blinding shift for multipole ell=3",
            "TUNIT7": "Blinding shift for multipole ell=4",
            "TUNIT8": "Number of modes averaged",
        }
    )

    if os.path.exists(path):
        os.remove(path)

    with fitsio.FITS(path, "rw") as f:
        f.write(None, header={"EXTNAME": "PRIMARY"})
        f.write(data, header=header, extname="SPECTRUM")
