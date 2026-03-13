from __future__ import annotations

# Standard library imports
import os
from os import PathLike
from typing import TYPE_CHECKING, Optional
from warnings import warn

# Third-party imports
import fitsio  # type: ignore [import-not-found]
import numpy as np
from cosmolib.data import (
    PowerSpectrumMultipoles,
    PowerSpectrumMultipolesCovariance,
    PowerSpectrumMultipolesMixingMatrix,
)

# Local library imports
from .._util import writer
from ._common import (
    check_input,
    get_cosmology_from_header,
    read_data_vectors,
    read_and_reshape_covariance_matrix,
    read_mixing_matrix,
)

if TYPE_CHECKING:
    from typing import Union

    try:
        from typing import TypeAlias
    except ImportError:
        from typing_extensions import TypeAlias

    _DictKey: TypeAlias = Union[str, int, tuple["_DictKey", ...]]


def power_spectrum_multipoles(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, Optional[PowerSpectrumMultipoles]]:
    """
    Reads the power spectrum Legendre multipoles from a LE3-format fits file
    and returns it as a cosmolib data structure.

    Parameters
    ----------
    path : str or PathLike
        Path to the input FITS file. If the string contains a `{}` placeholder,
        it will be formatted using the provided redshift labels.
    *redshifts: str
        Redshift labels used to format the input file name. Each value replaces
        the `{}` placeholder in `path`, generating one input file per redshift.
        If no labels are provided, `path` is assumed to be already complete.

    Returns
    -------
    results : dict[_DictKey, Optional[PowerSpectrumMultipoles]]
        Dictionary containing the power spectrum multipoles for each redshift
        bin pair. Keys are tuples of the form ``("SPE", "SPE", i, j)``. For
        diagonal pairs (``i == j``) the value is a `PowerSpectrumMultipoles`
        instance read from the corresponding FITS file, while off-diagonal
        entries are set to ``None``.
    """
    redshifts, nz = check_input(redshifts)
    results: dict[_DictKey, Optional[PowerSpectrumMultipoles]] = {}

    for i in range(nz):
        for j in range(nz):
            results[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        header, data = read_data_vectors(str(path).format(zlab), "SPECTRUM")
        zeff, fiducial_cosmology = get_cosmology_from_header(header)
        nbar = 1.0 / header["SN_VALUE"]
        multipoles = np.array([data[f"PK{ell}"] for ell in range(5)])

        results[("SPE", "SPE", i, i)] = PowerSpectrumMultipoles(
            k=data["K"],
            keff=data["K_EFF"],
            Nmodes=data["NUM_MOD"],
            multipoles=multipoles,
            fiducial_cosmology=fiducial_cosmology,
            zeff=zeff,
            nbar=nbar,
            Psn=header["SN_VALUE"],
        )

    return results


@writer(power_spectrum_multipoles)
def _(
    results: dict[_DictKey, PowerSpectrumMultipoles],
    path: Union[str, PathLike[str]],
    *redshifts: str,
) -> None:
    """
    Writes power spectrum Legendre multipoles from a cosmolib data structure to
    the LE3 FITS format.

    Parameters
    ----------
    results : dict
        Dictionary mapping keys in the euclidlib format ``("SPE", "SPE", i, j)``
        to cosmolib objects containing the data to write to FITS files.
    path : str or PathLike
        Path to the output FITS file. If the string contains a `{}` placeholder,
        it will be formatted using the provided redshift labels.
    *redshifts: str
        Redshift labels used to format the output file name. Each value replaces
        the `{}` placeholder in `path`, generating one output file per redshift.
        If no labels are provided, `path` is assumed to be already finalised.
    """
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

    redshifts, nz = check_input(redshifts)

    for i, zlab in enumerate(redshifts):
        obj = results[("SPE", "SPE", i, i)]
        out_path = str(path).format(zlab)

        nk = len(obj.k)
        data = np.zeros(nk, dtype=dtype)

        data["K"] = obj.k
        data["K_EFF"] = obj.keff if obj.keff is not None else obj.k
        for ell in range(5):
            data[f"PK{ell}"] = obj.multipoles[ell]
        if obj.Nmodes is not None:
            data["NUM_MOD"] = obj.Nmodes

        header = [
            {"name": "COMMENT", "value": " "},
            {"name": "COMMENT", "value": "----------- PowerSpectrum HDU ----------"},
            {"name": "COMMENT", "value": " "},
            {"name": "TELESCOP", "value": "EUCLID  "},
            {"name": "INSTRUME", "value": "cloelib + euclidlib"},
            {"name": "RUNTYPE", "value": "AUTO    "},
            {"name": "TUNIT1", "value": "Center of the k bin [h/Mpc]"},
            {"name": "TUNIT2", "value": "Effective k [h/Mpc]"},
            {"name": "TUNIT3", "value": "Multipole ell=0 [(Mpc/h)^3]"},
            {"name": "TUNIT4", "value": "Multipole ell=1 [(Mpc/h)^3]"},
            {"name": "TUNIT5", "value": "Multipole ell=2 [(Mpc/h)^3]"},
            {"name": "TUNIT6", "value": "Multipole ell=3 [(Mpc/h)^3]"},
            {"name": "TUNIT7", "value": "Multipole ell=4 [(Mpc/h)^3]"},
            {"name": "TUNIT8", "value": "Number of modes"},
            {"name": "COMMENT", "value": " "},
            {"name": "COMMENT", "value": "----------- Spectrum parameters ----------"},
            {"name": "COMMENT", "value": " "},
            {"name": "STAT", "value": "MULTIPOLES"},
            {"name": "Z_EFF", "value": obj.zeff},
            {"name": "NBAR", "value": obj.nbar},
            {"name": "SN_VALUE", "value": obj.Psn},
            {"name": "SCALE", "value": "LINEAR  "},
            {"name": "N_BIN", "value": nk},
            {"name": "K_MIN", "value": np.min(obj.k)},
            {"name": "K_MAX", "value": np.max(obj.k)},
            {"name": "DELTA_K", "value": obj.k[1] - obj.k[0]},
            {"name": "COMMENT", "value": " "},
            {"name": "COMMENT", "value": "----------- Fiducial cosmology ----------"},
            {"name": "COMMENT", "value": " "},
        ]
        for key, value in obj.fiducial_cosmology.items():
            header.append({"name": key, "value": value})

        if os.path.exists(out_path):
            os.remove(out_path)
        with fitsio.FITS(out_path, "rw") as f:
            f.write(None, header={"EXTNAME": "PRIMARY"})
            f.write(data, header=header, extname="SPECTRUM")


def power_spectrum_multipole_covariance(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, Optional[PowerSpectrumMultipolesCovariance]]:
    """
    Reads the covariance matrix of the power spectrum Legendre multipoles from
    a LE3-format fits file and returns it as a cosmolib data structure.

    Parameters
    ----------
    path : str or PathLike
        Path to the input FITS file. If the string contains a `{}` placeholder,
        it will be formatted using the provided redshift labels.
    *redshifts: str
        Redshift labels used to format the input file name. Each value replaces
        the `{}` placeholder in `path`, generating one input file per redshift.
        If no labels are provided, `path` is assumed to be already complete.

    Returns
    -------
    results : dict[_DictKey, Optional[PowerSpectrumMultipolesCovariance]]
        Dictionary containing the covariance matrix of the power spectrum
        multipoles for each redshift bin pair. Keys are tuples of the form
        ``("SPE", "SPE", i, j)``. For diagonal pairs (``i == j``) the value is
        a `PowerSpectrumMultipolesCovariance` instance read from the
        corresponding FITS file, while off-diagonal entries are set to ``None``.
    """
    redshifts, nz = check_input(redshifts)
    results: dict[_DictKey, Optional[PowerSpectrumMultipolesCovariance]] = {}

    for i in range(nz):
        for j in range(nz):
            results[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        k_values, covariance_blocks, zeff = read_and_reshape_covariance_matrix(
            path=str(path).format(zlab),
            type="SPECTRUM",
        )
        results[("SPE", "SPE", i, i)] = PowerSpectrumMultipolesCovariance(
            k=k_values, covariance=covariance_blocks, zeff=zeff
        )

    return results


def power_spectrum_multipole_mixing_matrix(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, Optional[PowerSpectrumMultipolesMixingMatrix]]:
    """
    Reads the mixing matrix of the power spectrum Legendre multipoles from
    a LE3-format fits file and returns it as a cosmolib data structure.

    Parameters
    ----------
    path : str or PathLike
        Path to the input FITS file. If the string contains a `{}` placeholder,
        it will be formatted using the provided redshift labels.
    *redshifts: str
        Redshift labels used to format the input file name. Each value replaces
        the `{}` placeholder in `path`, generating one input file per redshift.
        If no labels are provided, `path` is assumed to be already complete.

    Returns
    -------
    results : dict[_DictKey, Optional[PowerSpectrumMultipolesMixingMatrix]]
        Dictionary containing the mixing matrix of the power spectrum multipoles
        for each redshift bin pair. Keys are tuples of the form
        ``("SPE", "SPE", i, j)``. For diagonal pairs (``i == j``) the value is
        a `PowerSpectrumMultipolesMixingMatrix` instance read from the
        corresponding FITS file, while off-diagonal entries are set to ``None``.
    """
    even_multipoles = [0, 2, 4]

    redshifts, nz = check_input(redshifts)
    results: dict[_DictKey, Optional[PowerSpectrumMultipolesMixingMatrix]] = {}

    for i in range(nz):
        for j in range(nz):
            results[("SPE", "SPE", i, j)] = None

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

        results[("SPE", "SPE", i, i)] = PowerSpectrumMultipolesMixingMatrix(
            kout=kout, kin=kin, mixing=mixing_matrix_blocks, zeff=zeff
        )

    return results
