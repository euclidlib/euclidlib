from __future__ import annotations

# Standard library imports
import os
from os import PathLike
from typing import TYPE_CHECKING, Optional

# Third-party imports
import fitsio  # type: ignore [import-not-found]
import numpy as np
from cosmolib.data import (
    TwoPointCorrelationCartesian,
    TwoPointCorrelationPolar,
    TwoPointCorrelationMultipoles,
    TwoPointCorrelationMultipolesCovariance,
)

# Local library imports
from .._util import writer
from ._common import (
    check_input,
    get_cosmology_from_header,
    read_data_vectors,
    read_and_reshape_covariance_matrix,
    build_2d_correlation,
)

if TYPE_CHECKING:
    from typing import Union

    try:
        from typing import TypeAlias
    except ImportError:
        from typing_extensions import TypeAlias

    _DictKey: TypeAlias = Union[str, int, tuple["_DictKey", ...]]


def twopoint_correlation_cartesian(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, Optional[TwoPointCorrelationCartesian]]:
    """
    Reads the 2-dimensional cartesian 2PCF from a LE3-format fits file and
    returns it as a cosmolib data structure.

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
    results : dict[_DictKey, Optional[TwoPointCorrelationCartesian]]
        Dictionary containing the 2d cartesian 2PCF for each redshift bin pair.
        Keys are tuples of the form ``("SPE", "SPE", i, j)``. For diagonal
        pairs (``i == j``) the value is a `PowerSpectrumMultipoles` instance
        read from the corresponding FITS file, while off-diagonal entries are
        set to ``None``.
    """
    redshifts, nz = check_input(redshifts)
    results: dict[_DictKey, Optional[TwoPointCorrelationCartesian]] = {}

    for i in range(nz):
        for j in range(nz):
            results[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        header, data = read_data_vectors(str(path).format(zlab), "CORRELATION")
        zeff, fiducial_cosmology = get_cosmology_from_header(header)

        scale_1d = data["SCALE_1D"]
        scale_2d = data["SCALE_2D"]
        xi = data["XI"]
        s_perp, s_para, correlation = build_2d_correlation(scale_1d, scale_2d, xi)

        results[("SPE", "SPE", i, i)] = TwoPointCorrelationCartesian(
            s_perp=s_perp,
            s_para=s_para,
            correlation=correlation,
            fiducial_cosmology=fiducial_cosmology,
            zeff=zeff,
        )

    return results


def twopoint_correlation_polar(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, Optional[TwoPointCorrelationPolar]]:
    """
    Reads the 2-dimensional polar 2PCF from a LE3-format fits file and
    returns it as a cosmolib data structure.

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
    results : dict[_DictKey, Optional[TwoPointCorrelationPolar]]
        Dictionary containing the 2d polar 2PCF for each redshift bin pair.
        Keys are tuples of the form ``("SPE", "SPE", i, j)``. For diagonal
        pairs (``i == j``) the value is a `PowerSpectrumMultipoles` instance
        read from the corresponding FITS file, while off-diagonal entries are
        set to ``None``.
    """
    redshifts, nz = check_input(redshifts)
    results: dict[_DictKey, Optional[TwoPointCorrelationPolar]] = {}

    for i in range(nz):
        for j in range(nz):
            results[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        header, data = read_data_vectors(str(path).format(zlab), "CORRELATION")
        zeff, fiducial_cosmology = get_cosmology_from_header(header)

        s_1d = data["SCALE_1D"]
        mu_1d = data["SCALE_2D"]
        correlation_1d = data["XI"]
        s, mu, correlation = build_2d_correlation(s_1d, mu_1d, correlation_1d)

        results[("SPE", "SPE", i, i)] = TwoPointCorrelationPolar(
            s=s,
            mu=mu,
            correlation=correlation,
            fiducial_cosmology=fiducial_cosmology,
            zeff=zeff,
        )

    return results


def twopoint_correlation_multipoles(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, TwoPointCorrelationMultipoles]:
    """
    Reads the 2PCF Legendre multipoles from a LE3-format fits file and returns
    it as a cosmolib data structure.

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
    results : dict[_DictKey, Optional[TwoPointCorrelationMultipoles]]
        Dictionary containing the 2PCF multipoles for each redshift bin pair.
        Keys are tuples of the form ``("SPE", "SPE", i, j)``. For diagonal
        pairs (``i == j``) the value is a `PowerSpectrumMultipoles` instance
        read from the corresponding FITS file, while off-diagonal entries are
        set to ``None``.
    """
    redshifts, nz = check_input(redshifts)
    results: dict[_DictKey, Optional[TwoPointCorrelationMultipoles]] = {}

    for i in range(nz):
        for j in range(nz):
            results[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        header, data = read_data_vectors(str(path).format(zlab), "CORRELATION")
        zeff, fiducial_cosmology = get_cosmology_from_header(header)
        multipoles = np.array([data[f"XI{ell}"] for ell in range(5)])

        results[("SPE", "SPE", i, i)] = TwoPointCorrelationMultipoles(
            s=data["SCALE"],
            multipoles=multipoles,
            fiducial_cosmology=fiducial_cosmology,
            zeff=zeff,
        )

    return results


@writer(twopoint_correlation_multipoles)
def _(
    results: dict[_DictKey, TwoPointCorrelationMultipoles],
    path: Union[str, PathLike[str]],
    *redshifts: str,
) -> None:
    """
    Writes 2PCF Legendre multipoles from a cosmolib data structure to the LE3
    FITS format.

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
        ("SCALE", "f8"),
        ("XI0", "f8"),
        ("XI1", "f8"),
        ("XI2", "f8"),
        ("XI3", "f8"),
        ("XI4", "f8"),
    ]

    redshifts, nz = check_input(redshifts)

    for i, zlab in enumerate(redshifts):
        obj = results[("SPE", "SPE", i, i)]
        out_path = str(path).format(zlab)

        ns = len(obj.s)

        data = np.zeros(ns, dtype=dtype)
        data["SCALE"] = obj.s
        for ell in range(5):
            data[f"XI{ell}"] = obj.multipoles[ell]

        header = {
            "TELESCOP": "EUCLID  ",
            "INSTRUME": "LE3GC   ",
            "RUNTYPE": "AUTO    ",
            "Z_EFF": obj.zeff,
            "STAT": "MULTIPOLE",
            "BIN1TYPE": "LIN     ",
            "BIN1NUM": ns,
            "BIN1MIN": np.min(obj.s) if ns > 0 else 0.0,
            "BIN1MAX": np.max(obj.s) if ns > 0 else 0.0,
            "TUNIT1": "Center of the s bin (Mpc/h)",
            "TUNIT2": "Multipole ell=0",
            "TUNIT3": "Multipole ell=1",
            "TUNIT4": "Multipole ell=2",
            "TUNIT5": "Multipole ell=3",
            "TUNIT6": "Multipole ell=4",
            "COMMENT": "----------- COSMOLOGICAL PARAMETERS USED ----------",
        }
        header.update(obj.fiducial_cosmology)

        if os.path.exists(out_path):
            os.remove(out_path)
        with fitsio.FITS(out_path, "rw") as f:
            f.write(None, header={"EXTNAME": "PRIMARY"})
            f.write(data, header=header, extname="SPECTRUM")


def twopoint_correlation_multipole_covariance(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, Optional[TwoPointCorrelationMultipolesCovariance]]:
    """
    Reads the covariance matrix of the 2PCF Legendre multipoles from a
    LE3-format fits file and returns it as a cosmolib data structure.

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
    results : dict[_DictKey, Optional[TwoPointCorrelationMultipolesCovariance]]
        Dictionary containing the covariance matrix of the 2PCF multipoles for
        each redshift bin pair. Keys are tuples of the form
        ``("SPE", "SPE", i, j)``. For diagonal pairs (``i == j``) the value is
        a `PowerSpectrumMultipoles` instance read from the corresponding FITS
        file, while off-diagonal entries are set to ``None``.
    """
    redshifts, nz = check_input(redshifts)
    results: dict[_DictKey, Optional[TwoPointCorrelationMultipolesCovariance]] = {}

    for i in range(nz):
        for j in range(nz):
            results[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        s_values, covariance_blocks, zeff = read_and_reshape_covariance_matrix(
            path=str(path).format(zlab), type="CORRELATION"
        )

        results[("SPE", "SPE", i, i)] = TwoPointCorrelationMultipolesCovariance(
            s=s_values, covariance=covariance_blocks, zeff=zeff
        )

    return results
