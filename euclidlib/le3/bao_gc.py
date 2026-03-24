from __future__ import annotations

from os import PathLike
from typing import TYPE_CHECKING, Optional
from cosmolib.data import (
    BaryonAcousticOscillations,
    BaryonAcousticOscillationsCovariance,
)

from ._common import (
    check_input,
    get_cosmology_from_header,
    read_data_vectors,
    read_covariance_matrix,
)

from numpy.typing import NDArray

if TYPE_CHECKING:
    from typing import Any, Union

    try:
        from typing import TypeAlias
    except ImportError:
        from typing_extensions import TypeAlias

    _DictKey: TypeAlias = Union[str, int, tuple["_DictKey", ...]]


def BAO_alphas(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, BaryonAcousticOscillations]:
    """
    Returns alphas from BAO in cloe-compatible euclidlib data format

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
    results : dict[_DictKey, Optional[BaryonAcousticOscillations]]
    """
    redshifts, nz = check_input(redshifts)
    result: dict[_DictKey, Optional[BaryonAcousticOscillations]] = {}

    for i in range(nz):
        for j in range(nz):
            result[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        header, data = read_data_vectors(str(path).format(zlab), "BAO_ALPHAS")
        zeff, fiducial_cosmology = get_cosmology_from_header(header)

        alpha_par = data["ALPHA_PAR"]
        alpha_perp = data["ALPHA_PER"]
        alpha_iso = data["ALPHA_ISO"]
        alpha_ap = data["ALPHA_AP"]

        result[("SPE", "SPE", i, i)] = BaryonAcousticOscillations(
            alpha_par=alpha_par,
            alpha_perp=alpha_perp,
            alpha_iso=alpha_iso,
            alpha_ap=alpha_ap,
            fiducial_cosmology=fiducial_cosmology,
            zeff=zeff,
        )

    return result


def BAO_alphas_covariance(
    path: Union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, BaryonAcousticOscillationsCovariance]:
    """
    Returns the covariance for BAO alphas in cloe-compatible euclidlib data format

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
    results : dict[_DictKey, Optional[BaryonAcousticOscillationsCovariance]]
        Dictionary containing the covariance matrix of the BAO alphas for each
        redshift bin pair. Keys are tuples of the form ``("SPE", "SPE", i, j)``.
        For diagonal pairs (``i == j``) the value is a
        `BaryonAcousticOscillationsCovariance` instance read from the
        corresponding FITS file, while off-diagonal entries are set to ``None``.
    """
    redshifts, nz = check_input(redshifts)
    results: dict[_DictKey, Optional[BaryonAcousticOscillationsCovariance]] = {}

    for i in range(nz):
        for j in range(nz):
            results[("SPE", "SPE", i, j)] = None

    for i, zlab in enumerate(redshifts):
        header, data = read_covariance_matrix(
            path=str(path).format(zlab),
        )
        zeff, _ = get_cosmology_from_header(header, get_fiducial=False)
        correction_factor = header["CORR_FAC"]

        observables = sorted(set(data["ALPHA-I"]).union(set(data["ALPHA-J"])))

        covariance_blocks: dict[str, NDArray[Any]] = {}
        for row in data:
            key = f"{row['ALPHA-I']}-{row['ALPHA-J']}"
            covariance_blocks[key] = row["COVARIANCE"]

        results[("SPE", "SPE", i, i)] = BaryonAcousticOscillationsCovariance(
            observables=observables,
            covariance=covariance_blocks,
            correction_factor=correction_factor,
            zeff=zeff,
        )

    return results
