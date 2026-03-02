from __future__ import annotations

from os import PathLike

from cosmolib.data import BaryonAcousticOscillations


from ._common import (
    check_input,
    get_cosmology_from_header,
    read_data_vectors,
)

TYPE_CHECKING = True
if TYPE_CHECKING:
    from typing import Union

    try:
        from typing import TypeAlias
    except ImportError:
        from typing_extensions import TypeAlias

    _DictKey: TypeAlias = Union[str, int, tuple["_DictKey", ...]]


def get_BAO_alphas(
    path: union[str, PathLike[str]], *redshifts: str
) -> dict[_DictKey, TwoPointCorrelationCartesian]:
    """
    Returns alphas from BAO in cloe-compatible euclidlib data format
    """
    redshifts, nz = check_input(redshifts)
    result: dict[_DictKey] = {}

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
