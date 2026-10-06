"""
Reading routines for SHE shear calibration (multiplicative and additive
bias) products, i.e. ``she.biasParams`` FITS files.

The ``SHEBIASPARAMS`` table has one row per tomographic bin, with
``TOM_BIN_ID = 0`` holding the non-tomographic biases of the full sample.
The biases are stored in two formats:

``PARAMS`` / ``PARAM_COVARS``
    Separate biases per shear component, ``[m1, m2, c1, c2, NaN, NaN]``.
    The two trailing entries are reserved (e.g. for ``m12, m21``).
``EMPIRICAL_PARAMS`` / ``EMPIRICAL_PARAM_COVARS``
    Mean multiplicative bias across components,
    ``[(m1+m2)/2, (m1+m2)/2, c1, c2]``.

In both cases the covariance matrices currently only contain the
variances of each parameter; all other entries are NaN.
"""

from __future__ import annotations

from os import PathLike
from typing import TYPE_CHECKING

import fitsio  # type: ignore [import-not-found]
import numpy as np

if TYPE_CHECKING:
    from typing import Any


# names of the bias parameters read from the products, in file order
_NAMES = ("m1", "m2", "c1", "c2")


def multiplicative_shear_bias(
    path: str | PathLike[str],
    *,
    ext: str | int = "SHEBIASPARAMS",
    empirical: bool = True,
) -> dict[int, dict[str, Any]]:
    """Read shear calibration (m and c bias) parameters in Euclid format.

    Reads a ``she.biasParams`` product, as produced e.g. for the LensMC
    and MetaCal shear measurement methods.

    Parameters
    ----------
    path : str
        Path to a FITS file in Euclid format.
    ext : str or int, optional
        The FITS extension to read.
    empirical : bool, optional
        If true (the default), read the ``EMPIRICAL_PARAMS`` and
        ``EMPIRICAL_PARAM_COVARS`` columns, where ``m1 = m2`` is the mean
        multiplicative bias across the two shear components.  If false,
        read the ``PARAMS`` and ``PARAM_COVARS`` columns with separate
        ``m1`` and ``m2``; only the first four parameters are returned.

    Returns
    -------
    m_bias : dict of int and dict
        Dictionary where keys are tomographic bin IDs and values are
        dictionaries of shear bias parameters, so that e.g.
        ``m_bias[1]["m1"]`` is the ``m1`` bias of the first tomographic
        bin.  Bin ID 0 holds the parameters for the full
        (non-tomographic) sample.  The keys of each bin are:

        ``"m1"``, ``"m2"``
            Multiplicative biases of the two shear components.
        ``"c1"``, ``"c2"``
            Additive biases of the two shear components.
        ``"cov"``
            Covariance matrix of ``(m1, m2, c1, c2)``, shape ``(4, 4)``.
            Entries that are not provided are NaN (currently all
            off-diagonal entries).
        ``"num_objects"``
            Number of objects used for the calibration, or -1 if not
            available (e.g. for the full sample).

        The bias model is :math:`g_i^{\\rm obs} = (1 + m_i)
        g_i^{\\rm true} + c_i` for the two shear components
        :math:`i = 1, 2`.

    """

    if empirical:
        col_par, col_cov = "empirical_params", "empirical_param_covars"
    else:
        col_par, col_cov = "params", "param_covars"

    data = fitsio.read(path, ext=ext, lower=True)

    # check format
    fields = data.dtype.fields
    for col in ("tom_bin_id", "num_objects", col_par, col_cov):
        if col not in fields:
            msg = f"{path}: requires column {col.upper()}"
            raise ValueError(msg)
    npar = int(np.prod(fields[col_par][0].shape))
    if npar < len(_NAMES) or int(np.prod(fields[col_cov][0].shape)) != npar**2:
        msg = (
            f"{path}: column {col_par.upper()} must have at least {len(_NAMES)} "
            f"entries and {col_cov.upper()} must be its square"
        )
        raise ValueError(msg)

    npar_out = len(_NAMES)
    m_bias = {}
    for row in data:
        par = np.asarray(row[col_par], dtype=float).reshape(npar)
        cov = np.asarray(row[col_cov], dtype=float).reshape(npar, npar)
        m_bias[int(row["tom_bin_id"])] = {
            **{name: float(par[i]) for i, name in enumerate(_NAMES)},
            "cov": cov[:npar_out, :npar_out].copy(),
            "num_objects": int(row["num_objects"]),
        }

    return m_bias
