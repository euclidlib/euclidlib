from __future__ import annotations

from collections.abc import Mapping
from os import PathLike
from typing import Any

import fitsio  # type: ignore [import-not-found]
import numpy as np
from numpy.typing import NDArray


def redshift_distributions(
    path: str | PathLike[str],
    *,
    ext: str | int | None = None,
) -> tuple[NDArray[Any], Mapping[int, NDArray[Any]]]:
    """Read redshift distributions in Euclid format.

    Parameters
    ----------
    path : str
        Path to a FITS file in Euclid format.
    ext : str or int, optional
        The FITS extension to read.  By default, the first extension
        with data is used.

    Returns
    -------
    z : ndarray
        Redshift bin edges of the redshift distribution histograms.
    nzdict : dict of int and ndarray
        Dictionary where keys are tomographic bin IDs and values are the
        histograms of the redshift distributions.

    """

    # data and header from file
    data = fitsio.read(path, ext=ext, lower=True)

    # this is the fixed binning scheme used by PHZ
    z = np.arange(0.0, 6.001, 0.002)

    # check format
    shape = (z.size - 1,)
    if "n_z" not in data.dtype.fields or data.dtype.fields["n_z"][0].shape != shape:
        msg = f"{path}: requires column N_Z of shape {shape}"
        raise ValueError(msg)

    # read each n(z) histogram into a dict
    nz_dict = {}
    for row in data:
        bin_id, n_z = row["bin_id"], row["n_z"]
        nz_dict[bin_id] = n_z

    return z, nz_dict
