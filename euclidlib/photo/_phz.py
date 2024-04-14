from __future__ import annotations

from collections.abc import Mapping
from os import PathLike
from typing import TYPE_CHECKING

import fitsio  # type: ignore [import-not-found]
import numpy as np

if TYPE_CHECKING:
    from typing import Any
    from numpy.typing import ArrayLike, NDArray


def _hist2dist(x: NDArray[Any], y: NDArray[Any]) -> NDArray[Any]:
    """
    Convert histogram to distribution.
    """
    *dims, n = y.shape
    z = np.zeros((*dims, n + 1), dtype=y.dtype.base)
    np.cumsum(y, axis=-1, out=z[..., 1:])
    result: NDArray[Any] = np.gradient(z, x, axis=-1, edge_order=1)
    return result


# this is the dtype for n(z) data in the Euclid SGS format
NZ_DTYPE = [
    ("BIN_ID", ">i4"),
    ("MEAN_REDSHIFT", ">f4"),
    ("N_Z", ">f4", (3000,)),
]


def write_nz(
    path: str | PathLike[str],
    z: ArrayLike,
    nz: ArrayLike,
    *,
    weight_method: str = "NO_WEIGHT",
    bin_type: str = "TOM_BIN",
    hist: bool = False,
) -> None:
    """
    Write n(z) data in Euclid SGS format.
    """

    z = np.asanyarray(z)
    nz = np.asanyarray(nz)

    if z.ndim != 1:
        raise ValueError("redshift array must be 1D")

    if not hist and z.shape[-1] == nz.shape[-1]:
        pass
    elif hist and z.shape[-1] == nz.shape[-1] + 1:
        nz = _hist2dist(z, nz)
    else:
        raise ValueError("shape mismatch between redshifts and values")

    # PHZ uses a fixed binning scheme with z bins in [0, 6] and dz=0.002
    zbinedges = np.linspace(0.0, 6.0, 3001, dtype=np.float32)

    # compute the combined set of z grid points from data and binning
    zp = np.union1d(z, zbinedges)

    # create the output data in the correct format
    dims = nz.shape[:-1]
    nbin = np.prod(dims, dtype=int)
    out = np.empty(nbin, dtype=NZ_DTYPE)

    # convert every nz into the PHZ format
    for i, index in enumerate(np.ndindex(*dims)):
        # compute mean redshift
        mean_z = np.trapz(z * nz[index], z)

        # interpolate dndz onto the unified grid
        nzp = np.interp(zp, z, nz[index], left=0.0, right=0.0)

        # integrate the number of sources in each bin
        n_z = []
        for z1, z2 in zip(zbinedges, zbinedges[1:]):
            sel = (z1 <= zp) & (zp <= z2)
            n_z.append(np.trapz(nzp[sel], zp[sel]))

        # set the output data row
        out[i] = (i + 1, mean_z, n_z)

    # metadata
    header = {
        "WEIGHT_METHOD": weight_method,
        "BIN_TYPE": bin_type,
        "NBIN": nbin,
    }

    # write output data to FITS
    with fitsio.FITS(path, "rw", clobber=True) as fits:
        fits.write(out, extname="BIN_INFO", header=header)


def redshift_distributions(
    path: str | PathLike[str],
    *,
    ext: str | int | None = "BIN_INFO",
    hist: bool = False,
) -> tuple[NDArray[Any], Mapping[int, NDArray[Any]]]:
    """Read redshift distributions in Euclid format.

    Parameters
    ----------
    path : str
        Path to a FITS file in Euclid format.
    ext : str or int or None, optional
        The FITS extension to read.  If ``None``, the first extension
        with data is used.
    hist : bool, optional
        By default, the histograms are converted to distributions.  If
        true, return the redshift histograms unmodified.

    Returns
    -------
    z : ndarray
        Redshift values.
    nz: dict of int and ndarray
        Dictionary where keys are tomographic bin IDs and values are the
        redshift distributions.

    """

    # data and header from file
    data = fitsio.read(path, ext=ext, lower=True)

    # this is the fixed binning scheme used by PHZ
    z = np.linspace(0.0, 6.0, 3001)

    # check format
    shape = (z.size - 1,)
    if "n_z" not in data.dtype.fields or data.dtype.fields["n_z"][0].shape != shape:
        msg = f"{path}: requires column N_Z of shape {shape}"
        raise ValueError(msg)

    # read each n(z) histogram into a dict
    out = {}
    for row in data:
        bin_id, n_z = row["bin_id"], row["n_z"]
        if not hist:
            n_z = _hist2dist(z, n_z)
        out[bin_id] = n_z

    return z, out
