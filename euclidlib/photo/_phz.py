from __future__ import annotations

from collections.abc import Mapping
from os import PathLike
from datetime import datetime
from typing import TYPE_CHECKING

import fitsio  # type: ignore [import-not-found]
import numpy as np

from .._util import writer

if TYPE_CHECKING:
    from typing import Any
    from numpy.typing import ArrayLike, NDArray

if np.lib.NumpyVersion(np.__version__) >= "2.0.0b1":
    trapezoid = np.trapezoid
else:
    trapezoid = np.trapz  # type: ignore


def _hist2dist(x: NDArray[Any], y: NDArray[Any]) -> NDArray[Any]:
    """
    Convert histogram to distribution.
    """
    *dims, n = y.shape
    z = np.zeros((*dims, n + 1), dtype=y.dtype.base)
    np.cumsum(y, axis=-1, out=z[..., 1:])
    result: NDArray[Any] = np.gradient(z, x, axis=-1, edge_order=1)
    return result


def redshift_distributions(
    path: str | PathLike[str],
    *,
    ext: str | int | None = None,
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


@writer(redshift_distributions)
def _(
    path: str | PathLike[str],
    z: ArrayLike,
    nz: dict,
    *,
    weight_method: str = "NO_WEIGHT",
    bin_type: str = "TOM_BIN",
    hist: bool = False,
) -> None:
    """
    Write n(z) data in Euclid SGS format.  Supports both distributions
    (when *hist* is false, the default) and histograms (when *hist* is
    true).
    Parameters
    ----------
    path : str or PathLike[str]
        Path to the output FITS file.
    z : array-like
        Redshift values.  Must be 1D.
    nz : dictionary
        Redshift distributions or histograms as read by euclidlib.  Must be at least 1D.
        If *hist* is true, the last dimension must match the length of
        *z* minus one.  If *hist* is false, the last dimension must
        match the length of *z*.
    weight_method : str, optional
        Weight method used for the binning.  Default is "NO_WEIGHT".
    bin_type : str, optional
        Type of binning used for the redshift distributions.  Default is
        "TOM_BIN".
    hist : bool, optional
        If true, the input *nz* is interpreted as a histogram.  If
        false, it is interpreted as a distribution.  Default is false.
    """

    z = np.asanyarray(z)
    nz = np.asanyarray(np.array(list(nz.values())).T)
    nz = nz.T  # Flip axes if shape is (z, bins)
    nz = nz.reshape(-1, nz.shape[-1])
    nbin = nz.shape[0]
    # PHZ uses a fixed binning scheme with z bins in [0, 6] and dz=0.002
    zbinedges = np.linspace(0.0, 6.0, 3001)

    # create the output data in the correct format
    out = np.empty(
        nbin,
        dtype=[
            ("BIN_ID", ">i4"),
            ("MEAN_REDSHIFT", ">f4"),
            ("N_Z", ">f4", (3000,)),
        ],
    )

    # set increasing bin IDs
    out["BIN_ID"] = np.arange(1, nbin + 1)

    # convert every nz into the PHZ format
    if hist:
        # shorthand for the left and right z boundaries, respectively
        zl, zr = z[:-1], z[1:]
        # compute the mean redshifts
        mid_z = (zl + zr) / 2
        nz_bins = 0.5 * (nz[:, :-1] + nz[:, 1:]) 
        mid_z = mid_z[-nz_bins.shape[1]:]
        out["MEAN_REDSHIFT"] = np.sum(mid_z * nz_bins, axis=1) / np.sum(nz, axis=1)
        # compute resummed bin counts
        for j, (z1, z2) in enumerate(zip(zbinedges, zbinedges[1:])):
            frac = (np.clip(z2, zl, zr) - np.clip(z1, zl, zr)) / (zr - zl)
            frac = frac[-nz_bins.shape[1]:]
            out["N_Z"][:, j] = np.dot(nz_bins, frac)
    else:
        # integrate the n(z) over each histogram bin
        # compute mean redshifts
        out["MEAN_REDSHIFT"] = trapezoid(z * nz, z, axis=-1) / trapezoid(nz, z, axis=-1)

        # compute the combined set of z grid points from data and binning
        zp = np.union1d(z, zbinedges)

        # integrate over each bin
        for i in range(nbin):
            # interpolate dndz onto the unified grid
            nzp = np.interp(zp, z, nz[i], left=0.0, right=0.0)

            # integrate the distribution over each bin
            for j, (z1, z2) in enumerate(zip(zbinedges, zbinedges[1:])):
                sel = (z1 <= zp) & (zp <= z2)
                out["N_Z"][i, j] = trapezoid(nzp[sel], zp[sel])

    # metadata
    timestamp = datetime.now().isoformat(timespec="seconds")
    header = {
        "WEIGHT_METHOD": weight_method,
        "BIN_TYPE": bin_type,
        "NBIN": nbin,
        "HISTORY": f"Generated by euclidlib on {timestamp}",
    }

    # write output data to FITS
    with fitsio.FITS(path, "rw", clobber=True) as fits:
        fits.write(out, extname="BIN_INFO", header=header)
