from __future__ import annotations

from collections.abc import Mapping
from os import PathLike
from typing import TYPE_CHECKING

import fitsio  # type: ignore [import-not-found]
import numpy as np

from ._util import writer

if TYPE_CHECKING:
    from typing import Any, Literal
    from numpy.typing import NDArray

if np.lib.NumpyVersion(np.__version__) >= "2.0.0b1":
    trapezoid = np.trapezoid
else:
    trapezoid = np.trapz  # type: ignore


def redshift_distributions(
    path: str | PathLike[str],
    *,
    ext: str | int | None = None,
    version: Literal["q1", "dr1"] = "dr1",
) -> tuple[NDArray[Any], Mapping[int, NDArray[Any]]]:
    """Read redshift distributions in Euclid format.

    Parameters
    ----------
    path : str
        Path to a FITS file in Euclid format.

    Returns
    -------
    z : ndarray
        Redshift values.
    nz : dict of int and ndarray
        Dictionary where keys are tomographic bin IDs and values are the
        redshift distributions.

    Other Parameters
    ----------------
    ext : str or int or None, optional
        The FITS extension to read.  If ``None``, the first extension
        with data is used.
    version : {"q1", "dr1"}, optional
        File format version.  By default, the latest available version is used.

    """

    # data and header from file
    data, hdr = fitsio.read(path, ext=ext, lower=True, header=True)

    # check for new format
    if version == "q1":
        # this is the fixed binning scheme used by PHZ
        z = np.linspace(0.0, 6.0, 3000, endpoint=False)

        # check format
        shape = (z.size,)
        if "n_z" not in data.dtype.fields or data.dtype.fields["n_z"][0].shape != shape:
            msg = f"{path}: requires column N_Z of shape {shape}"
            raise ValueError(msg)

        # load n(z) histogram as distribution
        nz = {row["bin_id"]: row["n_z"] for row in data}

    elif version == "dr1":
        # get redshift grid from file
        z_step = hdr["Z_STEP"]
        z_size = hdr["Z_STEP_NUMBER"]
        z = z_step * np.arange(z_size)

        # read n(z) and cut to size
        nz = {row["bin_id"]: row["n_z"][:z_size] for row in data}

    else:
        raise ValueError(f"invalid version: {version}")

    return z, nz


@writer(redshift_distributions)
def _(
    path: str | PathLike[str],
    z: NDArray[Any],
    nz: Mapping[int, NDArray[Any]],
    *,
    weight_method: str = "NO_WEIGHT",
    bin_type: str = "TOM_BIN",
    version: Literal["dr1"] = "dr1",
) -> None:
    """
    Write n(z) data in Euclid SGS format.

    Parameters
    ----------
    path : str
        Path to a FITS file in Euclid format.
    z : ndarray
        Redshift values.
    nz : dict of int and ndarray
        Dictionary where keys are tomographic bin IDs and values are the
        redshift distributions.

    Other Parameters
    ----------------
    weight_method : str, optional
        Set weight method in FITS header.
    bin_type : str, optional
        Set bin type in FITS header.
    version : {"dr1"}, optional
        File format version.  By default, the latest available version is used.

    """

    if version not in ["dr1"]:
        raise ValueError(f"invalid version: {version}")

    if z.ndim != 1:
        raise ValueError("z array must be 1D")
    if z.size < 2:
        raise ValueError("z array must contain more than 1 value")

    # figure out the redshift spacing
    z_step = z[1] - z[0]
    if not np.allclose(z, z_step * np.arange(z.size)):
        raise ValueError("z array is not a regular grid")

    # get number of bins from nz
    nbin = len(nz)

    # create the output data in the correct format
    out = np.empty(
        nbin,
        dtype=[
            ("BIN_ID", ">i4"),
            ("MEAN_REDSHIFT", ">f4"),
            ("MEAN_REDSHIFT_ERR", ">f4"),
            ("VAR_REDSHIFT", ">f4"),
            ("VAR_REDSHIFT_ERR", ">f4"),
            ("N_Z", ">f8", (z.size,)),
        ],
    )

    # go through nz and set each as row in output
    for i, (bin_id, dist) in enumerate(nz.items()):
        # check array shape
        if dist.shape != z.shape:
            raise ValueError(f"shape mismatch for bin {bin_id}")

        # compute mean and variance of distribution
        norm = trapezoid(dist, z, axis=-1)
        mean = trapezoid(z * dist, z, axis=-1) / norm
        var = trapezoid((z - mean) ** 2 * dist, z, axis=-1) / norm

        # set data for row
        out[i]["BIN_ID"] = bin_id
        out[i]["MEAN_REDSHIFT"] = mean
        out[i]["MEAN_REDSHIFT_ERR"] = 0.0
        out[i]["VAR_REDSHIFT"] = var
        out[i]["VAR_REDSHIFT_ERR"] = 0.0
        out[i]["N_Z"] = dist

    # metadata
    header = {
        "WEIGHT_METHOD": weight_method,
        "BIN_TYPE": bin_type,
        "NBIN": nbin,
        "Z_STEP": z_step,
        "Z_STEP_NUMBER": z.size,
    }

    # write output data to FITS
    with fitsio.FITS(path, "rw", clobber=True) as fits:
        fits.write(None)
        fits.write_table(out, extname="BIN_INFO", header=header)
