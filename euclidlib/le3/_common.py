from __future__ import annotations

from warnings import warn
from os import PathLike
import os
import fitsio  # type: ignore [import-not-found]
import numpy as np


from numpy.typing import NDArray

TYPE_CHECKING = True
if TYPE_CHECKING:
    from typing import Any, Dict, Tuple, Union

    try:
        from typing import TypeAlias
    except ImportError:
        from typing_extensions import TypeAlias

    _DictKey: TypeAlias = Union[str, int, tuple["_DictKey", ...]]


def _verify_input_file(
    path: Union[str, PathLike[str]], check_extra_hdu: bool = False
) -> None:
    """
    Verifies that input file is a compatible LE3-GC fits file
    """
    hdul_target_length = 3 if check_extra_hdu else 2

    if not any((isinstance(path, str), isinstance(path, PathLike))):
        raise TypeError(
            "Provided fits file name must be a string or a PathLike object."
        )

    if not os.path.isfile(path):
        raise FileNotFoundError(f"Could not find file '{path}'.")

    try:
        with fitsio.FITS(path) as fits_input:
            assert len(fits_input) > hdul_target_length - 1
            assert "EXTNAME" in _get_hdu_header(fits_input[hdul_target_length - 1])
    except OSError:
        raise ValueError("Provided file is not a valid fits file.")
    except AssertionError:
        raise ValueError(
            "Provided fits file does not match the structure of a valid LE3-GC product."
        )
    except Exception:
        raise RuntimeError("Invalid file provided.")


def _get_hdu_header(hdu: fitsio.TableHDU) -> Dict[str, Any]:
    """
    Reads header from fits file HDU
    """
    head: Dict[str, Any] = dict(hdu.read_header())
    return head


def _get_hdu_data(hdu: fitsio.TableHDU) -> NDArray[Any]:
    """
    Reads data from fits file HDU
    """
    arr: NDArray[Any] = hdu.read()
    return arr


def get_cosmology_from_header(
    header: Dict[str, Any], get_fiducial: bool = True
) -> Tuple[float, Dict[str, float]]:
    """
    Extracts redshift and fiducial cosmology from a FITS header.
    """
    try:
        z_eff = header["Z_EFF"]
    except KeyError:
        warn("Effective redshift not specified in fits file header. Setting it to 0.")
        z_eff = 0.0

    if not get_fiducial:
        return z_eff
    else:
        fiducial_cosmology = {
            "OMEGA_M": header["OMEGA_M"],
            "OMEGA_R": header["OMEGA_R"],
            "OMEGA_B": header["OMEGA_B"],
            "OMEGA_V": header["OMEGA_V"],
            "OMEGA_K": header["OMEGA_K"],
            "HUBBLE": header["HUBBLE"],
            "INDEX_N": header["INDEX_N"],
            "SIGMA_8": header["SIGMA_8"],
            "W_STATE": header["W_STATE"],
            "N_EFF": header["N_EFF"],
            "T_CMB": header["T_CMB"],
        }
        return z_eff, fiducial_cosmology


def _read_le3_data(
    path: Union[str, PathLike[str]], ext_name: str, check_extra_hdu: bool = False
) -> Tuple[Dict[str, Any], NDArray[Any]]:
    """
    Reads data from a Euclid LE3 fits file
    """
    not_found_message = f"Invalid fits file provided, cannot locate {ext_name} data."

    _verify_input_file(path, False)  # Always check for at least 2 HDUs

    with fitsio.FITS(path) as fits_input:
        header = _get_hdu_header(fits_input[1])

        if ext_name in header["EXTNAME"]:
            data = _get_hdu_data(fits_input[1])
            return header, data

        if check_extra_hdu:
            warn(
                f"\n\033[0;31m[!]\033[0m Found {header['EXTNAME']}, not {ext_name}. "
                "Falling back to next HDU."
            )
            _verify_input_file(path, True)  # Now check for 3 HDUs
            header = _get_hdu_header(fits_input[2])
            if ext_name in header["EXTNAME"]:
                data = _get_hdu_data(fits_input[2])
                return header, data

    raise ValueError(not_found_message)


def _read_covariance_data(
    path: Union[str, PathLike[str]],
) -> Tuple[Dict[str, Any], NDArray[Any]]:
    """
    Reads covariance matrix data from Euclid LE3-CM-GC fits file
    """
    data = None
    header = None

    with fitsio.FITS(path) as fits_input:
        for hdu in fits_input:
            extname = hdu.get_extname() if hdu.get_extname() else ""
            if "COVARIANCE" in extname:
                data = hdu.read()
                header = _get_hdu_header(hdu)

    if data is None:
        raise ValueError("HDU 'COVARIANCE' not found in file.")

    return header, data


def build_covariance_matrix(
    data: NDArray[Any],
    mode: str,
) -> Tuple[NDArray[Any], NDArray[Any]]:
    """
    Builds the full covariance matrix from the table data.
    """
    even_multipoles = [0, 2, 4]
    mask = np.isin(data["MULTIPOLE-I"], even_multipoles) & np.isin(
        data["MULTIPOLE-J"], even_multipoles
    )
    data = data[mask]

    if mode == "PS":
        scale_prefix = "k"
    elif mode == "TPCF":
        scale_prefix = "s"
    else:
        raise ValueError(
            "Mode for building covariance must be picked from ('PS', 'TPCF')."
        )

    scale_i_col = f"{scale_prefix.upper()}I"
    scale_j_col = f"{scale_prefix.upper()}J"

    all_unique_scale_values = np.unique(data[scale_i_col])
    scale_to_idx = {s_val: i for i, s_val in enumerate(all_unique_scale_values)}
    n_scale = len(all_unique_scale_values)

    unique_ells = np.unique(data["MULTIPOLE-I"])
    ell_to_block_idx = {ell: i for i, ell in enumerate(unique_ells)}
    n_ells = len(unique_ells)

    total_size = n_ells * n_scale
    full_cov_matrix = np.zeros((total_size, total_size))

    for row in data:
        ell_i = row["MULTIPOLE-I"]
        ell_j = row["MULTIPOLE-J"]
        scale_i = row[scale_i_col]
        scale_j = row[scale_j_col]
        cov = row["COVARIANCE"]

        block_row_start = ell_to_block_idx[ell_i] * n_scale
        block_col_start = ell_to_block_idx[ell_j] * n_scale

        scale_idx_i = scale_to_idx[scale_i]
        scale_idx_j = scale_to_idx[scale_j]

        final_row = block_row_start + scale_idx_i
        final_col = block_col_start + scale_idx_j
        full_cov_matrix[final_row, final_col] = cov

    # Ensure the full matrix is symmetric
    full_cov_matrix = (
        full_cov_matrix + full_cov_matrix.T - np.diag(np.diag(full_cov_matrix))
    )

    return all_unique_scale_values, full_cov_matrix


def build_2d_correlation(
    s_1d: NDArray[Any],
    mu_1d: NDArray[Any],
    correlation_1d: NDArray[Any],
) -> Tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    Reshapes a 1D correlation array into a 2D matrix based on s and mu values.
    """
    unique_s, s_indices = np.unique(s_1d, return_inverse=True)
    unique_mu, mu_indices = np.unique(mu_1d, return_inverse=True)

    n_s = len(unique_s)
    n_mu = len(unique_mu)

    correlation_2d = np.zeros((n_s, n_mu))
    correlation_2d[s_indices, mu_indices] = correlation_1d

    return unique_s, unique_mu, correlation_2d
