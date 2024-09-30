"""

:Name: photo_data.py

:Description: This script contains classes that read and
process LE3 output files for photometric 2D statistics.

:Authors: Martin Kilbinger, Jingwei Wang

:Date: 2021

"""

import os
import glob

import re
import numpy as np
from astropy.io import fits
from astropy import units


def fit_bins(ells, es, nbins=None, bin=None):
    """Fit Bins.

    Rebin the ell-Cell data.

    :param ells: ells of the power spectrum
    :param es: Cells power spectrum
    :param nbins: defaut as 30
    :param bin: if given, the width of bins is set to "bin", no matter how much is the "nbins"
    :return: [ells, C_ells] after rebin

    """
    ma = max(ells)
    if bin is None:
        if nbins is None:
            nbins = 30
            wid = round(ma / nbins)  # the width of bins
        else:
            ma = max(ells)
            wid = round(ma / nbins)  # the width of bins
    else:
        wid = bin
    ourbin = np.arange(0, ma, wid)
    counts = np.zeros(len(ourbin))
    values = np.zeros(len(ourbin))
    for i in range(len(ells)):
        order = int(ells[i] // wid)
        values[order] = (values[order] * counts[order] + es[i]) / (counts[order] + 1)
        counts[order] += 1
    return ourbin[values > 0] + round(wid / 2), values[values > 0]


def parse_file_name(file_name):
    """Parse File Name.

    Parse and match input file name and return information split up into
    indivdual parts. See https://euclid.roe.ac.uk/dmsf/files/1695/view .

    Parameters
    ----------
    file_name : str
        input file name

    Raises
    ------
    ValueError :
        if match fails
    IndexError :
        if number of matches is incorrect

    Returns
    -------
    dict
        information fields in file name

    """
    fname_info = {}

    # Number of expected matching groups
    n_match = 4

    # Numbero of matches for time stamp
    n_match_time_stamp = 2

    pattern = "EUC_(\w{3})_(\S+)_(\S{18}).(\d{2}\.\d{2})\.fits"
    res_match = re.match(pattern, file_name)
    if res_match:

        if len(res_match.groups()) == n_match:
            fname_info["pf"] = res_match.group(1)
            fname_info["dat_prod"] = res_match.group(2)
            fname_info["time_stamp"] = res_match.group(3)
            fname_info["release"] = res_match.group(4)

            # Split up date
            pattern_time_stamp = "(\d{8})T(\d{6}\.\d{1})Z"
            res_match_time_stamp = re.match(
                pattern_time_stamp, fname_info["time_stamp"]
            )
            if res_match_time_stamp:
                if len(res_match_time_stamp.groups()) == n_match_time_stamp:
                    fname_info["date"] = res_match_time_stamp.group(1)
                    fname_info["time"] = res_match_time_stamp.group(2)
                else:
                    raise IndexError(
                        "Invalid number of matches in time stamp,"
                        + f" {len(res_match_time_stamp.groups())} instead of"
                        + f" {n_match_time_stamp}"
                    )
            else:
                raise ValueError("No match for time_stamp part in file name")

        else:
            raise IndexError(
                "Invalid number of matches in file name",
                +f" {len(res_match.groups())} instead of" + f" {n_match}",
            )
    # else:
    # raise ValueError('No match for file name')

    return fname_info


class PhotoDataError(Exception):
    """Photo Data Error

    Generic exception
    """

    pass


class PhotoDataReader(object):
    """Photo Data Reader

    Basic class to read LE3 output products

    Parameters
    ----------
    corr_str : str, optional, default='ShearShear'
        correlation type, one in 'ShearShear, 'PosShear', 'PosPos'

    Raises
    ------
    PhotoDataError
        for invalid correlation string

    """

    # relative precision for consistency checks
    prec_rel = 1e-4

    def __init__(self, corr_str="ShearShear"):

        # Real- or Fourier space
        self._space = None

        # Correlation type
        self._corr_str = corr_str
        if corr_str not in ["ShearShear", "PosShear", "PosPos"]:
            raise PhotoDataError(f"Invalid correlation string {corr_str}")

        # Estimator
        self._estim = None

        # Correlation data
        self._data = {}
        self._data[corr_str] = {}

        # Angular scales
        self._scales = {}

        # Unit of angular scales, e.g. 'arcmin'
        self._unit = None

    def __str__(self):
        s = ""
        s = f"{s}\nspace = {self._space}"
        s = f"{s}\nestimator = {self._estim}"
        s = f"{s}\ncorr_str = {self._corr_str}"
        for zcomb in self._data[self._corr_str]:
            s = f"{s}\n{zcomb}: {self._data[self._corr_str][zcomb]}"
        s = f"{s}\nscales = {self._scales}"
        s = f"{s}\nunit = {self._unit}"
        if self._space == "real":
            s = f"{s}\nonly_xi_p = {self._only_xi_p}"

        return s

    def init_from_fits(self, file_path, ignore_thunit=False):
        """Init From FITS.

        Initialise class instance from FITS input file information, without
        copying content. Useful if some (meta-)data is required
        before reading the file.

        Parameters
        ----------
        file_path : str
            input FITS file path
        ignore_thunit : bool, optional
            if False (default), ignore THUNIT in fits header to set
            units of angular scales.

        Raises
        ------
        FileNotFoundError
            if input file is not found
        KeyError
            if "THUNIT" key not found in header, with ignore_thunit=True

        Returns
        -------
        list
            list of HDUs

        """
        try:
            hdu_list = fits.open(file_path)
        except FileNotFoundError:
            raise (f"Input file {file_path} not found")

        # Get units of angular scales
        if not ignore_thunit:
            if "THUNIT" in hdu_list[0].header:
                unit_str = hdu_list[0].header["THUNIT"]
                self._unit = units.Unit(unit_str)
            else:
                raise KeyError(
                    f"Key THUNIT not found in header of input file {file_path}"
                )
        else:
            unit_str = None

        return hdu_list

    def check_consistency_scales(self, x, key):
        """Check Consistency Scales

        Check consistency of input and class data, of the scales column,
        regarding data vector lengths and element-wise differences.
        Use to check consistency between new and already read HDUs.

        Paramters
        ---------
        x : array(float)
            input data
        key : string
            name of scales key

        Raises
        ------
        PhotoDataError
            if data are not consistent

        """
        y = self._scales[key]

        # Check length
        if len(x) != len(y):
            raise PhotoDataError(
                f"column '{key}': inconsistent lengths " + f"{len(x)}, {len(y)}"
            )

        # Check values of elements
        for xx, yy in zip(x, y):
            if xx != 0:
                f = yy / xx - 1
            else:
                f = yy - xx

            if np.abs(f) > self.prec_rel:
                raise PhotoDataError(f"values {xx}, {yy} are inconsistent")

    @classmethod
    def check_redshifts(cls, z1, z2, corr_str="ShearShear"):
        """Check Redshifts.

        Check consistency of redshift indices.

        Parameters
        ----------
        z1, z2 : int
            redshift bin indices
        corr_str : string, optional, default='ShearShear'

        Returns
        -------
        res : bool
            True (False) if redshifts are valid (invalid)

        """
        # Redshifts cannot be negative
        if z1 < 0 or z2 < 0:
            return False

        # Second redshift cannot be larger than first except
        # for cross-correlation
        if z1 > z2 and corr_str != "PosShear":
            return False

        return True

    @classmethod
    def get_zcomb_str(cls, z1, z2):
        """Get Zcomb Str

        Return string representation redshift combination from two redshifts

        Parameters
        ----------
        z1, z2 : int or str
            redshift bin indices

        Returns
        -------
        zcomb : str
            redshift combination

        """
        return f"{z1}-{z2}"

    def get_zcomb_str_from_hdu(self, header):
        """Get Zcomb Str From HDU.

        Extract redshift combination from FITS header.

        Parameters
        ----------
        hdu : astropy.io.fits.header.Header
            FITS header

        Raises
        ------
        ValueError
            if instance estimator is not valid
        ValueError
            if header field with redshift information does
            not match expected pattern
        PhotoDataError
            if redshift combination is not valid

        Returns
        -------
        str
            redshift combination string

        """
        if self._estim == "XI":
            pattern_estim = ""
        elif self._estim == "BANDPOWER":
            pattern_estim = "PEB_"
        elif self._estim == "COSEBI":
            pattern_estim = "COSEBI_"
        else:
            raise ValueError(f'Invalid estimator "{self._estim}"')

        # Get redshift bins from bin combination
        pattern = f"{self._corr_str.upper()}2D_{pattern_estim}(\d+)_(\d+)"
        ext_name = header["EXTNAME"]
        match_res = re.match(pattern, ext_name)
        if not match_res or len(match_res.groups()) != 2:
            raise ValueError(
                f'HDU extension name "{ext_name}", did not match expected '
                + f'pattern "{pattern}"'
            )

        z1 = int(match_res.group(1))
        z2 = int(match_res.group(2))

        if not self.check_redshifts(z1, z2):
            raise PhotoDataError(f"Invalid redshifts {z1}, {z2}")

        # Create redshift combination and dictionary
        zcomb = self.get_zcomb_str(z1, z2)

        return zcomb

    @classmethod
    def parse_data_product(cls, dat_prod_str):
        """Parse Data Product.

        Parse data product string and extract information.

        Parameters
        ----------
        dat_prod_str : str
            input string

        Raises
        ------
        IndexEerror
            if number of matches is different from expected number
        ValueError
            if no match

        Returns
        -------
        dict
            data product information

        """
        dat_prod = {}

        # Number of expected matches
        n_match = 4

        pattern = "(\S+)-WL-(\S+?)-(\S+?)-(\S+)"
        res_match = re.match(pattern, dat_prod_str)
        if res_match:
            if len(res_match.groups()) == n_match:

                # '2PCF' or 'PK'
                dat_prod["space_abbr"] = res_match.group(1)

                # 'CS', 'GA', or 'SP' [TBC]
                dat_prod["corr_abbr"] = res_match.group(2)

                # 'XI', 'BANDPOWER', 'COSEBI', with or without suffix 'XICOV'
                dat_prod["estim"] = res_match.group(3)

                dat_prod["else"] = res_match.group(4)
            else:
                raise IndexError(
                    "Invalid number of matches in data product string",
                    +f" {len(res_match.groups())} instead of" + f" {n_match}",
                )
        else:
            raise ValueError("No match for data product string")

        return dat_prod

    @classmethod
    def get_dir_content(cls, dir_name):
        """Get Dir Content.

        Return information about LE3-PF-WL output files in directory.

        Parameters
        ----------
        dir_name : str
            directory name

        Returns
        -------
        list
            file name array
            file-name information array
            data product information array

        """
        file_info_arr = []

        file_paths = glob.glob(f"{dir_name}/*.fits")

        # Loop over all FITS files in directory
        for fp in file_paths:
            fname = os.path.basename(fp)

            # Get general SGS information from file name
            fname_info = parse_file_name(fname)

            if fname_info:
                # Get LE3-WL-specific information from data product field
                dat_prod = cls.parse_data_product(fname_info["dat_prod"])
                # print(dat_prod)

                # Copy all data
                file_info = {}
                file_info["fname"] = fname
                for key in fname_info:
                    file_info[key] = fname_info[key]
                for key in dat_prod:
                    file_info[key] = dat_prod[key]

                file_info_arr.append(file_info)

        return file_info_arr


class TpcfDataReader(PhotoDataReader):
    """Tpcf Data Reader

    Read data for photometric 2D real-space two-point correlation functions
    and variates.

    """

    def __init__(self, corr_str="ShearShear"):
        super().__init__(corr_str=corr_str)
        self._space = "real"
        self._estim = "XI"

    @classmethod
    def from_fits(cls, file_path, corr_str="ShearShear"):
        """From Fits

        Read from FITS intput file.

        Parameters
        ----------
        file_path : string
            input FITS file path
        corr_str : str, optional, default='ShearShear'
            correlation type, one in 'ShearShear, 'PosShear', 'PosPos'

        Raises
        ------
        FileNotFoundError
            for missing input ITS file
        PhotoDataError
            for invalid redshift bins

        Returns
        -------
        TpcfDataReader
            instance of this class

        """
        # Initialize instance
        tp = cls(corr_str=corr_str)

        # Initialize info from FITS files
        hdu_list = tp.init_from_fits(file_path)

        # short-hand for the data vector dictionary
        data = tp._data[tp._corr_str]

        tp._only_xi_p = False

        # Loop over HDU's = z-bin combinations
        for index in range(1, len(hdu_list)):

            # Get redshift bins from bin combination
            zcomb = tp.get_zcomb_str_from_hdu(hdu_list[index].header)
            if zcomb not in data:
                data[zcomb] = {}

            # Get angular scales
            theta = hdu_list[index].data["THETA"]
            if index == 1:
                # Copy fo class instance
                tp._scales["theta"] = theta
            else:
                # Check consistency of scales
                tp.check_consistency_scales(theta, "theta")

            # Get correlation function data
            if corr_str == "ShearShear":
                xi_p = hdu_list[index].data["XI_P"]
                xi_m = hdu_list[index].data["XI_M"]
                xi_x = hdu_list[index].data["XI_X"]

                # 2PCF data vector is (xip, xim), technically this is not
                # a pure E-mode
                data[zcomb]["E-E"] = np.concatenate((xi_p, xi_m), axis=None)
                data[zcomb]["E-B"] = xi_x

            elif corr_str == "PosShear":
                gamma_t = hdu_list[index].data["GAMMA_T"]
                gamma_x = hdu_list[index].data["GAMMA_X"]

                data[zcomb]["E-E"] = gamma_t
                data[zcomb]["E-B"] = gamma_x

            elif corr_str == "PosPos":
                w = hdu_list[index].data["WTHETA"]

                data[zcomb]["E-E"] = w

        # Duplicate angular scale vector to match data vector
        # (xip, xim)
        if corr_str == "ShearShear":
            tp._scales["theta"] = np.concatenate((theta, theta))

        # For all correlation types, theta_uniq contains uniq input angular scales
        tp._scales["theta_uniq"] = theta

        return tp

    @classmethod
    def from_arrays(cls, input_data, zcomb_list, corr_str="ShearShear"):
        """From Arrays

        Read from numpy arrays.

        Parameters
        ----------
        input_data : array(nzomb, 2, float)
            input correlation data, nzcomb bins of (theta, xi)
        zcomb_list : array of string
            redshift combinations
        corr_str : bool, optional, default='ShearShear'
            one in 'ShearShear', 'ShearPos', 'PosPos'

        Raises
        ------
        FileNotFoundError
            for missing input ITS file
        PhotoDataError
            for invalid redshift bins
        """

        # Initialize instance
        tp = cls(corr_str=corr_str)

        tp._only_xi_p = False

        data = tp._data[tp._corr_str]

        # Create dictionary for data at each redshift combination
        for index, zcomb in enumerate(zcomb_list):
            if zcomb not in data:
                data[zcomb] = {}

            # Get angular scales
            theta = input_data[index][0]
            if index == 0:
                tp._scales["theta"] = theta
            else:
                # Check consistency of scales
                tp.check_consistency_scales(theta, "theta")

            # Get correlation
            if corr_str == "ShearShear":
                xi_p = input_data[index][1]
                if len(input_data[index]) > 2:
                    xi_m = input_data[index][2]
                    data[zcomb]["E-E"] = np.concatenate((xi_p, xi_m), axis=None)
                else:
                    data[zcomb]["E-E"] = xi_p
                    tp._only_xi_p = True
                if len(input_data[index]) == 4:
                    xi_x = input_data[index][3]
                    data[zcomb]["E-B"] = xi_x
            else:
                raise PhotoDataError("Only ShearShear is implemented as yet")

        # Duplicate angular scale vector to match data vector
        # (xip, xim)
        if corr_str == "ShearShear" and not tp._only_xi_p:
            data._scales["theta"] = np.concatenate((theta, theta))
            data._scales["theta_uniq"] = theta

        return tp

    def get_xi_component(self, component, zcomb):
        """get xi component
        Return one of the three components of the shear-shear
        correlation function, xi+, xi_-, xi_x.

        Parameters
        ----------
        component : str
            'p', '+' for xi_p
            'm', '/'-' for xi_m
            'x' for xi_x
        zcomb : string
            redshift combination

        Returns
        -------
        xi_comp : array of float
            values of the xi component

        Raises
        ------
        ValueError :
            for invalid component on input
        """

        # Number of unique angular scales = divider in
        # xi array between components
        num_scales = len(self._scales["theta_uniq"])

        # Get correlation data
        key = "ShearShear"
        if key not in self._data:
            raise KeyError(f"{key} correlation not found in data")
        data = self._data[key]
        data = self._data[key]

        if component in ("p", "+"):
            # First half of E-E array
            xi_comp = data[zcomb]["E-E"][:num_scales]
        elif component in ("m", "-"):
            # Second half of E-E array
            xi_comp = data[zcomb]["E-E"][num_scales:]
        elif component == "x":
            # Parity E-B array
            xi_comp = data[zcomb]["E-B"]
        else:
            raise ValueError("Invalid component 'f{component}'")

        return xi_comp

    def get_gamma_component(self, component, zcomb):
        """get gamma component
        Return one of the two components of the position-shear
        correlation function <gamma_t>, <gamma_x>

        Parameters
        ----------
        component : str
            't' for <gamma_t>
            'x' for <gamma_x>
        zcomb : string
            redshift combination

        Returns
        -------
        gamma_comp : array of float
            values of the gamma component

        Raises
        ------
        ValueError :
            for invalid component on input
        KeyError :
            if 'PosShear' is not in data
        """

        # Get correlation data
        key = "PosShear"
        if key not in self._data:
            raise KeyError(f"{key} correlation not found in data")
        data = self._data[key]

        if component == "t":
            gamma_comp = data[zcomb]["E-E"]
        elif component == "x":
            gamma_comp = data[zcomb]["E-B"]
        else:
            raise ValueError("Invalid component 'f{component}'")

        return gamma_comp

    def get_w(self, zcomb):
        """get w
        Return the position-position
        correlation function w

        Parameters
        ----------
        zcomb : string
            redshift combination

        Returns
        -------
        w : array of float
            values of w

        Raises
        ------
        KeyError :
            if 'PosPos' is not in data
        """

        # Get correlation data
        key = "PosPos"
        if key not in self._data:
            raise KeyError(f"{key} correlation not found in data")
        data = self._data[key]

        w = data[zcomb]["E-E"]

        return w


class CellDataReader(PhotoDataReader):
    """Cell Data Reader

    Read data for photometric 2D angular power spectrum functions and variates.

    """

    def __init__(self, corr_str="ShearShear"):
        super().__init__(corr_str=corr_str)
        self._space = "fourier"
        self._estim = "BANDPOWER"

    @classmethod
    def from_fits(cls, file_path, corr_str="ShearShear"):
        """From Fits

        Read from FITS intput file, created by the LE3 2PCF-WL
        processing function.

        Parameters
        ----------
        file_path : string
            input FITS file path
        corr_str : str, optional, default='ShearShear'
            correlation type, one in 'ShearShear, 'PosShear', 'PosPos'

        Raises
        ------
        FileNotFoundError
            for missing input ITS file
        PhotoDataError
            for invalid redshift bins

        """
        # Initialize instance
        tp = cls(corr_str=corr_str)

        # Initialize info from FITS files
        hdu_list = tp.init_from_fits(file_path)

        # Short-cut to instance data
        data = tp._data[tp._corr_str]

        # Loop over HDU's = z-bin combinations
        for index in range(1, len(hdu_list)):

            # Short-cut to input data
            data_inp = hdu_list[index].data

            # Get redshift bins from bin combination
            zcomb = tp.get_zcomb_str_from_hdu(hdu_list[index].header)
            if zcomb not in data:
                data[zcomb] = {}

            # Get ell modes
            ell = data_inp["ELL"]
            ell_lower = data_inp["ELL_LOWER"]
            ell_upper = data_inp["ELL_UPPER"]
            if index == 1:
                # Copy fo class instance
                tp._scales["ell"] = ell
                tp._scales["ell_lower"] = ell_lower
                tp._scales["ell_upper"] = ell_upper
            else:
                # Check consistency of modes
                tp.check_consistency_scales(ell, "ell")
                tp.check_consistency_scales(ell_lower, "ell_lower")
                tp.check_consistency_scales(ell_upper, "ell_upper")

            # Get correlation
            if corr_str == "ShearShear":
                data[zcomb]["E-E"] = data_inp["PE"]
                data[zcomb]["B-B"] = data_inp["PB"]
                data[zcomb]["E-B"] = data_inp["PEB"]

            elif corr_str == "PosShear":
                data[zcomb]["E"] = data_inp["PE"]
                data[zcomb]["B"] = data_inp["PB"]

            else:
                raise ValueError(f"Correlation {corr_str} not yet implemented")

        return tp

    @classmethod
    def from_fits_camb(cls, file_path, corr_str="ShearShear"):
        """From Fits CAMB

        Read from FITS intput file, created by CAMB.
        processing function.

        Parameters
        ----------
        file_path : string
            input FITS file path
        corr_str : str, optional, default='ShearShear'
            correlation type, one in 'ShearShear, 'PosShear', 'PosPos'

        Raises
        ------
        FileNotFoundError
            for missing input ITS file
        PhotoDataError
            for invalid redshift bins

        """
        # Initialize instance
        tp = cls(corr_str=corr_str)

        # Initialize info from FITS files; ignore missing THUNIT
        hdu_list = tp.init_from_fits(file_path, ignore_thunit=True)
        hdu_num = 1

        # Short-cut to instance data
        data = tp._data[tp._corr_str]

        # Loop over data columns
        keys = hdu_list[hdu_num].data.dtype.names
        for key in keys:
            data_inp = hdu_list[hdu_num].data

            if key == "ell":
                # Get ell modes
                ell = data_inp["ell"]
                tp._scales["ell"] = ell
            else:
                # key is the redshift combination
                zcomb = key

                if zcomb not in data:
                    data[zcomb] = {}

                # Get correlation
                if corr_str == "ShearShear":
                    data[zcomb]["E-E"] = data_inp[zcomb]
                else:
                    raise ValueError(f"Correlation {corr_str} not yet implemented")

        return tp

    @classmethod
    def from_arrays(cls, input_data, zcomb_list, corr_str="ShearShear"):
        """From Arrays.

        Read from numpy arrays, created by the make glass jar library.

        Parameters
        ----------
        input_data : array(nzomb, 2, float)
            input power spectrum data, nzcomb bins of (ell, C_EE),
            or (ell, C_EE, C_BB, C_EB),
            or (ell, C_EE, C_BB, C_EB, Noise)
        zcomb_list : array of string
            redshift combinations
        corr_str : bool, optional, default='ShearShear'
            one in 'ShearShear', 'ShearPos', 'PosPos'

        Raises
        ------
        FileNotFoundError
            for missing input ITS file
        PhotoDataError
            for invalid redshift bins

        """
        # Initialize instance
        tp = cls(corr_str=corr_str)

        data = tp._data[tp._corr_str]

        # Create dictionary for data at each redshift combination
        for index, zcomb in enumerate(zcomb_list):

            # Make sure input data has correct number of entries
            if len(input_data[index]) not in (2, 4, 5):
                raise PhotoDataError(
                    f"Input data at index={index}, zcomb={zcomb} has length"
                    + f" {len(input_data[index])}, needs"
                    + " to be 2, 4, or 5"
                )

            if zcomb not in data:
                data[zcomb] = {}

            # Get angular scales
            ell = input_data[index][0]
            if index == 0:
                tp._scales["ell"] = ell
            else:
                # Check consistency of scales
                tp.check_consistency_scales(ell, "ell")

            # Get correlation
            if corr_str == "ShearShear":
                data[zcomb]["E-E"] = input_data[index][1]
                if len(input_data[index]) in (4, 5):
                    data[zcomb]["B-B"] = input_data[index][2]
                    data[zcomb]["E-B"] = input_data[index][3]

                if len(input_data[index]) == 5:
                    data[zcomb]["N"] = input_data[index][4]
            else:
                raise PhotoDataError("Only ShearShear is implemented as yet")

        return tp

    @classmethod
    def get_key(cls, corr_str, component):
        """Get Key.

        Return dictionary key of power-spectrum component.

        Parameters
        ----------
        corr_str : str
            correlation type
        component : str
            power-spectrum component

        Returns
        -------
        str
            key

        """
        if corr_str == "ShearShear":
            if component == "E":
                key = "E-E"
            elif component == "B":
                key = "B-B"
            elif component == "EB":
                key = "E-B"
        elif corr_str == "PosShear":
            key = component

        return key

    def get_cell_component(self, component, zcomb, dimensionless=False):
        """Get Cell Component.

        Return one of the three components of the shear-shear
        power spectrum, P_E, P_B, P_EB.

        Parameters
        ----------
        component : str
            'E' for P_E
            'B', for P_B
            'EB' for P_EB
        zcomb : string
            redshift combination
        dimensionless : bool, optional
            if ``True`` (default), return dimensionless power spectrum,
            i.e. ell^2 / (2 pi) * C_ell; default is ``False``

        Returns
        -------
        ndarray
            values of the cell component

        Raises
        ------
        ValueError :
            for invalid component on input
        ValueError :
            if key or correlation string not found in instance data

        """
        # Get correlation data
        if self._corr_str not in self._data:
            raise KeyError(f"{self._corr_str} correlation not found in data")

        # Check valid input component strings
        if component not in ("E", "B", "EB"):
            raise ValueError("Invalid component 'f{component}'")

        # Short-cut to data
        data = self._data[self._corr_str]

        key = self.get_key(self._corr_str, component)

        if key not in data[zcomb]:
            raise KeyError(f"Data for component 'f{component}' not found")

        cell_comp = data[zcomb][key]

        if dimensionless:
            pre_factor = self.get_pre_factor_dimensionless()
        else:
            pre_factor = 1

        return pre_factor * cell_comp

    def get_pre_factor_dimensionless(self):
        return self._scales["ell"] ** 2 / (2 * np.pi)


class COSEBIsDataReader(PhotoDataReader):
    """COSEBIs Data Reader

    Read data for photometric 2D cosmic shear COSEBIs.

    """

    def __init__(self):
        super().__init__(corr_str="ShearShear")
        self._space = "real"
        self._estim = "COSEBI"

    @classmethod
    def from_fits(cls, file_path):
        """From Fits.

        Read from FITS intput file.

        Parameters
        ----------
        file_path : string
            input FITS file path

        Raises
        ------
        FileNotFoundError
            for missing input ITS file
        PhotoDataError
            for invalid redshift bins

        """
        # The COSEBIs are only defined for shear-shear correlations
        # corr_str = "ShearShear"

        # Initialize instance
        tp = cls()

        # Initialize info from FITS files
        hdu_list = tp.init_from_fits(file_path)

        # short-hand for the data vector dictionary
        data = tp._data[tp._corr_str]

        # Loop over HDU's = z-bin combinations
        for index in range(1, len(hdu_list)):

            # Get redshift bins from bin combination
            zcomb = tp.get_zcomb_str_from_hdu(hdu_list[index].header)
            if zcomb not in data:
                data[zcomb] = {}

            # Get COSEBI mode numbers (integer)
            n_mode = hdu_list[index].data["MODE"]
            if index == 1:
                # Copy fo class instance
                tp._scales["n_mode"] = n_mode
            else:
                # Check consistency of modes
                tp.check_consistency_scales(n_mode, "n_mode")

            # Get correlation
            E_cosebis = hdu_list[index].data["EE"]
            B_cosebis = hdu_list[index].data["BB"]
            EB_cosebis = hdu_list[index].data["EB"]

            data[zcomb]["E-E"] = E_cosebis
            data[zcomb]["B-B"] = B_cosebis
            data[zcomb]["E-B"] = EB_cosebis

        return tp

    def get_cosebis_component(self, component, zcomb):
        """Get Cosebis Component.

        Return one of the three components of the (shear-shear)
        COSEBI modes E, B, EB.

        Parameters
        ----------
        component : str
            'E' for E COSEBI modes
            'B', for B COSEBI modes
            'EB' for EB COSEBI modes
        zcomb : string
        redshift combination

        Returns
        -------
        ndarray
            values of the COSEBIs component

        Raises
        ------
        ValueError :
            for invalid component on input

        """
        # Get correlation data
        key = "ShearShear"
        if key not in self._data:
            raise KeyError(f"{key} correlation not found in data")
        data = self._data[key]
        data = self._data[key]

        if component == "E":
            cosebis_comp = data[zcomb]["E-E"]
        elif component == "B":
            cosebis_comp = data[zcomb]["B-B"]
        elif component == "EB":
            cosebis_comp = data[zcomb]["E-B"]
        else:
            raise ValueError("Invalid component 'f{component}'")

        return cosebis_comp
