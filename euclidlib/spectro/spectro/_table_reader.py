from astropy.io import fits
from astropy.table import Table


class DpdLE3GCTable:

    def __init__(self, fits_file, hdu):
        """
        Constructor of DpdLE3GCTable, that reads fits file, converts table to pandas and saves header
        :param fits_file: full path of the fits file
        :param hdu: hdu number
        """

        self.table = None
        self.header = None

        self.read(fits_file, hdu)

    def read(self, fits_file, hdu):
        """
        Read fits file, convert table to pandas and save header
        :param fits_file: full path of the fits file
        :param hdu: hdu number
        """
        self.table = Table.read(fits_file, hdu=hdu).to_pandas()
        self.header = fits.open(fits_file)[hdu].header
