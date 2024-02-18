import xml.etree.ElementTree as ElT
import os
from ._table_reader import DpdLE3GCTable
from astropy.io import fits

class DpdLE3GCReader:

    def __init__(self):
        """
        Constructor of the class DpdLE3GCReader
        """

        self._workdir = None
        self._head_dict = None
        self._data_dict = None
        self._table_dict = None
        self._tables = {}

    @property
    def header_dictionary(self):
        return self._head_dict

    @property
    def data_dictionary(self):
        return self._data_dict

    @property
    def table_dictionary(self):
        return self._table_dict

    @property
    def header_keywords(self):
        return self._head_dict.keys()

    @property
    def data_keywords(self):
        return self._data_dict.keys()

    @property
    def table_keywords(self):
        return self._table_dict.keys()

    @property
    def table_names(self):
        return self.list_tables()

    @property
    def tables(self):
        return self._tables

    @property
    def workdir(self):
        return self._workdir

    def parse_xml(self, xml_file):
        """
        Parse the xml file in input, extracting the working directory,
        all the keywords and the paths of the tables, in format of a fits file
        :param xml_file: full path to the xml file
        """

        # Set the working directory
        self._workdir = os.path.dirname(os.path.abspath(xml_file)) + "/"

        # Parse the XML file
        root = ElT.parse(xml_file).getroot()

        # Read 'Headers' keywords
        self._head_dict = {}
        for child in root.find("Header"):
            self._head_dict[child.tag] = root.find("Header").find(child.tag).text

            # Assign Attribute for each key
            setattr(self, child.tag, self._head_dict[child.tag])

        # Read 'Data' Keywords
        self._data_dict = {}
        self._table_dict = {}

        for child in root.find("Data"):
            if "File" in child.tag:
                # Extract the file name
                self._table_dict[child.tag] = root.find("Data").find(child.tag).find("DataContainer").find(
                    "FileName").text
                setattr(self, child.tag, self._table_dict[child.tag])
            else:
                self._data_dict[child.tag] = root.find("Data").find(child.tag).text

                # Assign Attribute for each key
                setattr(self, child.tag, self._data_dict[child.tag])

    def check_table_names(self):
        """
        Check if the table names have been parsed from xml file
        """

        if self._table_dict is None:
            raise Exception("No table name found, please run parse_xml() first or check your .xml file.")

    def list_tables(self):
        """
        List the tables
        """
        self.check_table_names()

        return [tt.replace("File", "") for tt in self.table_keywords]

    def read_tables(self):
        """
        Read the tables from fits file into pandas DataFrames
        """
        for tt in self.table_keywords:
            hdulist = fits.open(self._workdir + self._table_dict[tt])
            n_hdu = len(hdulist)-1
            for i in range(n_hdu):
                hdu = i+1
                ext_name = hdulist[hdu].header["EXTNAME"]
                nn = tt.replace("File", "")+f"_{ext_name}"
                self._tables[nn] = DpdLE3GCTable(self._workdir + self._table_dict[tt], hdu)
                setattr(self, nn, self._tables[nn])

    def plot(self, table_name, x, y, **kwargs):
        """
        Plot the selected table using the chosen column
        If x is None, plot against the table index
        :param table_name: name of the table
        :type table_name: ``str``
        :param x: name of the variable to use as x
        :type x: ``str``
        :param y: name of the variable to use as y
        :type y: ``str``
        (see https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.plot.html)
        :return matplotlib.pyplot.axes object
        """
        if x is None:
            return getattr(self, table_name).table.reset_index().plot(kind="scatter", x="index", y=y, **kwargs)

        return getattr(self, table_name).table.plot(kind="scatter", x=x, y=y, **kwargs)

    def header(self, table_name):
        """
        Get the header of the FITS file
        """
        return getattr(self, table_name).header

