import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


class Cube:
    def __init__(self, data, wvl_axis, wcs):
        self.data = data
        self.wvl_axis = wvl_axis
        self.wcs = wcs

    @classmethod
    def read(cls, filename, extname='SCI'):
        """Read in file to create cube object

        Args:
            filename (str): path to input file
        """
        # get the data and header
        data, header = fits.getdata(filename, extname=extname)

        # get wcs object
        wcs = WCS(header)

        # Convert data to Jy if necessary
        if header['BUNIT'] == "MJy/sr":
            pix_area = header["PIXAR_SR"]
            data = data * 1.0e6 * pix_area
            header['BUNIT'] = "Jy"

        # get wavelength axis
        start_wvl = header["CRVAL3"]
        wvl_step = header["CDELT3"]
        nchan = header["NAXIS3"]
        wvl_axis = (np.arange(nchan) * wvl_step) + start_wvl

        return cls(data=data, wvl_axis=wvl_axis, wcs=wcs)




