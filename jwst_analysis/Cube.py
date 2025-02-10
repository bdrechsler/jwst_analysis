import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.signal import correlate2d

class Cube:
    def __init__(self, data, header, wvl_axis, wcs):

        self.data = data
        self.header = header
        self.wvl_axis = wvl_axis
        self.wcs = wcs

    def combine_cubes(self, new_cube):
        """Align and combine the cube with another cube

        Args:
            new_cube (Cube): Second cube object

        Returns:
            combined_cube (cube): combined cube object
        """
        # take median to get representative image
        img1 = np.nan_to_num(np.nanmedian(self.data, axis=0))
        img2 = np.nan_to_num(np.nanmedian(new_cube.data, axis=0))

        # take cross correlation to find offset
        corr = correlate2d(img1, img2, mode='full', boundary='fill',
                           fillvalue=0)
        # find the peak of the correlation
        y_max, x_max = np.unravel_index(np.argmax(corr), corr.shape)
        # calculate the offset
        y_offset = (img1.shape[0] - 1) - y_max
        x_offset = (img1.shape[1]) - x_max

        # align two images and use entire FOV
        # get dimensions of original image
        ny, nx = img1.shape
        # get dimensions of final image
        nx_big = np.abs(x_offset) + nx
        ny_big = np.abs(y_offset) + ny

        # determine where to place the new/ old arrays based on the direction of
        # the offset, also how much to shift the reference pixels
        if x_offset < 0 and y_offset < 0:
            old_x_bounds = (0, nx)
            old_y_bounds = (0, ny)
            new_x_bounds = (np.abs(x_offset), nx_big)
            new_y_bounds = (np.abs(y_offset), ny_big)
            x_shift, y_shift = -x_offset, -y_offset
        elif x_offset < 0 and y_offset > 0:
            old_x_bounds = (0, nx)
            old_y_bounds = (np.abs(y_offset), ny_big)
            new_x_bounds = (np.abs(x_offset), nx_big)
            new_y_bounds = (0, ny)
            x_shift, y_shift = -x_offset, 0
        elif x_offset > 0 and y_offset < 0:
            old_x_bounds = (np.abs(x_offset), nx_big)
            old_y_bounds = (0, ny)
            new_x_bounds = (0, nx)
            new_y_bounds = (np.abs(y_offset), ny_big)
            x_shift, y_shift = 0, -y_offset
        elif x_offset > 0 and y_offset > 0:
            old_x_bounds = (np.abs(x_offset), nx_big)
            old_y_bounds = (np.abs(y_offset), ny_big)
            new_x_bounds = (0, nx)
            new_y_bounds = (0, ny)
            x_shift, y_shift = 0, 0

        # get bounds of overlap region
        x1 = np.abs(x_offset) + 1
        x2 = -np.abs(x_offset) - 1
        y1 = np.abs(y_offset) + 1
        y2 = -np.abs(y_offset) - 1

        # create footprint for the combined cube
        combined_cube = np.zeros((self.data.shape[0], ny_big, nx_big))

        # shift each channel
        for i in range(self.data.shape[0]):
            # get new and old channel
            old_chan = np.nan_to_num(self.data[i])
            new_chan = np.nan_to_num(new_cube.data[i])

            big_chan = np.zeros((ny_big, nx_big))
            # add the old channel
            big_chan[old_y_bounds[0]:old_y_bounds[1], old_x_bounds[0]:
                    old_x_bounds[1]] += old_chan
            # add the new channel
            big_chan[new_y_bounds[0]: new_y_bounds[1], new_x_bounds[0]:
                    new_x_bounds[1]] += new_chan
            # divide overlap region by 2 to not double count (average)
            big_chan[y1:y2, x1:x2] /= 2
            # set 0s to nans
            big_chan[big_chan == 0] = np.nan
            combined_cube[i] = big_chan

        # update the header
        new_header = new_cube.header
        combined_header = new_header.copy()
        # get reference pixel position
        ref_1 = new_header['CRPIX1']
        ref_2 = new_header['CRPIX2']
        # set new dimensions
        combined_header['NAXIS1'] = ny_big
        combined_header['NAXIS2'] = nx_big
        # set new reference pixel position
        combined_header['CRPIX1'] = ref_1 + x_shift
        combined_header['CRPIX2'] = ref_2 + y_shift

        combined_wcs = WCS(combined_header)

        # create a return a new cube object
        return Cube(data=combined_cube, header=combined_header,
                    wvl_axis=new_cube.wvl_axis, wcs=combined_wcs)

    @classmethod
    def read(cls, filename, extname='SCI'):
        """Read in file to create cube object

        Args:
            filename (str): path to input file
        """
        # get the data and header
        data, header = fits.getdata(filename, extname=extname, header=True)

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

        return cls(data=data, header= header, wvl_axis=wvl_axis, wcs=wcs)

    def write(self, filename):
        hdu = fits.PrimaryHDU(data=self.data, header=self.header)
        hdu.writeto(filename, overwrite=True)








