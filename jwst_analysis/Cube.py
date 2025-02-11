import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.constants import c
from scipy.signal import correlate2d
import bettermoments as bm
from .Line import Line
from .Spectrum import Spectrum

class Cube:
    def __init__(self, data, header, wvl_axis, wcs):

        self.data = data
        self.header = header
        self.wvl_axis = wvl_axis
        self.wcs = wcs

        # check if a line has been attached to this cube
        header_keys = list(header.keys())
        if 'REST_WVL' in header_keys:
            self.line = Line(rest_wvl=header['REST_WVL'],
                             species=header['SPECIES'],
                             transition=header['TRANSITION'],
                             line_width=header['LINE_WIDTH'])

            # calculate the velocity axis
            wvl = np.array(wvl_axis) * u.um
            rest_wvl = self.line.rest_wvl * u.um
            vel = c * (wvl - rest_wvl) / rest_wvl
            self.vel_axis = vel.to(u.km/u.s).value

    @staticmethod
    def combine_cubes(old_cube, new_cube):
        """Align and combine the cube with another cube

        Args:
            new_cube (Cube): Second cube object

        Returns:
            combined_cube (cube): combined cube object
        """
        # take median to get representative image
        img1 = np.nan_to_num(np.nanmedian(old_cube.data, axis=0))
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
        combined_cube = np.zeros((old_cube.data.shape[0], ny_big, nx_big))

        # shift each channel
        for i in range(old_cube.data.shape[0]):
            # get new and old channel
            old_chan = np.nan_to_num(old_cube.data[i])
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

    def attach_line(self, line):
        """Attach a line object to this spectral cube

        Args:
            line (Line): The line object to denote what line this cube
            corresponds to
        """
        self.line = line
        self.header['REST_WVL'] = line.wvl
        self.header['SPECIES'] = line.species
        self.header['TRANSITION'] = line.transition
        self.header['LINE_WIDTH'] = line.line_width

    def spectral_region(self, center_wvl, region_width):
        """Get a spectral region (slab) of the spectral cube

        Args:
            center_wvl (float): central wavelength of the region [um]
            region_width (float): width of the spectral region [um]

        Returns:
            spectral_region (Cube): region of spectral cube within spectral
            window
        """
        wvl_min = center_wvl - (region_width / 2.)
        wvl_max = center_wvl + (region_width / 2.)

        region_data = self.data[(self.wvl_axis > wvl_min) &
                              (self.wvl_axis < wvl_max)]
        region_wvl = self.wvl_axis[(self.wvl_axis > wvl_min) &
                              (self.wvl_axis < wvl_max)]

        return Cube(data=region_data, header=self.header, wvl_axis=region_wvl,
                    wcs=self.wcs)

    def cont_sub(self):
        """Continuum subtract the spectral cube and create a continuum
        cube and continuum subtracted cube objects

        Returns:
            cont_cube (Cube): Continuum Cube
            cont_sub_cube (Cube): Contiuum subtracted Cube
        """
        # check if a line has been attached
        if 'REST_WVL' not in list(self.header.keys()):
            print("Attatch a line to this spectral cube before continuum \
                  subtracting.")
            return

        # create holder arrays for continuum and continuum
        # subtracted cubes
        cont_data = np.zeros(self.data.shape)
        cont_sub_data = np.zeros(self.data.shape)
        # perform the continuum subtraction for each pixel
        nchan, ni, nj = self.data.shape
        for i in range(ni):
            for j in range(nj):
                # fit the continuum to i, j spectrum
                pixel_spectrum = Spectrum(wvl_axis=self.wvl_axis,
                                          flux=self.data[:,i,j])
                continuum = pixel_spectrum.fit_continuum(self.line)
                # get the continuum subtractede flux
                cont_sub_flux = pixel_spectrum.flux - continuum.flux
                # set the i, j spectrum for the cont and cont sub cubes
                cont_data[:, i, j] = continuum.flux
                cont_sub_data[:, i, j] = cont_sub_flux

        # create continuum and continuum subtracted cube objects
        cont_cube = Cube(data=cont_data, header=self.header,
                         wvl_axis=self.wvl_axis, wcs=self.wcs)
        cont_sub_cube = Cube(data=cont_sub_data, header=self.header,
                         wvl_axis=self.wvl_axis, wcs=self.wcs)

        return cont_cube, cont_sub_cube

    def collapse(self, moment):
        """Collapse the spectral cube to create a moment map.
        There must be a line attached to the cube.

        Args:
            moment (int): which method should be used to collapse
            the data (0,1,8)

        Returns:
            M (np.ndarray): The resulting moment map
            dM (np.ndarray): Error on the map
        """

        # spectrally smooth the data
        smoothed_data = bm.smooth_data(data=self.data, smooth=3, polyorder=0)
        # estimate the noise
        rms = bm.estimate_RMS(data=smoothed_data, N=10)
        # create a threshold mask
        mask = bm.get_threshold_mask(data=smoothed_data, clip=2.0,
                                     smooth_threshold_mask=2)
        masked_data = smoothed_data * mask

        # collapse based on the desired moment map
        if moment == 0:
            M, dM = bm.collapse_zeroth(velax=self.vel_axis, data=masked_data,
                                       rms=rms)
        elif moment == 1:
            M, dM = bm.collapse_first(velax=self.vel_axis, data=masked_data,
                                       rms=rms)
        elif moment == 8:
            M, dM = bm.collapse_eighth(velax=self.vel_axis, data=masked_data,
                                       rms=rms)

        return M, dM



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








