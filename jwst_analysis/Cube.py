import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clip
import astropy.units as u
from astropy.constants import c
from scipy.signal import correlate2d
import photutils.aperture as ap
from .Line import Line
from .Spectrum import Spectrum

class Cube:
    def __init__(self, data, header, wvl_axis):

        self.data = data
        self.header = header
        self.wvl_axis = wvl_axis
        self.wcs = WCS(header)

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

        # create a return a new cube object
        return Cube(data=combined_cube, header=combined_header,
                    wvl_axis=new_cube.wvl_axis)

    def attach_line(self, line):
        """Attach a line object to this spectral cube

        Args:
            line (Line): The line object to denote what line this cube
            corresponds to
        """
        self.line = line
        self.header['REST_WVL'] = line.rest_wvl
        self.header['SPECIES'] = line.species
        self.header['TRANSITION'] = line.transition
        self.header['LINE_WIDTH'] = line.line_width

    def extract_spectrum(self, aperture_dict):
        r"""
        Extract 1D spectrum from a dictionary of photutil apertures
        and populate the spectra and aperture attributes

        Args:
            aperture_dict (dict): Dictionary to define apertures used to extract spectra.
            In the format: {"name": aperture}, aperture can be a sky
            or pixel aperture

        Returns:
            spectrum_dict (dict): Dictionary of extracted spectra

        """

        # iterate through all the apertures in the dictionary
        for name, aperture in aperture_dict.items():
            # check if provided aperture is a pixel or sky aperture and convert it accordingly
            pixel_aps = {}
            sky_aps = {}
            if isinstance(aperture, ap.EllipticalAperture) or isinstance(
                aperture, ap.CircularAperture
            ):
                pixel_aps[name] = aperture
                sky_aps[name] = aperture.to_sky(self.wcs.celestial)
            elif isinstance(aperture, ap.SkyEllipticalAperture) or isinstance(
                aperture, ap.SkyCircularAperture
            ):
                pixel_aps[name] = aperture.to_pixel(self.wcs.celestial)
                sky_aps[name] = aperture

        spectrum_dict = {}
        # iterate through all of the pixel apertures
        for name, pix_ap in pixel_aps.items():
            # create a mask
            mask = pix_ap.to_mask(method="exact")

            # initialize spectrum array
            spectrum_flux = np.zeros(len(self.data))

            # extract the 1D spectrum
            for i in range(len(spectrum_flux)):
                # get data of current channel
                chan = self.data[i]
                # extract the data in the aperture
                ap_data = mask.get_values(chan)
                # sum to get value for spectrum
                spectrum_flux[i] = np.nansum(ap_data)

            # add to the dictionary of 1D spectra
            spectrum = Spectrum(wvl_axis=self.wvl_axis, flux=spectrum_flux)
            spectrum_dict[name] = spectrum

        # sum left and right apertures to get total spectrum
        spectrum_dict["total"] = spectrum_dict["left"] + spectrum_dict["right"]
        return spectrum_dict

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

        # update header with new wvl range
        header = self.header.copy()
        header['CRVAL3'] = wvl_min
        header['NAXIS3'] = len(region_wvl)

        return Cube(data=region_data, header=header, wvl_axis=region_wvl)

    def line_cube(self, line, n_lws):
        """Get a spectral region around a given line

        Args:
            line (Line): Line to center spectral region on
            n_lws (float): size of spectral region in units of line widths

        Returns:
            spectral region (Cube): spectral cube around the given line with the
            line object attached
        """
        # get region around the line
        line_region = self.spectral_region(center_wvl=line.rest_wvl,
                                           region_width=n_lws * line.line_width)
        # attach the line
        line_region.attach_line(line=line)
        return line_region


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
                         wvl_axis=self.wvl_axis)
        cont_sub_cube = Cube(data=cont_sub_data, header=self.header,
                         wvl_axis=self.wvl_axis)

        return cont_cube, cont_sub_cube

    def collapse(self, moments, sigma_thresh=3.):
        """Collapse the data cube to create moment maps

        Args:
            moments (list): the moments the method will return (0,1,8)

        Returns:
            moment_list: list of the requested moment maps
        """

        # check if a line has been attached
        if 'REST_WVL' not in list(self.header.keys()):
            print("Attatch a line to this spectral cube before making maps")
            return

        dv = self.vel_axis[1] - self.vel_axis[0] # velocity spacing
        M0 = np.nansum(self.data, axis=0) * dv
        # estimate noise on M0
        # use first and last two channels to estimate sigma for a
        # line free channel
        sigmas = [np.nanstd(self.data[i]) for i in range(-2, 2)]
        sigma_chan = np.mean(sigmas)
        # error propogation to get error on M0
        sigma_M0 = np.sqrt(len(self.data)) * sigma_chan * dv
        thresh = 3 * sigma_M0
        # create a mask for the M1 map
        mask = np.ones(M0.shape)
        mask[(M0 > -thresh) & (M0 < thresh)] = 0

        # next, calculate 1st moment map
        # create a velocity cube, same shape as data where each channel is just
        # filled with that channels velocity
        vel_cube = np.zeros(self.data.shape)
        for i in range(len(self.data)):
            vel_cube[i] = np.full(M0.shape, self.vel_axis[i])
        M1 = (np.nansum(vel_cube * self.data, axis=0) * dv) / M0
        M1 *= mask

        # last, get moment 8 map (peak intensity)
        M8 = np.nanmax(self.data, axis=0)
        # return the requested maps
        return_list = []
        if 0 in moments:
            return_list.append(M0)
        if 1 in moments:
            return_list.append(M1)
        if 8 in moments:
            return_list.append(M8)

        return return_list



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

        return cls(data=data, header= header, wvl_axis=wvl_axis)

    def write(self, filename):
        hdu = fits.PrimaryHDU(data=self.data, header=self.header)
        hdu.writeto(filename, overwrite=True)