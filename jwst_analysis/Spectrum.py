import numpy as np
from scipy.ndimage import median_filter
from astropy.modeling import models, fitting

class Spectrum:
    def __init__(self, wvl_axis, flux):
        self.wvl_axis = wvl_axis
        self.flux = flux

    def __add__(self, other):
        if (self.wvl_axis != other.wvl_axis).any():
            print("Spectral axis are not matched, cannot add.")
            return
        else:
            summed_flux = self.flux + other.flux
            return Spectrum(wvl_axis=self.wvl_axis, flux=summed_flux)


    def fit_continuum(self, line=None, med_kernel=3, fit_order=3):
        """Fit the continuum of the spectrum, ignore line emission

        Args:
            line (Line): the Line to ignore. Defaults to None
            med_kernel (odd int): Size of kernel used in
            median smoothing. Defaults to 3.
            fit_order (int): Order of polynomial to fit
            to the continuum. Defaults to 3.

        Returns:
            continuum_spectrum (Spectrum): Spectrum of the fitted continuum
        """
        # define the bounds of the line
        line_min = line.rest_wvl - (line.line_width / 2.)
        line_max = line.rest_wvl + (line.line_width / 2.)

        # get the wvl and flux for the continuum region
        cont_wvl = self.wvl_axis[(self.wvl_axis < line_min) |
                                 (self.wvl_axis > line_max)]
        cont_flux = self.flux[(self.wvl_axis < line_min) |
                                 (self.wvl_axis > line_max)]

        # median smooth then fit a polynomial to the continuum spectrum
        cont_flux_smooth = median_filter(cont_flux, size=med_kernel)
        fit_params = np.polyfit(cont_wvl, cont_flux_smooth, fit_order)
        model_cont_flux = np.poly1d(fit_params)

        return Spectrum(wvl_axis=self.wvl_axis,
                        flux=model_cont_flux(self.wvl_axis))

    def fit_gaussian(self, line):
        # initial guess of parameters
        amp = np.max(self.flux)
        mean = line.rest_wvl
        stddev = line.line_width / 3.
        g_init = models.Gaussian1D(amplitude=amp, mean=mean, stddev=stddev)
        # fit the gaussian to the data
        fitter = fitting.TRFLSQFitter()
        g = fitter(g_init, self.wvl_axis, self.flux)
        # return the gaussian
        return g
