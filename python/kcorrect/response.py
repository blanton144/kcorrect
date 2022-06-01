#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: filter.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import os
import numpy as np
import kcorrect.utils
import scipy.interpolate as interpolate
import scipy.integrate as integrate
import pydl.pydlutils.yanny as yanny
import fitsio


# Class to define a singleton
class FilterDictSingleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(FilterDictSingleton,
                                        cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class FilterDict(dict, metaclass=FilterDictSingleton):
    """Dictionary of all filters (singleton)
"""
    def __init__(self):
        return

    def load_filter(self, filter=None, reload=False):
        """Load filter into dictionary

        Parameters
        ----------

        filter : str
            filter to load

        reload : bool
            if True, reload the filter if already in FilterDict (default False)
"""
        if((filter in self) & (reload is False)):
            return
        self[filter] = Filter()
        self[filter].frompar(filename='{n}.par'.format(n=filter))
        return


class Filter(object):
    """Astronomical bandpass description

    Parameters
    ----------

    wave : ndarray of np.float32
        wavelength grid in Angstroms

    response : ndarray of np.float32
        response function

    Notes
    -----

    The response should be the relative response of the system
    (atmosphere, telescope, detector, etc) to a photon at each 
    given wavelength entering the Earth's atmosphere (for a ground
    based telescope) or the telescope aperture (space based).
"""
    def __init__(self, filename=None, wave=None, response=None):
        self.wave = wave
        self.response = response
        self.filename = filename
        if(self.filename is not None):
            self.frompar(filename)
        else:
            if(self.wave is not None):
                self.nwave = len(self.wave)
                self.set_interp()
        return

    def response_dtype(self):
        """Returns numpy dtype for SED"""
        response_dtype = np.dtype([('wave', np.zeros(self.nwave, dtype=np.float32)),
                                   ('response', np.zeros(self.nwave, dtype=np.float32))])
        return(response_dtype)

    def set_interp(self):
        """Sets attribute interp to interpolation function"""
        if((self.wave is None) | (self.response is None)):
            self.interp = None
            return
        self.interp = interpolate.interp1d(self.wave, self.response,
                                           kind='cubic',
                                           bounds_error=False,
                                           fill_value=0.)
        return

    def fromfits(self, filename=None, ext='RESPONSE'):
        """Read response from FITS files

        Parameters
        ----------

        filename : str
            input file name

        ext : str or int
            extension to read from
"""
        filter = fitsio.read(filename, ext=ext)
        self.nwave = len(filter)
        self.wave = filter['wave']
        self.response = filter['response']
        self.set_interp()
        return

    def frompar(self, filename=None):
        """Read response from Yanny parameter file

        Parameters
        ----------

        filename : str
            input file name

        Notes
        -----

        If an absolute path, reads that. If not, looks relative
        to $KCORRECT_DIR/python/kcorrect/data/filters
"""
        if(os.path.isabs(filename)):
            infilename = filename
        else:
            infilename = os.path.join(os.getenv('KCORRECT_DIR'),
                                      'python', 'kcorrect', 'data',
                                      'filters', filename)

        par = yanny.yanny(infilename)
        name = par.tables()[0]
        parstr = par[name]
        self.nwave = len(parstr)
        self.wave = parstr['lambda']
        self.response = parstr['pass']
        self.set_interp()
        return

    def tofits(self, filename=None, ext='FLUX', clobber=True):
        """Write SED to FITS files

        filename : str
            output file name

        ext : str or int
            extension to write to
        
        clobber : bool
            whether to clobber the existing file or add an HDU
"""
        out = np.zeros(1, self.response_dtype())
        out['wave'] = self.wave
        out['response'] = self.flux

        fitsio.write(filename, out, clobber=clobber)
        return

    def project(self, sed=None, wave=None, flux=None):
        """Project spectrum onto filter

        Parameters
        ----------

        sed : kcorrect.SED object
            spectrum in kcorrect format

        wave : ndarray of np.float32
            wavelength grid (used only if sed and func not set)

        flux : ndarray of np.float32
            flux grid (used only if sed and func not set)

        Returns
        -------

        maggies : np.float32
            [nsed] nmaggies associated with spectrum through bandpass

        Notes
        -----

        Fluxes should be in erg cm^{-2} s^{-1} Angstrom^{-1}

        Assumes AB calibration.
"""
        if(sed is None):
            sed_wave = wave
            interp = interpolate.interp1d(wave, flux, kind='cubic',
                                          bounds_error=False, fill_value=0.)
            nsed = 1
        else:
            sed_wave = sed.wave
            interp = sed.interp
            nsed = sed.nsed

        # Find SED wavelengths to integrate over
        keep = (sed_wave >= self.wave[0]) & (sed_wave <= self.wave[-1])
        ikeep = np.where(keep)[0]
        if(ikeep[0] > 0):
            keep[ikeep[0] - 1] = 1
        if(ikeep[-1] < len(sed_wave) - 1):
            keep[ikeep[-1] + 1] = 1

        # Find full grid of wavelengths for integration
        integrate_wave = np.unique(np.append(sed_wave[keep], self.wave))

        # Interpolate to grid
        integrate_sed = interp(integrate_wave)
        integrate_response = self.interp(integrate_wave)

        # Perform integration for numerator
        numer = np.zeros(nsed, dtype=np.float32)
        if(nsed == 1):
            integrand_numer = integrate_sed * integrate_response * integrate_wave
            numer = integrate.trapezoid(integrate_wave, integrand_numer)
        else:
            for ised in np.arange(nsed, dtype=int):
                integrand_numer = integrate_sed[ised, :] * integrate_response * integrate_wave
                numer[ised] = integrate.trapezoid(integrate_wave,
                                                  integrand_numer)

        # Perform integration for denominator
        integrand_denom = (kcorrect.utils.sed_ab(integrate_wave) *
                           integrate_response * integrate_wave)
        denom = integrate.trapezoid(integrate_wave, integrand_denom)

        # AB maggies are projection of SED onto filter divided by same
        # projection for the AB source.
        maggies = np.squeeze(numer / denom)

        return(maggies)
