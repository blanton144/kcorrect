#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: response.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import os
import re

import astropy.io.ascii
import astropy.io.fits
import astropy.table
import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import scipy.optimize as optimize

import kcorrect
import kcorrect.template
import kcorrect.utils


def all_responses(response_dir=os.path.join(kcorrect.KCORRECT_DIR, 'data',
                                            'responses'),
                  check_validity=False):
    """List all responses available
    
    Parameters
    ----------

    response_dir : str
        path to directory containing responses

    check_validity : bool
        if True, check the validity of the files

    Returns
    -------

    responses : list of str
        response names for available bandpasses
    
    Notes
    -----

    Returns all base names of files with the suffix '.dat' in the
    response_dir directory.

    By default, response_dir is the "data/responses" directory within
    the kcorrect Python distribution. 

    If the user specifies response_dir and check_validity is False,
    there is no guarantee that the response files in the specified
    directory are valid!

    If check_validity is True, the responses are also loaded into the
    ResponseDict() singleton.
"""
    rdir = os.path.join(response_dir)
    files = os.listdir(rdir)
    responses = []
    f = ResponseDict()
    for filename in files:
        if(os.path.isfile(os.path.join(rdir, filename))):
            m = re.match('^(.*)\\.dat$', filename)
            if(m is not None):
                response = m.group(1)
                valid = True
                if(check_validity):
                    try:
                        f.load_response(response)
                    except:
                        valid = False
                if(valid):
                    responses.append(response)
    return(responses)


# Class to define a singleton
class ResponseDictSingleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(ResponseDictSingleton,
                                        cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class ResponseDict(dict, metaclass=ResponseDictSingleton):
    """Dictionary of all responses (singleton)
"""
    def __init__(self):
        return

    def load_response(self, response=None, reload=False):
        """Load response into dictionary

        Parameters
        ----------

        response : str
            response to load

        reload : bool
            if True, reload the response if already in ResponseDict (default False)
"""
        if((response in self) & (reload is False)):
            return
        self[response] = Response()
        self[response].fromdat(filename='{n}.dat'.format(n=response))
        return


class Response(object):
    """Astronomical bandpass description

    Parameters
    ----------

    filename : str
        file name to read from

    wave : ndarray of np.float32
        wavelength grid in Angstroms

    response : ndarray of np.float32
        response function

    Attributes
    ----------

    filename : str
        source filename, or None

    fwhm, fwhm_low, fwhm_hight : np.float32
        FWHM of response, with low and high wavelength limits (Angstroms)

    interp : scipy.interpolate.interp1d object
        interpolation object

    lambda_eff : np.float32
        effective wavelength in Angstroms

    nwave : int
        number of wavelength samples

    response : ndarray of np.float32
        response function

    solar_magnitude : np.float32
        absolute magnitude of Sun through filter, or None

    solar_sed : kcorrect.template.SED object
        SED associated with Sun for solar_magnitude

    vega2ab : np.float32
        magnitude offset from Vega to AB (m_AB - m_Vega), or None

    vega_sed : kcorrect.template.SED object
        SED associated with Vega for vega_magnitude

    wave : ndarray of np.float32
        wavelength grid in Angstroms

    Notes
    -----

    The response should be the relative response of the system
    (atmosphere, telescope, detector, etc) to a photon at each 
    given wavelength entering the Earth's atmosphere (for a ground
    based telescope) or the telescope aperture (space based).

    If a file is given it is assumed to be in fixed_width 
    format, readable and writeable by astropy.io.ascii

    The wavelengths are sorted on input so the attribute wave is
    always increasing.

    The attribute interp() takes wavelength in Angstroms as its one
    positional argument.

    The effective wavelength is defined as described in Blanton &
    Roweis (2007).

    """
    def __init__(self, filename=None, wave=None, response=None):
        if((wave is not None) &
           (response is not None)):
            isort = np.argsort(wave)
            self.wave = wave[isort]
            self.response = response[isort]
        else:
            self.wave = wave
            self.response = response
        self.filename = filename
        self.solar_sed = None
        self.solar_magnitude = None
        self.vega_sed = None
        self.vega_magnitude = None
        self.lambda_eff = None
        self.interp = None
        self.fwhm = None
        self.fwhm_low = None
        self.fwhm_high = None
        if(self.filename is not None):
            self.fromdat(filename)
        else:
            if(self.wave is not None):
                self.nwave = len(self.wave)
                self._setup()
        return

    def _setup(self):
        """Some initial setup after an input"""
        self.set_interp()
        self.set_lambda_eff()
        self.set_solar_magnitude()
        self.set_vega2ab()
        self.set_fwhm_limits()
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
        response_hdulist = astropy.io.fits.open(filename)
        response = response_hdulist[ext]
        self.nwave = len(response)
        isort = np.argsort(response['wave'])
        self.wave = response['wave'][isort]
        self.response = response['response'][isort]
        self._setup()
        return

    def fromdat(self, filename=None):
        """Read response from fixed_width file

        Parameters
        ----------

        filename : str
            input file name

        Notes
        -----

        If an absolute path, reads that. If not, looks relative
        to KCORRECT_DIR/python/kcorrect/data/responses
"""
        if(os.path.isabs(filename)):
            infilename = filename
        else:
            infilename = os.path.join(kcorrect.KCORRECT_DIR, 'data',
                                      'responses', filename)

        if(os.path.exists(infilename) is False):
            raise ValueError("No response file: {f}".format(f=infilename))
        dat = astropy.io.ascii.read(infilename, format='fixed_width')
        self.nwave = len(dat)
        isort = np.argsort(dat['lambda'])
        self.wave = dat['lambda'][isort]
        self.response = dat['pass'][isort]
        self._setup()
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

        out_table = astropy.table.Table(out)

        out_table.writeto(filename, overwrite=clobber)
        return

    def project(self, sed=None, wave=None, flux=None):
        """Project spectrum onto response

        Parameters
        ----------

        sed : kcorrect.template.SED object
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

        If "flux" and "wave" are specified, then wave must be
        a 1-dimensional array, and flux must be a 1-dimensional
        or 2-dimensional array, with the last axis corresponding
        to wavelength and with the same number of elements as
        wave.

        Assumes AB calibration.

        If the bandpass is outside the range of the solar model, 0 is returned.
"""
        if(sed is None):
            if((wave is None) or (flux is None)):
                raise ValueError("must specify sed, or wave and flux")
            wave = np.float32(wave)
            flux = np.float32(flux)
            if(wave.ndim != 1):
                raise ValueError("wave must be 1-D array")
            if(flux.shape[-1] != wave.shape[0]):
                raise ValueError("last axis of flux must match wave")
            if(flux.ndim > 2):
                raise ValueError("flux must be 1-D or 2-D array")
            sed_wave = wave
            if(flux.ndim == 1):
                interp = interpolate.interp1d(wave, flux, kind='cubic',
                                              bounds_error=False, fill_value=0.)
                nsed = 1
            else:
                nsed = flux.shape[0]
                interp = interpolate.interp1d(wave, flux, kind='cubic',
                                              bounds_error=False, fill_value=0.)
        else:
            sed_wave = sed.wave
            interp = sed.interp
            nsed = sed.nsed

        # Find SED wavelengths to integrate over
        keep = (sed_wave >= self.wave[0]) & (sed_wave <= self.wave[-1])
        ikeep = np.where(keep)[0]
        if(len(ikeep) == 0):
            return(0.)

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

        # AB maggies are projection of SED onto response divided by same
        # projection for the AB source.
        maggies = np.squeeze(numer / denom)

        return(maggies)

    def set_lambda_eff(self):
        """Set effective wavelength

        Notes
        -----

        Sets attribute lambda_eff
"""
        # Just use original grid; good enough.
        wave = self.wave
        response = self.response

        # Perform integration for numerator
        integrand_numer = np.log(wave) * response / wave
        numer = integrate.trapezoid(wave, integrand_numer)

        # Perform integration for denominator
        integrand_denom = response / wave
        denom = integrate.trapezoid(wave, integrand_denom)

        # Set effective wavelength
        self.lambda_eff = np.exp(numer / denom)

        return

    def set_fwhm_limits(self):
        """Set limits for FWHM

        Notes
        -----

        Sets attributes fwhm_low, fwhm_high, fwhm.

        fwhm_low is the lowest wavelength value for which the response
        reaches 50% maximum when starting from the low end.

        fwhm_high is the highest wavelength value for which the response
        reaches 50% maximum when starting from the high end.

        fwhm is (fwhm_high - fwhm_low)
"""
        # Just use original grid; good enough.
        wave = self.wave.copy()
        iresponse = self.interp(wave)
        maxresponse = iresponse.max()
        iresponse = iresponse / maxresponse

        # Find lower
        iupper = np.where(iresponse > 0.5)[0][0]
        if(iupper == 0):
            fwhm_low = wave[iupper]
        else:
            ilower = iupper - 1
            fwhm_low = optimize.brentq(lambda x : (self.interp(x) / maxresponse - 0.5),
                                       wave[ilower], wave[iupper])

        # Find upper
        ilower = np.where(iresponse >= 0.5)[0][-1]
        if(ilower == len(iresponse) - 1):
            fwhm_high = wave[-1]
        else:
            iupper = ilower + 1
            fwhm_high = optimize.brentq(lambda x : (self.interp(x) / maxresponse - 0.5),
                                        wave[ilower], wave[iupper])

        self.fwhm_low = fwhm_low
        self.fwhm_high = fwhm_high
        self.fwhm = fwhm_high - fwhm_low
        return

    def set_solar_magnitude(self):
        """Set absolute magnitude of Sun through filter

        Notes
        -----

        Uses lcbsun.ori model from Lejeune et al (1997)

        If the response function is outside the model wavelength range,
        solar_magnitude is set to None.
"""
        if(self.solar_sed is None):
            sunfile = os.path.join(kcorrect.KCORRECT_DIR, 'data', 'basel',
                                   'lcbsun.ori')
            info, wave, flux = kcorrect.utils.read_basel(filename=sunfile)

            # Now convert to Angstroms and erg/cm^2/s/Ang at 10 pc
            radius = 6.960e+10
            wave = wave * 10.  # nm to Angstrom
            pctocm = 3.086e+18
            cspeed = 2.99792e+18   # Ang/s
            for unit in range(flux.shape[0]):
                flux[unit, :] = np.pi * 4. * flux * cspeed / wave**2
                flux = flux * (radius / (10. * pctocm))**2
            self.solar_sed = kcorrect.template.SED(wave=wave, flux=flux)
            self.solar_sed.info = info

        solar_maggies = self.project(sed=self.solar_sed)

        if(solar_maggies > 0):
            self.solar_magnitude = - 2.5 * np.log10(solar_maggies)
        else:
            self.solar_magnitude = None

        return

    def set_vega2ab(self):
        """Set Vega to AB magnitude conversion

        Notes
        -----

        Uses lcbvega.ori model from Lejeune et al (1997)

        If the response function is outside the model wavelength range,
        vega2ab is set to None.
"""
        if(self.vega_sed is None):
            vegafile = os.path.join(kcorrect.KCORRECT_DIR, 'data', 'basel',
                                    'lcbvega.ori')
            info, wave, flux = kcorrect.utils.read_basel(filename=vegafile)

            # Conversion to match Hayes et al. 1985
            radius = 1.91144e+11  # Backed out to get normalization right
            dvega = 7.68  # Vega is 7.68 pc
            wave = wave * 10.  # nm to Angstrom
            pctocm = 3.086e+18
            cspeed = 2.99792e+18   # Ang/s
            for unit in range(flux.shape[0]):
                flux[unit, :] = np.pi * 4. * flux * cspeed / wave**2
                flux = flux * (radius / (dvega * pctocm))**2
            self.vega_sed = kcorrect.template.SED(wave=wave, flux=flux)
            self.vega_sed.info = info

        vega_maggies = self.project(sed=self.vega_sed)
        if(vega_maggies > 0):
            self.vega2ab = - 2.5 * np.log10(vega_maggies)
        else:
            self.vega2ab = None
        return
