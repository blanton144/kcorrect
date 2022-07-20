#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: response.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import os
import re
import numpy as np
import kcorrect
import kcorrect.utils
import kcorrect.template
import scipy.interpolate as interpolate
import scipy.integrate as integrate
import pydl.pydlutils.yanny as yanny
import fitsio


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

    Returns all base names of files with the suffix '.par' in the
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
            m = re.match('^(.*)\\.par$', filename)
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
        self[response].frompar(filename='{n}.par'.format(n=response))
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

    If a file is given it is assumed to be in Yanny parameter
    format, readable and writeable by the class
    pydl.pydlutils.yanny.yanny

    The attribute interp() takes wavelength in Angstroms as its one
    positional argument.

    The effective wavelength is defined as described in Blanton &
    Roweis (2007).
"""
    def __init__(self, filename=None, wave=None, response=None):
        self.wave = wave
        self.response = response
        self.filename = filename
        self.solar_sed = None
        self.solar_magnitude = None
        self.vega_sed = None
        self.vega_magnitude = None
        self.lambda_eff = None
        self.interp = None
        if(self.filename is not None):
            self.frompar(filename)
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
        response = fitsio.read(filename, ext=ext)
        self.nwave = len(response)
        self.wave = response['wave']
        self.response = response['response']
        self._setup()
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
        to KCORRECT_DIR/python/kcorrect/data/responses
"""
        if(os.path.isabs(filename)):
            infilename = filename
        else:
            infilename = os.path.join(kcorrect.KCORRECT_DIR, 'data',
                                      'responses', filename)

        par = yanny.yanny(infilename)
        name = par.tables()[0]
        parstr = par[name]
        self.nwave = len(parstr)
        isort = np.argsort(parstr['lambda'])
        self.wave = parstr['lambda'][isort]
        self.response = parstr['pass'][isort]
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

        fitsio.write(filename, out, clobber=clobber)
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

        Assumes AB calibration.

        If the bandpass is outside the range of the solar model, 0 is returned.
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
