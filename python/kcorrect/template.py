#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: template.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import numpy as np
import scipy.interpolate as interpolate
import fitsio


class SED(object):
    """Spectral energy distribution(s)

    Parameters
    ----------

    filename : str
        name of FITS file to read from

    wave : ndarray of np.float32
        rest frame wavelength grid in Angstroms

    flux : ndarray of np.float32
        [nsed, nwave] rest frame flux grid in erg/cm^2/s/A

    Attributes
    ----------

    filename : str
        name of FITS file associated with this SED

    flux : ndarray of np.float32
        [nsed, nwave] flux grid in erg/cm^2/s/A

    nsed : np.float32
        number of SEDs

    nwave : np.int32
        number of wavelengths in grid
    
    redshift : np.float32
        redshift of SED

    restframe_flux : ndarray of np.float32
        [nsed, nwave] rest frame flux grid in erg/cm^2/s/A

    restframe_wave : ndarray of np.float32
        [nwave] rest frame wavelength grid in Angstroms

    wave : ndarray of np.float32
        [nwave] wavelength grid in Angstroms

    Notes
    -----

    If filename is set, overrides wave and flux.
"""
    def __init__(self, filename=None, wave=None, flux=None, ext='FLUX'):
        self.restframe_wave = wave
        self.restframe_flux = flux
        self.filename = filename
        if(self.filename is not None):
            self.fromfits(filename, ext=ext)
            return
        self._setup()
        return

    def _setup(self):
        """Set up after restframe_wave and restframe_flux are read in"""
        self.redshift = 0.
        if((self.restframe_wave is not None) &
           (self.restframe_flux is not None)):
            self.nwave = len(self.restframe_wave)
            if(len(self.restframe_flux.shape) == 1):
                self.restframe_flux = self.restframe_flux.reshape(1, self.nwave)
            else:
                self.restframe_flux = self.restframe_flux
            self.nsed = self.restframe_flux.shape[0]
        self.wave = self.restframe_wave
        self.flux = self.restframe_flux
        self.set_interp()
        return

    def sed_dtype(self):
        """Returns numpy dtype for SED"""
        sed_dtype = np.dtype([('wave', np.float32, self.nwave),
                              ('flux', np.float32, (self.nsed,
                                                    self.nwave))])
        return(sed_dtype)

    def set_interp(self):
        """Sets attribute interp to interpolation function"""
        if((self.wave is None) | (self.flux is None)):
            self.interp = None
            return
        self.interp = interpolate.interp1d(self.wave, self.flux,
                                           kind='cubic',
                                           bounds_error=False,
                                           fill_value=0.)
        return

    def fromfits(self, filename=None, ext='FLUX'):
        """Read SED from FITS files

        Parameters
        ----------

        filename : str
            input file name

        ext : str or int
            extension to read from

        Notes
        -----

        The FITS table should have two columns:

           wave - an [nwave] array of wavelengths in Angstrom
           flux - an [nsed, nwave] array of fluxes

        Only imports the first row of the FITS table.
"""
        sed = fitsio.read(filename, ext=ext)
        self.restframe_wave = np.squeeze(sed['wave'])
        if(len(sed['flux'][0].shape) > 1):
            self.restframe_flux = sed['flux'][0]
        else:
            self.restframe_flux = sed['flux'].reshape(1, len(sed['flux'][0]))
        self._setup()
        return

    def tofits(self, filename=None, ext='FLUX', clobber=True):
        """Write SED to FITS files

        Parameters
        ----------

        filename : str
            output file name

        ext : str or int
            extension to write to
        
        clobber : bool
            whether to clobber the existing file or add an HDU

        Notes
        -----

        Writes a one-row FITS table with two columns:

           wave - an [nwave] array of restframe wavelengths in Angstrom
           flux - an [nsed, nwave] array of restframe fluxes

        Only imports the first row of the FITS table.
"""
        out = np.zeros(1, self.sed_dtype())
        out['wave'] = self.restframe_wave
        out['flux'] = self.restframe_flux

        fitsio.write(filename, out, extname=ext, clobber=clobber)
        return

    def set_redshift(self, redshift=0.):
        """Set redshift of SED

        Parameters
        ----------

        redshift : np.float32
            redshift to shift to

        Notes
        -----

        Conserves integral of flux.
"""
        self.wave = self.restframe_wave * (1. + redshift)
        self.flux = self.restframe_flux / (1. + redshift)
        self.redshift = redshift
        self.set_interp()
        return
