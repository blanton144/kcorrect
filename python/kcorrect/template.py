#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: template.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import numpy as np
import scipy.interpolate as interpolate
import astropy.io.fits as fits


class SED(object):
    """Spectral energy distribution(s)

    Parameters
    ----------

    filename : str
        name of FITS file to read from

    wave : ndarray of np.float32
        rest frame wavelength grid in Angstroms

    flux : ndarray of np.float32
        [nsed, nwave] rest frame flux grid in erg/cm^2/s/A at 10pc

    Attributes
    ----------

    filename : str
        name of FITS file associated with this SED

    flux : ndarray of np.float32
        [nsed, nwave] flux grid in erg/cm^2/s/A at 10pc and currently set redshift

    info : dict
        dictionary for storing assorted metadata associated with spectra

    nsed : np.float32
        number of SEDs

    nwave : np.int32
        number of wavelengths in grid

    redshift : np.float32
        redshift of SED

    restframe_flux : ndarray of np.float32
        [nsed, nwave] rest frame flux grid in erg/cm^2/s/A at 10pc

    restframe_wave : ndarray of np.float32
        [nwave] rest frame wavelength grid in Angstroms

    wave : ndarray of np.float32
        [nwave] wavelength grid in Angstroms

    Notes
    -----

    The fluxes are a bit funnily defined, in terms of the flux that
    the galaxy would have at 10pc, analogous to an absolute magnitude.
    When the redshift is applied, the bolometric flux is conserved
    (i.e. there is no luminosity distance applied, it is a pure 
    redshifting of the spectrum).

    If filename is set, overrides wave and flux.

    """
    def __init__(self, filename=None, wave=None, flux=None, ext='FLUX'):
        self.restframe_wave = wave
        self.restframe_flux = flux
        self.filename = filename
        self.info = dict()
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
        sed_hdus = fits.open(filename)
        sed = sed_hdus[ext].data
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

        hdu = fits.BinTableHDU(out, name=ext)
        hdu.writeto(filename, overwrite=True)
        return

    def set_redshift(self, redshift=0.):
        """Set redshift of SED

        Parameters
        ----------

        redshift : np.float32
            redshift to shift to

        Notes
        -----

        Conserves bolometric integral of flux.
"""
        self.wave = self.restframe_wave * (1. + redshift)
        self.flux = self.restframe_flux / (1. + redshift)
        self.redshift = redshift
        self.set_interp()
        return


class Template(SED):
    """Spectral energy distribution template(s)

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

    info : dict
        dictionary for storing assorted metadata associated with spectra

    intsfh : ndarray of np.float32
        [nsed] integrated star formation history in solar masses

    m300 : ndarray of np.float32
        [nsed] stars formed last 300 million years in solar masses

    m1000 : ndarray of np.float32
        [nsed] stars formed last billion years in solar masses

    mets : ndarray of np.float32
        [nsed] metallicity in current stars and stellar remnants

    mremain : ndarray of np.float32
        [nsed] remaining mass in stars and remnants in solar masses

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

    The file should have the HDUs:

       FLUX : an ndarray with two tags: 
         'wave' : an [nwave]-array of np.float32 with wavelength in Ang.
         'flux' : an [nsed, nwave]-array of np.float32 with flux in erg/s/cm^2/A
       METS : an [nsed]-array with metallicity
       INTSFH : an [nsed]-array with integrated SF in solar units
       MREMAIN : an [nsed]-array with current stellar mass in solar units
       M300 : an [nsed]-array with mass formed within 300 My in solar units
       M1000 : an [nsed]-array with mass formed within 1 Gy in solar units

"""
    def __init__(self, filename=None, wave=None, flux=None, ext='FLUX'):
        super().__init__(filename=filename)

        hdul = fits.open(filename)
        self.info['filename'] = filename
        self.intsfh = hdul['INTSFH'].data
        self.mremain = hdul['MREMAIN'].data
        self.mets = hdul['METS'].data
        self.m300 = hdul['M300'].data
        self.m1000 = hdul['M1000'].data

        return
