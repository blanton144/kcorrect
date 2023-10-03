#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: template.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interpolate


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

    ext : str
        extension from which to read FLUX in FITS file

    Attributes
    ----------

    binimage : bool
        if True, files use WAVE and FLUX extensions

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

    If binimage is True, then instead of a FLUX HDU table, there should
    be WAVE and FLUX HDUs with binary images.

"""
    def __init__(self, filename=None, wave=None, flux=None, ext='FLUX'):
        self.restframe_wave = wave
        self.restframe_flux = flux
        self.filename = filename
        self.info = dict()
        self.binimage = False
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
        self.binimage = ('FLUX' in sed_hdus) & ('WAVE' in sed_hdus)
        if(self.binimage):
            self.restframe_wave = sed_hdus['WAVE'].data
            self.restframe_flux = sed_hdus['FLUX'].data
        else:
            sed = sed_hdus[ext].data
            self.restframe_wave = sed['wave']
            self.restframe_flux = sed['flux']
        if(len(self.restframe_flux.shape) > 1):
            nsed = self.restframe_flux.shape[-2]
        else:
            nsed = 1
        self.restframe_wave = np.squeeze(self.restframe_wave)
        nwave = len(self.restframe_wave)
        self.restframe_flux = self.restframe_flux.reshape(nsed, nwave)
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

            wave : an [nwave] array of restframe wavelengths in Angstrom

            flux : an [nsed, nwave] array of restframe fluxes

        If binimage is set for this object, instead write
        two HDUs.
"""
        if(self.binimage):
            hdul = fits.HDUList()
            whdu = fits.ImageHDU(self.restframe_wave, name='WAVE')
            hdul.append(whdu)
            fhdu = fits.ImageHDU(self.restframe_flux, name='FLUX')
            fhdu.header['BINIMAGE'] = 'T'
            hdul.append(fhdu)
            hdul.writeto(filename, overwrite=clobber)
        else:
            out = np.zeros(1, self.sed_dtype())
            out['wave'] = self.restframe_wave
            out['flux'] = self.restframe_flux

            hdu = fits.BinTableHDU(out, name=ext)
            hdu.writeto(filename, overwrite=clobber)
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

    def plot(self, indx=None, wavelim=None):
        """Simple matplotlib plot of template(s)

        Parameters
        ----------

        indx : np.int32 or ndarray of np.int32
            index of template(s) to plot

        wavelim : list or ndarray of np.float32
            [2] wavelength limits of plot

        Notes
        -----

        Plots log of wavelength and log of lambda * flux_lambda.
"""
        if(indx is None):
            indx = np.arange(self.nsed, dtype=int)

        indxs = np.int32(indx)
        if(indxs.ndim == 0):
            indxs = np.array([indxs], dtype=np.int32)

        for i in indxs:
            ok = self.flux[i, :] > 0.
            plt.plot(np.log10(self.wave[ok]),
                     np.log10(self.wave[ok] * self.flux[i, ok]),
                     linewidth=1, alpha=0.5)

        if(wavelim is not None):
            iwave = np.where((self.wave > wavelim[0]) &
                             (self.wave < wavelim[1]))[0]
            plt.xlim(np.log10(np.float32(wavelim)))
            lfl = self.flux * np.outer(np.ones(self.flux.shape[0],
                                               dtype=np.float32),
                                       self.wave)
            yvals = lfl[:, iwave][indxs, :].flatten()
            yvals = yvals[yvals > 0.]
            ymin = np.log10(yvals.min()) - 0.1
            ymax = np.log10(yvals.max()) + 0.1
            plt.ylim([ymin, ymax])

        plt.xlabel('$\\log_{10}$ $\\lambda$  (Angstroms)')
        plt.ylabel('$\\log_{10}$ $\\lambda f_\\lambda$ (erg s$^{-1}$ cm$^{-2}$ M$_{\\odot}^{-1}$ at 10 pc)')
        return


class Template(SED):
    """Spectral energy distribution template(s)

    Parameters
    ----------

    filename : str
        name of FITS file to read from

    binimage : bool
        if True, read in WAVE and FLUX extensions as binary images

    ext : str
        extension from which to read FLUX in FITS file

    Attributes
    ----------

    binimage : bool
        if True, read in WAVE and FLUX extensions as binary images

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
       M50 : an [nsed]-array with mass formed within 50 My in solar units
       M300 : an [nsed]-array with mass formed within 300 My in solar units
       M1000 : an [nsed]-array with mass formed within 1 Gy in solar units

    If binimage is True, then instead of a FLUX HDU table, there should 
    be WAVE and FLUX HDUs with binary images.
"""
    def __init__(self, filename=None, ext='FLUX'):
        super().__init__(filename=filename, ext=ext)

        hdul = fits.open(filename)
        self.info['filename'] = filename
        self.intsfh = hdul['INTSFH'].data
        self.mremain = hdul['MREMAIN'].data
        self.mets = hdul['METS'].data
        self.m300 = hdul['M300'].data
        self.m1000 = hdul['M1000'].data
        try:
            self.m50 = hdul['M50'].data
        except:
            self.m50 = self.m300 * 0.

        return

    def tofits(self, filename=None, clobber=True):
        """Write template set to FITS

        Parameters
        ----------

        filename : str
            file to write to

        clobber : bool
            if True, overwrite existing file

"""
        super().tofits(filename=filename, clobber=clobber)
        hdul = fits.open(filename, mode='update')
        hdu = fits.ImageHDU(self.intsfh, name='INTSFH')
        hdul.append(hdu)
        hdu = fits.ImageHDU(self.mremain, name='MREMAIN')
        hdul.append(hdu)
        hdu = fits.ImageHDU(self.mets, name='METS')
        hdul.append(hdu)
        hdu = fits.ImageHDU(self.m300, name='M300')
        hdul.append(hdu)
        hdu = fits.ImageHDU(self.m50, name='M50')
        hdul.append(hdu)
        hdu = fits.ImageHDU(self.m1000, name='M1000')
        hdul.append(hdu)
        hdul.close()
