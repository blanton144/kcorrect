#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: kcorrect.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import os
import numpy as np
import kcorrect.fitter
import kcorrect.template
import astropy.cosmology
import astropy.units


class Kcorrect(kcorrect.fitter.Fitter):
    """K-correction object

    Parameters
    ----------

    abcorrect : bool
        correct maggies to AB (default False)

    responses : list of str
        names of input responses to base SED on

    templates : list of kcorrect.template.SED
        templates to use (if None uses v4 default template set)

    responses_out : list of str
        output responses for K-corrections (default to "responses")

    responses_map : list of str
        input responses to use for K-corrections (default to "responses")

    redshift_range : list of np.float32
        minimum and maximum redshifts (default [0., 2.])

    nredshift : int or np.int32
        number of redshifts in interpolation grid (default 4000)

    cosmo : astropy.cosmology.FLRW-like object
        object with distmod() method (default Planck18)

    Attributes
    ----------

    abcorrect : bool
        correct maggies to AB

    Amatrix : scipy.interpolate.interp1d object
        interpolation function for each template and input response

    AmatrixOut : scipy.interpolate.interp1d object
        interpolation function for each template and output response

    cosmo : astropy.cosmology.FLRW-like object
        object with luminosity_distance() method

    imap : ndarray of np.int32
        for each responses_map element, its corresponding index in responses

    nredshift : int or np.int32
        number of redshifts in interpolation grid

    redshift_range : list of np.float32
        minimum and maximum redshifts

    redshifts : ndarray of np.float32
        redshifts in grid

    responses : list of str
        [Nin] names of input responses to use

    responses_map : list of str
        [Nout] input responses to use for K-corrections

    responses_out : list of str
        [Nout] output responses for K-corrections

    templates : kcorrect.template.Template object
        templates to use

    Notes
    -----

    K-corrections are magnitude shifts to account for the difference
    in the observed bandpass R and a desired output bandpass Q,
    denoted K_QR (see Blanton & Roweis 2007). The default behavior
    finds the K-corrections between each input bandpass with itself,
    i.e. K_QQ.

    In detail:

    "responses" corresponds to the observed bandpasses R.

    "responses_out" corresponds to the observed bandpasses Q; it
    defaults to "responses".

    "responses_map" is the same length as "responses_out" and
    defines which bandpasses in "responses" to use for each
    output bandpass in "responses out". It defaults to "responses".

    Amatrix accepts a redshift as its argument and returns a matrix of
    shape [nresponses, ntemplates]. This matrix can be dotted into a
    set of coefficients for each template, and the result will be the
    observed bandpass maggies at the desired redshift for an SED
    corresponding to the coefficients. A ValueError results if the
    input redshift is outside redshift_range.

    AmatrixOut is similar but returns a [nresponses_out, ntemplates]
    matrix for the output bandpasses.
"""
    def __init__(self, responses=None, templates=None, responses_out=None,
                 responses_map=None, redshift_range=[0., 2.], nredshift=4000,
                 abcorrect=False, cosmo=None):

        # Read in templates
        if(templates is None):
            filename = os.path.join(os.getenv('KCORRECT_DIR'), 'python',
                                    'kcorrect', 'data', 'templates',
                                    'kcorrect-default-v4.fits')
            templates = kcorrect.template.Template(filename=filename)

        # Initatialize using Fitter initialization
        super().__init__(responses=responses, templates=templates,
                         redshift_range=redshift_range,
                         nredshift=nredshift, abcorrect=abcorrect)

        # Set up the Amatrix for the input responses
        self.set_Amatrix()

        # Set up the output responses
        if(responses_out is None):
            responses_out = self.responses
        if(responses_map is None):
            responses_map = self.responses

        if(len(responses_map) != len(responses_out)):
            raise ValueError("responses_map must have the same number of elements as responses_out")

        self.responses_out = responses_out
        self.responses_map = responses_map

        # Get index map to calculate kcorrection
        self.imap = np.zeros(len(self.responses_map), dtype=int)
        for i, response in enumerate(self.responses_map):
            try:
                self.imap[i] = self.responses.index(response)
            except ValueError:
                raise ValueError("responses_map must contain only responses defined in responses")

        # Set up the AmatrixOut for the output responses
        if(self.responses_out == self.responses):
            self.AmatrixOut = self.Amatrix
        else:
            self.AmatrixOut = self._calc_Amatrix(responses=self.responses_out)

        # Initialize cosmology used for derived properties and absmag
        if(cosmo is not None):
            self.cosmo = cosmo
        else:
            self.cosmo = astropy.cosmology.Planck18

        return

    def derived(self, redshift=None, coeffs=None, band_shift=0.):
        """Return derived quantities based on coefficients

        Parameters
        ----------

        redshift : ndarray of np.float32, or np.float32
            redshift

        coeffs : ndarray of np.float32
            coefficients for each template for each object

        band_shift : np.float32
            band shift to apply for output response

        Returns
        -------

        derived : dict()
            dictionary with derived quantities (see below)

        Notes
        -----

        The derived dictionary contains the following keys with the
        associated quantities:

            'mremain'  : ndarray of np.float32, or np.float32
                current stellar mass in solar masses

            'intsfh' : ndarray of np.float32, or np.float32
               current stellar mass in solar masses

            'mtol' : ndarray of np.float32, or np.float32
               mass-to-light ratio in each output band

            'b300' : ndarray of np.float32, or np.float32
               current (< 300 Myr) over past star formation

            'b1000' : ndarray of np.float32, or np.float32
               current (< 1 Gyr) over past star formation

            'metallicity' :
               metallicity in current stars

        All of these quantities should be taken with extreme caution
        and not accepted literally. After all, they are just the result
        of a 5-template fit to a few bandpasses. See Moustakas et al.
        (2013) for a comparison of the masses with other estimators.
"""
        (array, n, redshift, d1, d2,
         coeffs) = self._process_inputs(redshift=redshift, coeffs=coeffs)

        dm = self.cosmo.distmod(redshift).to_value(astropy.units.mag)
        dfactor = 10.**(0.4 * dm)

        intsfh = coeffs.dot(self.templates.intsfh) * dfactor
        mremain = coeffs.dot(self.templates.mremain) * dfactor
        metals = (coeffs.dot(self.templates.mremain * self.templates.mets) *
                  dfactor)

        metallicity = metals / mremain

        m300 = coeffs.dot(self.templates.m300) * dfactor
        b300 = m300 / intsfh
        m1000 = coeffs.dot(self.templates.m1000) * dfactor
        b1000 = m1000 / intsfh

        # Doesn't do solar yet
        f = kcorrect.response.ResponseDict()
        rmaggies = self.reconstruct_out(redshift=redshift,
                                        coeffs=coeffs,
                                        band_shift=band_shift)
        for ir, response in enumerate(self.responses):
            solar = 10.**(- 0.4 * f[response].solar_magnitude)
            rmaggies[..., ir] = rmaggies[..., ir] * dfactor / solar
        mtol = (np.outer(coeffs.dot(self.templates.mremain),
                         np.ones(len(self.responses), dtype=np.float32)) /
                rmaggies)

        outdict = dict()
        outdict['mremain'] = mremain
        outdict['intsfh'] = intsfh
        outdict['mtol'] = mtol
        outdict['b300'] = b300
        outdict['b1000'] = b1000
        outdict['metallicity'] = metallicity

        return(outdict)

    def reconstruct_out(self, redshift=None, coeffs=None, band_shift=0.):
        """Reconstruct output maggies associated with coefficients

        Parameters
        ----------

        redshift : np.float32
            redshift

        coeffs : ndarray of np.float32
            coefficients

        band_shift : np.float32
            blueshift to apply to reconstructed bandpasses

        Returns
        -------

        maggies : ndarray of np.float32
            AB maggies in each output band
"""
        return(self._reconstruct(Amatrix=self.AmatrixOut, redshift=redshift,
                                 coeffs=coeffs, band_shift=band_shift))

    def kcorrect(self, redshift=None, coeffs=None, band_shift=0.):
        """Return K-correction in all bands

        Parameters
        ----------

        redshift : ndarray of np.float32, or np.float32
            redshift for K-correction

        coeffs : ndarray of np.float32
            coefficients for each template for each object

        band_shift : np.float32
            shift to apply for output responses

        Returns
        -------

        kcorrect : ndarray of np.float32
            K-correction from input to output magnitudes
"""
        maggies_in = self.reconstruct(redshift=redshift, coeffs=coeffs)
        maggies_out = self.reconstruct_out(redshift=0. * redshift,
                                           coeffs=coeffs,
                                           band_shift=band_shift)

        kcorrect = - 2.5 * np.log10(maggies_out / maggies_in[..., self.imap])
        return(kcorrect)

    def absmag(self, maggies=None, ivar=None, redshift=None, coeffs=None,
               band_shift=0.):
        """Return absolute magnitude in output bands

        Parameters
        ----------

        redshift : ndarray of np.float32, or np.float32
            redshift(s) for K-correction

        maggies : ndarray of np.float32
            fluxes of each band in maggies

        ivar : ndarray of np.float32
            inverse variance of each band

        coeffs : ndarray of np.float32
            coefficients for each template for each object

        band_shift : np.float32
            shift to apply for output responses

        Returns
        -------

        absmag : ndarray of np.float32
            AB absolute magnitude in each band for each object

        Notes
        -----

        Returns the K-corrected absolute magnitude.

        Depends on having run fit_coeffs on a consistent set of
        maggies and ivars. If ivar=0 or the maggies are negative 
        for any band, it uses the reconstructed absolute magnitude.

        Determines the distance modulus with the object's "cosmo.distmod()"
        method. By default this is the Planck18 cosmology. This use
        
        Calls to_ab() method on input maggies to convert to AB.
"""
        (array, n, redshift, maggies, ivar,
         coeffs) = self._process_inputs(redshift=redshift, maggies=maggies,
                                        ivar=ivar, coeffs=coeffs)

        dm = self.cosmo.distmod(redshift).to_value(astropy.units.mag)
        k = self.kcorrect(redshift=redshift, coeffs=coeffs,
                          band_shift=band_shift)
        omaggies = self.reconstruct_out(redshift=redshift, coeffs=coeffs,
                                        band_shift=band_shift)
        use_maggies = maggies
        ibad = np.where((use_maggies <= 0.) | (ivar <= 0.))
        use_maggies[ibad] = omaggies[ibad]

        mags = - 2.5 * np.log10(use_maggies)

        if(array):
            dm = np.outer(dm, np.ones(len(self.responses), dtype=np.float32))

        absmag = mags - dm - k

        return(absmag)


class KcorrectSDSS(Kcorrect):
    """K-correction object for SDSS data

    Parameters
    ----------

    abcorrect : bool
        correct maggies to AB (default True)

    templates : list of kcorrect.template.SED
        templates to use (if None uses v4 default template set)

    responses : list of str
        names of input responses to base SED on (default to SDSS)

    responses_out : list of str
        output responses for K-corrections (default to "responses")

    responses_map : list of str
        input responses to use for K-corrections (default to "responses")

    redshift_range : list of np.float32
        minimum and maximum redshifts (default [0., 2.])

    nredshift : int or np.int32
        number of redshifts in interpolation grid (default 4000)

    cosmo : astropy.cosmology.FLRW-like object
        object with distmod() method (default Planck18)

    Attributes
    ----------

    abcorrect : bool
        correct maggies to AB

    Amatrix : scipy.interpolate.interp1d object
        interpolation function for each template and input response

    AmatrixOut : scipy.interpolate.interp1d object
        interpolation function for each template and output response

    cosmo : astropy.cosmology.FLRW-like object
        object with luminosity_distance() method

    imap : ndarray of np.int32
        for each responses_map element, its corresponding index in responses

    nredshift : int or np.int32
        number of redshifts in interpolation grid

    redshift_range : list of np.float32
        minimum and maximum redshifts

    redshifts : ndarray of np.float32
        redshifts in grid

    responses : list of str
        [Nin] names of input responses to use

    responses_map : list of str
        [Nout] input responses to use for K-corrections

    responses_out : list of str
        [Nout] output responses for K-corrections

    templates : kcorrect.template.Template object
        templates to use

    Notes
    -----

    responses defaults to ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0']

    This class provides the method fit_coeffs_asinh() to use SDSS-style
    asinh magnitudes (these are the magnitudes that the SDSS imaging
    reports).

    The to_ab() method is applied to the maggies input for absmag() and
    fit_coeffs() and fit_coeffs_asinh(), which adjusts from the SDSS system
    to the AB system.
"""
    def __init__(self, responses=['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0',
                                  'sdss_z0'], templates=None,
                 responses_out=None, responses_map=None,
                 redshift_range=None, nredshift=None, cosmo=None,
                 abcorrect=True):

        # Read in templates
        if(templates is None):
            filename = os.path.join(os.getenv('KCORRECT_DIR'), 'python',
                                    'kcorrect', 'data', 'templates',
                                    'kcorrect-default-v4.fits')
            templates = kcorrect.template.Template(filename=filename)

        # Initatialize using Kcorrect initialization
        super().__init__(responses=responses, templates=templates,
                         redshift_range=redshift_range,
                         nredshift=nredshift, abcorrect=abcorrect)
        return

    def to_ab(self, maggies=None, ivar=None):
        """Convert input maggies to AB

        Parameters
        ----------

        maggies : ndarray of np.float32
            array of fluxes in standard SDSS system
        
        ivar : ndarray of np.float32
            inverse variances in standard SDSS system (optional)

        Returns
        -------

        ab_maggies : ndarray of np.float32
            array of fluxes converted to AB

        ab_ivar : ndarray of np.float32
            inverse variances converted to AB (if ivar input)

        Notes
        -----

        Calls kcorrect.utils.sdss_ab_correct(), which does the following:

        Uses the AB conversions produced by D. Eisenstein, in his
        message sdss-calib/1152

            u(AB,2.5m) = u(database, 2.5m) - 0.036
            g(AB,2.5m) = g(database, 2.5m) + 0.012
            r(AB,2.5m) = r(database, 2.5m) + 0.010
            i(AB,2.5m) = i(database, 2.5m) + 0.028
            z(AB,2.5m) = z(database, 2.5m) + 0.040

        fit_coeffs() and absmag() call this on their inputs.
"""
        if(ivar is not None):
            maggies = kcorrect.utils.sdss_ab_correct(maggies=maggies,
                                                     ivar=ivar)
            return(maggies, ivar)
        else:
            maggies, ivar = kcorrect.utils.sdss_ab_correct(maggies=maggies,
                                                           ivar=ivar)
            return(maggies)

    def fit_coeffs_asinh(self, redshift=None, mag=None, mag_err=None,
                         ext=None):
        """Fit coefficients to asinh mags

        Parameters
        ----------

        redshift : np.float32 or ndarray of np.float32
            [N] or scalar redshift(s)

        mag : ndarray of np.float32
            [N, 5] or [5] asinh magnitudes of each SDSS band

        mag_err : ndarray of np.float32
            [N, 5] or [5] inverse variance of each band

        extinction : ndarray of np.float32
            [N, 5] or [5] Galactic extinction in each band

        Returns
        -------

        coeffs : ndarray of np.float32
            coefficients for each template

        Notes
        -----

        Converts mag, mag_err, and extinction to extinction-corrected
        maggies and ivar, and then calls to_ab() method to create
        AB maggies and ivar.

        If redshift is an array, even with just one element, coeffs is
        returned as an [nredshift, ntemplate] array.

        Otherwise coeffs is returned as an [ntemplate] array.
"""
        if(redshift is None):
            raise TypeError("Must specify redshift to fit coefficients")

        maggies, ivar = kcorrect.utils.sdss_asinh_to_maggies(maggies=mag,
                                                             ivar=mag_err,
                                                             extinction=extinction)

        coeffs = self.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)
        return(coeffs)
