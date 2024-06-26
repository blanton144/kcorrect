#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: kcorrect.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import os

import astropy.cosmology
import astropy.io.fits as fits
import astropy.units
import numpy as np

import kcorrect
import kcorrect.fitter
import kcorrect.template


class Kcorrect(kcorrect.fitter.Fitter):
    """K-correction object

    Parameters
    ----------

    abcorrect : bool
        correct maggies to AB (default False)

    filename : str
        input file to define kcorrect object (overrides responses, responses_out, responses_map, templates, redshift_range, nredshift, abcorrect)

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
    output bandpass in "responses_out". It defaults to "responses".

    Amatrix accepts a redshift as its argument and returns a matrix of
    shape [nresponses, ntemplates]. This matrix can be dotted into a
    set of coefficients for each template, and the result will be the
    observed bandpass maggies at the desired redshift for an SED
    corresponding to the coefficients. A ValueError results if the
    input redshift is outside redshift_range.

    AmatrixOut is similar but returns a [nresponses_out, ntemplates]
    matrix for the output bandpasses.

    Once defined, a Kcorrect object can be output to a FITS file with 
    the method tofits(), and reimported with fromfits(). The filename
    parameter expects this format.
"""
    def __init__(self, filename=None, responses=None, templates=None,
                 responses_out=None, responses_map=None,
                 redshift_range=[0., 2.], nredshift=4000,
                 abcorrect=False, cosmo=None):

        if(filename is not None):
            self.fromfits(filename=filename)
        else:
            # Read in templates
            if(templates is None):
                tfilename = os.path.join(kcorrect.KCORRECT_DIR, 'data',
                                         'templates',
                                         'kcorrect-default-v4.fits')
                templates = kcorrect.template.Template(filename=tfilename)

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

            # Set up the AmatrixOut for the output responses
            if(self.responses_out == self.responses):
                self.AmatrixOut = self.Amatrix
            else:
                self.AmatrixOut = self._calc_Amatrix(responses=self.responses_out)

        # Get index map to calculate kcorrection
        self.imap = np.zeros(len(self.responses_map), dtype=int)
        for i, response in enumerate(self.responses_map):
            try:
                self.imap[i] = self.responses.index(response)
            except ValueError:
                raise ValueError("responses_map must contain only responses defined in responses")

        # Initialize cosmology used for derived properties and absmag
        if(cosmo is not None):
            self.cosmo = cosmo
        else:
            self.cosmo = astropy.cosmology.Planck18

        return

    def derived(self, redshift=None, coeffs=None, band_shift=0., distance=None):
        """Return derived quantities based on coefficients

        Parameters
        ----------

        redshift : ndarray of np.float32, or np.float32
            [ngalaxy] redshift

        coeffs : ndarray of np.float32
            [ngalaxy, ntemplates] coefficients for each template for each object

        band_shift : np.float32
            band shift to apply for output response

        distance : ndarray of np.float32 or np.float32
            [ngalaxy] luminosity distance in Mpc

        Returns
        -------

        derived : dict()
            dictionary with derived quantities (see below)

        Notes
        -----

        If distance is not specified, then it is derived from
        the redshift assuming the cosmology in the cosmo attribute.

        The derived dictionary contains the following keys with the
        associated quantities:

        'mremain'  : ndarray of np.float32, or np.float32
             [ngalaxy] current stellar mass in solar masses

        'intsfh' : ndarray of np.float32, or np.float32
             [ngalaxy] current stellar mass in solar masses

        'mtol' : ndarray of np.float32, or np.float32
             [ngalaxy] mass-to-light ratio in each output band

        'b50' : ndarray of np.float32, or np.float32
             [ngalaxy] current (< 50 Myr) over past star formation

        'b300' : ndarray of np.float32, or np.float32
             [ngalaxy] current (< 300 Myr) over past star formation

        'b1000' : ndarray of np.float32, or np.float32
             [ngalaxy] current (< 1 Gyr) over past star formation

        'metallicity' :
             [ngalaxy] metallicity in current stars

        All of these quantities should be taken with extreme caution
        and not accepted literally. After all, they are just the result
        of a template fit to a few bandpasses. 

        For the default template set, see Moustakas et al. (2013) for
        a comparison of the masses with other estimators.

        For responses where solar_magnitude is not defined, mtol is in
        per maggy units (not per solar luminosity)
"""
        (array, n, redshift, d1, d2,
         coeffs) = self._process_inputs(redshift=redshift, coeffs=coeffs)

        if(distance is not None):
            (array, n, distance, d1, d2,
             d3) = self._process_inputs(redshift=distance, coeffs=coeffs)
            dfactor = (distance / 1.e-5)**2  # factor relative to 10 pc
        else:
            dm = self.cosmo.distmod(redshift).to_value(astropy.units.mag)
            dfactor = 10.**(0.4 * dm)

        intsfh = coeffs.dot(self.templates.intsfh) * dfactor
        mremain = coeffs.dot(self.templates.mremain) * dfactor
        metals = (coeffs.dot(self.templates.mremain * self.templates.mets) *
                  dfactor)
        m50 = coeffs.dot(self.templates.m50) * dfactor
        m300 = coeffs.dot(self.templates.m300) * dfactor
        m1000 = coeffs.dot(self.templates.m1000) * dfactor

        if(array):
            ok = mremain > 0.
            metallicity = np.zeros(len(redshift), dtype=np.float32)
            b50 = np.zeros(len(redshift), dtype=np.float32)
            b300 = np.zeros(len(redshift), dtype=np.float32)
            b1000 = np.zeros(len(redshift), dtype=np.float32)
            metallicity[ok] = metals[ok] / mremain[ok]
            b50[ok] = m50[ok] / intsfh[ok]
            b300[ok] = m300[ok] / intsfh[ok]
            b1000[ok] = m1000[ok] / intsfh[ok]
        else:
            if(mremain > 0.):
                metallicity = metals / mremain
                b50 = m50 / intsfh
                b300 = m300 / intsfh
                b1000 = m1000 / intsfh
            else:
                metallicity = np.float32(0.)
                b50 = np.float32(0.)
                b300 = np.float32(0.)
                b1000 = np.float32(0.)

        f = kcorrect.response.ResponseDict()
        if(array):
            zero_redshift = np.zeros(len(redshift), dtype=np.float32)
        else:
            zero_redshift = np.float32(0.)
        rmaggies_solar = self.reconstruct_out(redshift=zero_redshift,
                                              coeffs=coeffs,
                                              band_shift=band_shift)

        for ir, response in enumerate(self.responses_out):
            solar = 10.**(- 0.4 * (f[response].solar_magnitude + dm))
            rmaggies_solar[..., ir] = rmaggies_solar[..., ir] / solar

        mtol = np.zeros(rmaggies_solar.shape, dtype=np.float32)
        ok = rmaggies_solar > 0.
        if(array):
            mtol[ok] = (np.outer(mremain,
                                 np.ones(len(self.responses), dtype=np.float32))[ok] /
                        rmaggies_solar[ok])
        else:
            mtol[ok] = mremain / rmaggies_solar[ok]

        outdict = dict()
        outdict['mremain'] = mremain
        outdict['intsfh'] = intsfh
        outdict['mtol'] = mtol
        outdict['b50'] = b50
        outdict['b300'] = b300
        outdict['b1000'] = b1000
        outdict['metallicity'] = metallicity

        return(outdict)

    def derived_mc(self, redshift=None, coeffs_mc=None, band_shift=0., distance=None):
        """Return derived quantities based on coefficients

        Parameters
        ----------

        redshift : ndarray of np.float32, or np.float32
            [ngalaxy] redshift

        coeffs_mc : ndarray of np.float32
            [ngalaxy, ntemplates, mc] coefficients for each template for each object

        band_shift : np.float32
            band shift to apply for output response

        distance : ndarray of np.float32 or np.float32
            [ngalaxy] luminosity distance in Mpc

        Returns
        -------

        derived : dict()
            dictionary with derived quantities (see below)

        Notes
        -----

        Relies on repeated calls to the kcorrect.kcorrect.Kcorrect.derived() method.

        The derived dictionary contains the following keys with the
        associated quantities:

        'mremain'  : ndarray of np.float32, or np.float32
             [ngalaxy, mc] current stellar mass in solar masses

        'intsfh' : ndarray of np.float32, or np.float32
             [ngalaxy, mc] current stellar mass in solar masses

        'mtol' : ndarray of np.float32, or np.float32
             [ngalaxy, mc] mass-to-light ratio in each output band

        'b50' : ndarray of np.float32, or np.float32
             [ngalaxy] current (< 50 Myr) over past star formation

        'b300' : ndarray of np.float32, or np.float32
             [ngalaxy, mc] current (< 300 Myr) over past star formation

        'b1000' : ndarray of np.float32, or np.float32
             [ngalaxy, mc] current (< 1 Gyr) over past star formation

        'metallicity' :
             [ngalaxy, mc] metallicity in current stars
"""
        outdict = None
        mc = coeffs_mc.shape[-1]
        for imc in np.arange(mc, dtype=np.int32):
            coeffs_curr = coeffs_mc[..., imc]
            derived = self.derived(redshift=redshift, band_shift=band_shift,
                                   distance=distance, coeffs=coeffs_curr)
            if(outdict is None):
                outdict = dict()
                for k in derived:
                    shp = derived[k].shape
                    shpmc = shp + (mc,)
                    outdict[k] = np.zeros(shpmc, dtype=np.float32)
            for k in derived:
                outdict[k][..., imc] = derived[k]
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
        (array, n, redshift, d1,
         d2, coeffs) = self._process_inputs(redshift=redshift, maggies=None,
                                            ivar=None, coeffs=coeffs)

        # maggies associated with bandpass R at observed z
        maggies_in = self.reconstruct(redshift=redshift, coeffs=coeffs)

        # maggies associated with bandpass Q at z=0
        maggies_out = self.reconstruct_out(redshift=0. * redshift,
                                           coeffs=coeffs,
                                           band_shift=band_shift)

        # Now carefully take the K-correction (avoiding cases
        # where the coefficients are zero)
        kcorrect = np.zeros(maggies_out.shape, dtype=np.float32)
        if(array):
            # check if coefficients are ever zero
            notzero = coeffs.sum(axis=-1) > 0
            nz_maggies_in = maggies_in[notzero, :]
            kcorrect[notzero, :] = (- 2.5 *
                                    np.log10(nz_maggies_in[:, self.imap] /
                                             maggies_out[notzero, :]))
        else:
            notzero = coeffs.sum(axis=-1) > 0
            if(notzero):
                kcorrect = - 2.5 * np.log10(maggies_in[self.imap] /
                                            maggies_out)

        return(kcorrect)

    def absmag(self, maggies=None, ivar=None, redshift=None, coeffs=None,
               band_shift=0., distance=None, reconstruct=False, limit=False,
               kcorrect=False):
        """Return absolute magnitude in output bands

        Parameters
        ----------

        redshift : ndarray of np.float32, or np.float32
            [ngalaxy] redshift(s) for K-correction

        maggies : ndarray of np.float32
            [ngalaxy, nbands] fluxes of each band in maggies

        ivar : ndarray of np.float32
            [ngalaxy, nbands] inverse variance of each band

        coeffs : ndarray of np.float32
            [ngalaxy, ntemplates] coefficients for each template for each object

        band_shift : np.float32
            shift to apply for output responses

        distance : ndarray of np.float32 or np.float32
            [ngalaxy] distance in Mpc (or None)

        reconstruct : bool
            if set, return absmag_reconstruct

        limit : bool
            if set, return absmag_limit

        kcorrect : bool
            if set, return kcorrect

        Returns
        -------

        absmag : ndarray of np.float32
            [ngalaxy, nbands] AB absolute magnitude in each band for each object

        absmag_reconstruct : ndarray of np.float32
            [ngalaxy, nbands] reconstructed AB absolute magnitude from SED fit (if reconstruct set)
        absmag_limit : ndarray of np.float32
            [ngalaxy, nbands] 1-sigma bright limit on absolute mag from error

        kcorrect : ndarray of np.float32
            [ngalaxy, nbands] K-corrections used (if kcorrect is True)

        Notes
        -----

        If distance is None, the distance is derived from the redshift
        using the cosmo attribute.

        Depends on having run fit_coeffs on a consistent set of
        maggies and ivars. 

        Determines the distance modulus with the object's "cosmo.distmod()"
        method. By default this is the Planck18 cosmology. 

        Returns the K-corrected absolute magnitude (or -9999 if there
        ivar=0, magggies are negative, or there is no valid value)

        The "reconstructed" absolute magnitude is calculated from the
        kcorrect nonnegative fit instead of being a K-correction of the data.
        These values are -9999 if there is no valid value (usually meaning
        all of the coefficients input are zero).

        The absolute magnitude "limit" is set to the 1-sigma limit from 
        the error in each band (using the K-correction from the real fit).
        It is returned for every object with ivar != 0 and positive
        coefficients, and is -9999 otherwise.

        If abcorrect is True, calls to_ab() method on input maggies
        to convert to AB.
"""
        (array, n, redshift, maggies, ivar,
         coeffs) = self._process_inputs(redshift=redshift, maggies=maggies,
                                        ivar=ivar, coeffs=coeffs)

        if(distance is not None):
            (array, n, distance, d1, d2,
             d3) = self._process_inputs(redshift=distance, coeffs=coeffs)
            dfactor = (distance / 1.e-5)**2  # factor relative to 10 pc
            dm = 2.5 * np.log10(dfactor)
        else:
            dm = self.cosmo.distmod(redshift).to_value(astropy.units.mag)

        k = self.kcorrect(redshift=redshift, coeffs=coeffs,
                          band_shift=band_shift)

        use_maggies = maggies[..., self.imap]
        use_ivar = ivar[..., self.imap]

        gd = np.where((use_maggies > 0.) & (use_ivar > 0.))
        bd = np.where((use_maggies <= 0.) | (use_ivar <= 0.))

        mags = np.zeros(use_maggies.shape, dtype=np.float32)
        mags[gd] = - 2.5 * np.log10(use_maggies[gd])

        absmag = np.zeros(use_maggies.shape, dtype=np.float32)

        if(array):
            dm = np.outer(dm, np.ones(len(self.responses_out),
                                      dtype=np.float32))
            absmag[gd] = mags[gd] - dm[gd] - k[gd]
        else:
            absmag[gd] = mags[gd] - dm - k[gd]

        absmag[bd] = - 9999.

        if(reconstruct):
            omaggies = self.reconstruct_out(redshift=redshift, coeffs=coeffs,
                                            band_shift=band_shift)
            gd = np.where(omaggies > 0.)
            bd = np.where(omaggies <= 0.)
            omags = np.zeros(omaggies.shape, dtype=np.float32)
            omags[gd] = - 2.5 * np.log10(omaggies[gd])
            absmag_reconstruct = np.zeros(omaggies.shape, dtype=np.float32) - 9999.
            if(array):
                absmag_reconstruct[gd] = omags[gd] - dm[gd] - k[gd]
            else:
                absmag_reconstruct[gd] = omags[gd] - dm - k[gd]
            absmag_reconstruct[bd] = - 9999.

        if(limit):
            absmag_limit = np.zeros(use_maggies.shape, dtype=np.float32) - 9999.
            gd = np.where(use_ivar > 0.)
            bd = np.where(use_ivar <= 0.)
            
            lmags = np.zeros(use_maggies.shape, dtype=np.float32)
            lmags[gd] = - 2.5 * np.log10(1. / np.sqrt(use_ivar[gd]))

            if(array):
                absmag_limit[gd] = lmags[gd] - dm[gd] - k[gd]
            else:
                absmag_limit[gd] = lmags[gd] - dm - k[gd]
            absmag_limit[bd] = - 9999.

        return_list = [absmag]
        if(reconstruct):
            return_list.append(absmag_reconstruct)
        if(limit):
            return_list.append(absmag_limit)
        if(kcorrect):
            return_list.append(k)
        if(len(return_list) == 1):
            return(return_list[0])
        else:
            return(tuple(return_list))

    def absmag_mc(self, maggies_mc=None, ivar=None, redshift=None,
                  coeffs_mc=None, band_shift=0., distance=None, reconstruct=False,
                  kcorrect=False):
        """Return absolute magnitude in output bands for Monte Carlo results

        Parameters
        ----------

        redshift : ndarray of np.float32, or np.float32
            [ngalaxy] redshift(s) for K-correction

        maggies_mc : ndarray of np.float32
            [ngalaxy, nbands, mc] fluxes of each band in maggies

        ivar : ndarray of np.float32
            [ngalaxy, nbands] inverse variance of each band

        coeffs_mc : ndarray of np.float32
            [ngalaxy, ntemplates, mc] coefficients for each template for each object

        band_shift : np.float32
            shift to apply for output responses

        distance : ndarray of np.float32 or np.float32
            [ngalaxy] distance in Mpc (or None)

        reconstruct : bool
            if set, return absmag_reconstruct

        kcorrect : bool
            if set, return kcorrect

        Returns
        -------

        absmag : ndarray of np.float32
            [ngalaxy, nbands, mc] AB absolute magnitude in each band for each object

        absmag_reconstruct : ndarray of np.float32
            [ngalaxy, nbands, mc] reconstructed AB absolute magnitude from SED fits (if reconstruct set)

        kcorrect : ndarray of np.float32
            [ngalaxy, nbands, mc] K-corrections used (if kcorrect is True)

        Notes
        -----

        Relies on multiple calls to kcorrect.kcorrect.Kcorrect.absmag() method.
"""
        mc = coeffs_mc.shape[-1]
        absmag_mc = np.zeros(maggies_mc.shape, dtype=np.float32)
        if(reconstruct):
            absmag_reconstruct_mc = np.zeros(maggies_mc.shape, dtype=np.float32)
        if(kcorrect):
            kcorrect_mc = np.zeros(maggies_mc.shape, dtype=np.float32)
        for imc in np.arange(mc, dtype=np.int32):
            out = self.absmag(redshift=redshift, maggies=maggies_mc[..., imc],
                              ivar=ivar, coeffs=coeffs_mc[..., imc], band_shift=band_shift,
                              distance=distance, reconstruct=reconstruct, kcorrect=kcorrect)
            if((reconstruct is False) & (kcorrect is False)):
                absmag_mc[..., imc] = out
            else:
                absmag_mc[..., imc] = out[0]
                iout = 1
                if(reconstruct):
                    absmag_reconstruct_mc[..., imc] = out[iout]
                    iout += 1
                if(kcorrect):
                    kcorrect_mc[..., imc] = out[iout]
                    iout += 1

        return_list = [absmag_mc]
        if(reconstruct):
            return_list.append(absmag_reconstruct_mc)
        if(kcorrect):
            return_list.append(kcorrect_mc)
        if(len(return_list) == 1):
            return(return_list[0])
        else:
            return(tuple(return_list))

        raise RuntimeError("Should not reach this point")

    def tofits(self, filename=None):
        """Output calculated information to FITS

        Parameters
        ----------

        filename : str
            name of output file

        Notes
        -----

        Overwrites file if it already exists.

        Stores HDUs:

        * DATA : single row table with: 'nredshift', 'redshift_range', 'abcorrect'
        * REDSHIFTS : redshift grid
        * RESPONSES : array of input response names
        * RESPONSES_OUT : array of output response names
        * RESPONSES_MAP : array of response mapping
        * A : A-matrix data
        * AOUT : A-matrix output data
"""
        self.templates.tofits(filename=filename)
        
        hdul = fits.open(filename, mode='update')
        data_dtype = np.dtype([('nredshift', np.int32),
                               ('redshift_range', np.int32, 2),
                               ('abcorrect', bool)])
        data = np.zeros(1, dtype=data_dtype)
        data['nredshift'] = self.nredshift
        data['redshift_range'] = self.redshift_range
        data['abcorrect'] = self.abcorrect
        hdu = fits.BinTableHDU(data, name='DATA')
        hdul.append(hdu)
        hdu = fits.ImageHDU(self.redshifts, name='REDSHIFTS')
        hdul.append(hdu)

        responses_dtype = np.dtype([('responses', np.compat.unicode, 200)])
        responses = np.zeros(len(self.responses), dtype=responses_dtype)
        responses['responses'] = np.array(self.responses)
        hdu = fits.BinTableHDU(responses, name='RESPONSES')
        hdul.append(hdu)

        responses_dtype = np.dtype([('responses_out', np.compat.unicode, 200)])
        responses_out = np.zeros(len(self.responses_out), dtype=responses_dtype)
        responses_out['responses_out'] = np.array(self.responses_out)
        hdu = fits.BinTableHDU(responses_out, name='RESPONSES_OUT')
        hdul.append(hdu)

        responses_dtype = np.dtype([('responses_map', np.compat.unicode, 200)])
        responses_map = np.zeros(len(self.responses_map), dtype=responses_dtype)
        responses_map['responses_map'] = np.array(self.responses_map)
        hdu = fits.BinTableHDU(responses_map, name='RESPONSES_MAP')
        hdul.append(hdu)

        A = self.Amatrix(self.redshifts)
        hdu = fits.ImageHDU(A, name='A')
        hdul.append(hdu)

        Aout = self.AmatrixOut(self.redshifts)
        hdu = fits.ImageHDU(Aout, name='AOUT')
        hdul.append(hdu)

        hdul.writeto(filename, overwrite=True)
        return

    def fromfits(self, filename=None):
        """Input from FITS file

        Parameters
        ----------

        filename : str
            name of input file

        Notes
        -----

        Expects format of tofits() output.
"""
        self.templates = kcorrect.template.Template(filename=filename)

        hdul = fits.open(filename)

        data = hdul['DATA'].data
        self.nredshift = data['nredshift']
        self.redshift_range = data['redshift_range']
        self.abcorrect = data['abcorrect']

        self.redshifts = hdul['REDSHIFTS'].data
        self.responses = list(hdul['RESPONSES'].data['responses'])
        self.responses_map = list(hdul['RESPONSES_MAP'].data['responses_map'])
        self.responses_out = list(hdul['RESPONSES_OUT'].data['responses_out'])

        f = kcorrect.response.ResponseDict()
        for response in self.responses:
            f.load_response(response)

        A = hdul['A'].data
        Aout = hdul['Aout'].data
        self.Amatrix = self._interpolate_Amatrix(redshifts=self.redshifts, A=A)
        self.AmatrixOut = self._interpolate_Amatrix(redshifts=self.redshifts,
                                                    A=Aout)

        hdul.close()
        return


class KcorrectSDSS(Kcorrect):
    """K-correction object for SDSS data

    Parameters
    ----------

    abcorrect : bool
        correct maggies to AB (default True)

    templates : list of kcorrect.template.SED
        templates to use (if None uses v4 default template set)

    responses : list of str
        names of input responses to base SED on (default to SDSS ugriz)

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

    The 'responses' input defaults to ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0']

    This class provides the method fit_coeffs_asinh() to use SDSS-style
    asinh magnitudes (these are the magnitudes that the SDSS imaging
    reports).

    If abcorrect is True, the to_ab() method is applied to the maggies
    input for absmag() and fit_coeffs() and fit_coeffs_asinh(), which
    adjusts from the SDSS system to the AB system.
"""
    def __init__(self, responses=['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0',
                                  'sdss_z0'], templates=None,
                 responses_out=None, responses_map=None,
                 redshift_range=[0., 2.], nredshift=4000,
                 abcorrect=True, cosmo=None):

        # Initatialize using Kcorrect initialization
        super().__init__(responses=responses, templates=templates,
                         redshift_range=redshift_range,
                         responses_out=responses_out,
                         responses_map=responses_map,
                         nredshift=nredshift, abcorrect=abcorrect,
                         cosmo=cosmo)
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

        ::

            u(AB,2.5m) = u(database, 2.5m) - 0.036
            g(AB,2.5m) = g(database, 2.5m) + 0.012
            r(AB,2.5m) = r(database, 2.5m) + 0.010
            i(AB,2.5m) = i(database, 2.5m) + 0.028
            z(AB,2.5m) = z(database, 2.5m) + 0.040

        fit_coeffs() and absmag() call this on their inputs if abcorrect is True.
"""
        if(ivar is not None):
            maggies, ivar = kcorrect.utils.sdss_ab_correct(maggies=maggies,
                                                           ivar=ivar)
            return(maggies, ivar)
        else:
            maggies = kcorrect.utils.sdss_ab_correct(maggies=maggies,
                                                     ivar=ivar)
            return(maggies)

    def fit_coeffs_asinh(self, redshift=None, mag=None, mag_err=None,
                         extinction=None):
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
        maggies and ivar, and then (if abcorrect is True) calls
        to_ab() method to create AB maggies and ivar.

        If redshift is an array, even with just one element, coeffs is
        returned as an [nredshift, ntemplate] array.

        Otherwise coeffs is returned as an [ntemplate] array.
"""
        if(redshift is None):
            raise TypeError("Must specify redshift to fit coefficients")

        (maggies,
         ivar) = kcorrect.utils.sdss_asinh_to_maggies(mag=mag,
                                                      mag_err=mag_err,
                                                      extinction=extinction)

        coeffs = self.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)
        return(coeffs)


class KcorrectGST(Kcorrect):
    """K-correction object for GALEX, SDSS, 2MASS data

    Parameters
    ----------

    abcorrect : bool
        correct maggies to AB (default True)

    templates : list of kcorrect.template.SED
        templates to use (if None uses v4 default template set)

    responses : list of str
        names of input responses to base SED on (default to FNugrizJHK)

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

    The 'responses' input defaults to ['galex_FUV', 'galex_NUV', 'sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0', 'twomass_J', 'twomass_H', 'twomass_Ks']

    abcorrect is by default False and the input maggies are assumed to
    be AB.  If abcorrect is set to True, the to_ab() method is applied
    to the maggies input for absmag() and fit_coeffs() and
    fit_coeffs_asinh(), which adjusts from the SDSS system to the AB
    system. However, there is no change applied to the 2MASS or
    or GALEX inputs.
"""
    def __init__(self, responses=['galex_FUV', 'galex_NUV', 'sdss_u0', 'sdss_g0',
                                  'sdss_r0', 'sdss_i0', 'sdss_z0', 'twomass_J',
                                  'twomass_H', 'twomass_Ks'],
                 templates=None, responses_out=None, responses_map=None,
                 redshift_range=[0., 2.], nredshift=4000,
                 abcorrect=False, cosmo=None):

        # Initatialize using Kcorrect initialization
        super().__init__(responses=responses, templates=templates,
                         redshift_range=redshift_range,
                         responses_out=responses_out,
                         responses_map=responses_map,
                         nredshift=nredshift, abcorrect=abcorrect,
                         cosmo=cosmo)
        return

    def to_ab(self, maggies=None, ivar=None):
        """Convert FNugrizJHK input maggies to AB

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

        Leaves FN and JHK alone, and fixes ugriz with
        kcorrect.utils.sdss_ab_correct(), which does the
        following:

        Uses the AB conversions produced by D. Eisenstein, in his
        message sdss-calib/1152

        ::

            u(AB,2.5m) = u(database, 2.5m) - 0.036
            g(AB,2.5m) = g(database, 2.5m) + 0.012
            r(AB,2.5m) = r(database, 2.5m) + 0.010
            i(AB,2.5m) = i(database, 2.5m) + 0.028
            z(AB,2.5m) = z(database, 2.5m) + 0.040

        fit_coeffs() and absmag() call this on their inputs if abcorrect is True.
"""
        if(ivar is not None):
            (smaggies,
             sivar) = kcorrect.utils.sdss_ab_correct(maggies=maggies[..., 2:7],
                                                     ivar=ivar[..., 2:7])
            maggies[..., 2:7] = smaggies
            ivar[..., 2:7] = sivar
            return(maggies, ivar)
        else:
            smaggies = kcorrect.utils.sdss_ab_correct(maggies=maggies[..., 2:7],
                                                      ivar=ivar[..., 2:7])
            maggies[..., 2:7] = smaggies
            return(maggies)
