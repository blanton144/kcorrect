#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: kcorrect.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import os
import numpy as np
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import kcorrect.fitter
import astropy.cosmology


class Kcorrect(kcorrect.fitter.Fitter):
    """K-correction object

    Parameters
    ----------

    responses : list of str
        names of input responses to base SED on

    responses_out : list of str
        output responses for K-corrections (default to "responses")

    responses_map : list of str
        input responses to use for K-corrections (default to "responses")

    nredshift : int or np.int32
        number of redshifts in interpolation grid

    redshift_range : list of np.float32
        minimum and maximum redshifts (default [0., 1.])

    templates : list of kcorrect.template.SED
        templates to use (if None uses v4 default template set)

    Attributes
    ----------

    responses : list of str
        names of input responses to use

    responses_out : list of str
        output responses for K-corrections

    responses_map : list of str
        input responses to use for K-corrections

    nredshift : int or np.int32
        number of redshifts in interpolation grid

    redshift_range : list of np.float32
        minimum and maximum redshifts (default [0., 1.])

    redshifts : ndarray of np.float32
        redshifts in grid

    templates : list of kcorrect.template.SED
        templates to use

    Notes
    -----

    "response_map" defines what the calculated K-correction means. The
    keys are the output magnitude responses "Q" and the values are
    the input magnitude responses "R$, such that the output K-correction
    is K_QR per Blanton & Roweis (2007).

    The default "response_map" is just mapping each item in "responses"
    to itself, i.e. K_QQ.
"""
    def __init__(self, responses=None, responses_out=None, responses_map=None,
                 templates=None, redshift_range=[0., 2.], nredshift=4000):
        if(templates is None):
            filename = os.path.join(os.getenv('KCORRECT_DIR'), 'python',
                                    'kcorrect', 'data', 'templates',
                                    'kcorrect-default-v4.fits')
            templates = kcorrect.template.Template(filename=filename)
        super().__init__(responses=responses, template=templates,
                         redshift_range=redshift_range,
                         nredshift=nredshift)
        self.set_Amatrix()
        if(responses_out is None):
            responses_out = self.responses
        if(responses_map is None):
            responses_map = self.responses
        self.responses_out = responses_out
        self.responses_map = responses_map
        if(self.responses_out == self.responses):
            self.AmatrixOut = self.Amatrix
        else:
            self.AmatrixOut = self._calc_Amatrix(responses=self.responses_out)
        self.imap = np.zeros(len(self.responses_map), dtype=int)
        for i, response in enumerate(self.responses_map):
            self.imap[i] = self.responses.index(response)
        return

    def derived(self, redshift=None, coeffs=None):
        """Return derived quantities based on coefficients

        Parameters
        ----------

        redshift : ndarray of np.float32, or np.float32
            redshift

        coeffs : ndarray of np.float32
            coefficients for each template for each object

        omega0 : np.float32
            z=0 matter density for cosmology

        omegal0 : np.float32
            z=0 cosmological constant for cosmology

        Returns
        -------

        mremain : ndarray of np.float32, or np.float32
           current stellar mass in solar masses

        intsfh : ndarray of np.float32, or np.float32
           current stellar mass in solar masses

        mtol : ndarray of np.float32, or np.float32
           mass-to-light ratio in each band

        b300 : ndarray of np.float32, or np.float32
           current (< 300 Myr) over past star formation

        b1000 : ndarray of np.float32, or np.float32
           current (< 1 Gyr) over past star formation

        metallicity :
           metallicity in current stars
"""
        dfactor = 10.**(0.4 * distmod(redshift))

        intsfh = coeffs.dot(self.templates.intsfh) * dfactor
        mremain = coeffs.dot(self.templates.mremain) * dfactor
        metals =(coeffs.dot(self.templates.mremain * self.templates.mets) *
                 dfactor)

        metallicity = metals / mremain

        m300 = coeffs.dot(self.templates.m300) * dfactor
        b300 = m300 / intsfh
        m1000 = coeffs.dot(self.templates.m1000) * dfactor
        b1000 = m1000 / intsfh

        rmaggies = self.templates.reconstruct(redshift=redshift,
                                              coeffs=coeffs,
                                              band_shift=band_shift)
        mtol = coeffs.dot(self.templates.mremain) / 

        return(mremain, intsfh, mtol, b300, b1000, metallicity)

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
            maggies in each output band
"""
        return(self._reconstruct(Amatrix=self.AmatrixOut, redshift=redshift,
                                 coeffs=coeffs, band_shift=band_shift)

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
        maggies_out = self.reconstruct_out(redshift=redshift, coeffs=coeffs,
                                           band_shift=band_shift)

        kcorrect = maggies_out / maggies[..., imap]
        return(kcorrect)

    def absmag(self, redshift=None, coeffs=None,
               band_shift=0., omega0=None, omegal0=None,
               H0=None):
        return
