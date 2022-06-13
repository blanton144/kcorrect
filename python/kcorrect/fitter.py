#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: fitter.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import numpy as np
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import kcorrect.response


class Fitter(object):
    """Nonnegative SED fitting object

    Parameters
    ----------

    responses : list of str
        names of responses to use

    nredshift : int or np.int32
        number of redshifts in interpolation grid

    redshift_range : list of np.float32
        minimum and maximum redshifts (default [0., 1.])

    templates : list of kcorrect.template.SED
        templates to use

    Attributes
    ----------

    Amatrix : scipy.interpolate.interp1d object
        interpolator to produce A matrix

    responses : list of str
        names of responses to use

    nredshift : int or np.int32
        number of redshifts in interpolation grid

    redshift_range : list of np.float32
        minimum and maximum redshifts (default [0., 1.])

    redshifts : ndarray of np.float32
        redshifts in grid

    templates : list of kcorrect.template.SED
        templates to use
"""
    def __init__(self, responses=None, templates=None, redshift_range=[0., 1.],
                 nredshift=2000):
        self.responses = responses
        self.templates = templates
        self.nredshift = nredshift
        self.redshift_range = redshift_range
        self.redshifts = (self.redshift_range[0] +
                          (self.redshift_range[1] - self.redshift_range[0]) *
                          np.arange(nredshift, dtype=np.float32) /
                          np.float32(nredshift - 1))
        return

    def _calc_Amatrix(self, responses=None):
        """Create an A matrix and return it (don't set attribute)"""
        # Make sure responses are loaded
        f = kcorrect.response.ResponseDict()
        for response in responses:
            f.load_response(response)

        # Create rmatrix
        rmatrix = np.zeros((self.nredshift,
                            len(responses),
                            self.templates.nsed), dtype=np.float32)
        for iz, z in enumerate(self.redshifts):
            self.templates.set_redshift(redshift=z)
            for ir, response in enumerate(responses):
                rmatrix[iz, ir, :] = f[response].project(sed=self.templates)

        # Now create Amatrix interpolator
        Amatrix = interpolate.interp1d(self.redshifts, rmatrix, axis=0)

        return(Amatrix)

    def set_Amatrix(self):
        """Set Amatrix, interpolator for the design matrix"""
        self.Amatrix = self._calc_Amatrix(responses=self.responses)

        return

    def _fit_coeffs(self, redshift=None, maggies=None, ivar=None):
        """Fit coefficients to single object

        Parameters
        ----------

        redshift : np.float32
            redshift

        maggies : ndarray of np.float32
            fluxes of each band in maggies

        ivar : ndarray of np.float32
            inverse variance of each band

        Returns
        -------

        coeffs : ndarray of np.float32
            coefficients for each template
"""
        default_zeros = np.zeros(self.templates.nsed, dtype=np.float32)

        inverr = np.sqrt(ivar)
        inverr_matrix = np.transpose(np.tile(inverr, (self.templates.nsed, 1)))
        try:
            A = inverr_matrix * self.Amatrix(redshift)
        except ValueError:
            return(default_zeros)
        b = maggies * inverr

        coeffs, rnorm = optimize.nnls(A, b)

        return(coeffs)

    def fit_coeffs(self, redshift=None, maggies=None, ivar=None):
        """Fit coefficients

        Parameters
        ----------

        redshift : np.float32 or ndarray of np.float32
            redshift(s)

        maggies : ndarray of np.float32
            fluxes of each band in maggies

        ivar : ndarray of np.float32
            inverse variance of each band

        Returns
        -------

        coeffs : ndarray of np.float32
            coefficients for each template

        Notes
        -----

        If redshift is an array, even with just one element, coeffs is
        returned as an [nredshift, ntemplate] array.

        Otherwise coeffs is returned as an [ntemplate] array.
"""
        if(redshift is None):
            raise TypeError("Must specify redshift to fit coefficients")

        # Check a bunch of things about the input arrays
        try:
            n = len(redshift)
            array = True
        except TypeError:
            n = 1
            array = False
        if(n == 1):
            nr = maggies.size
        else:
            if(len(maggies.shape) != 2):
                raise IndexError("maggies not 2D")
            nr = maggies.shape[1]
            nz = maggies.shape[0]
            if(nz != n):
                raise IndexError("maggies and redshift differ in number objects")

        if(nr != len(self.responses)):
            raise IndexError("Number of maggies differs from number of responses")
        if(maggies.size != ivar.size):
            raise IndexError("Number of maggies differs from number of ivar")
        if(maggies.shape != ivar.shape):
            raise IndexError("maggies and ivar differ in shape")

        # Call single
        if(n == 1):
            coeffs = self._fit_coeffs(redshift=np.squeeze(redshift),
                                      maggies=np.squeeze(maggies),
                                      ivar=np.squeeze(ivar))
            if(array):
                coeffs = coeffs.reshape(1, len(coeffs))
            return(coeffs)

        # Loop for multiple
        coeffs = np.zeros((n, self.templates.nsed), dtype=np.float32)
        for i, r in enumerate(redshift):
            coeffs[i, :] = self._fit_coeffs(redshift=r, maggies=maggies[i, :],
                                            ivar=ivar[i, :])

        return(coeffs)

    def _reconstruct(self, Amatrix=None, redshift=None, coeffs=None, band_shift=0.):
        """Reconstruct maggies associated with coefficients

        Parameters
        ----------

        Amatrix : scipy.interpolate.interp1d
            interpolator to use for Amatri

        redshift : np.float32
            redshift

        coeffs : ndarray of np.float32
            coefficients

        band_shift : np.float32
            blueshift to apply to reconstructed bandpasses

        Returns
        -------

        maggies : ndarray of np.float32
            maggies in each band
"""
        default_zeros = np.zeros(len(self.responses), dtype=np.float32)

        # Consider blueward shift of bandpass due to redshift
        # of observation and due to band_shift
        shift = (1. + redshift) * (1. + band_shift) - 1.

        # Calculate maggies
        try:
            A = self.Amatrix(shift)
        except ValueError:
            return(default_zeros)
        maggies = A.dot(coeffs)

        # For band_shift !=0, require this normalization given that
        # AB source is not altered.
        maggies = maggies / (1. + band_shift)

        return(maggies)

    def reconstruct(self, redshift=None, coeffs=None, band_shift=0.):
        """Reconstruct maggies associated with coefficients

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
            maggies in each band
"""
        return(self._reconstruct(Amatrix=self.Amatrix, redshift=redshift,
                                 coeffs=coeffs, band_shift=band_shift))
