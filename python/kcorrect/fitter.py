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

    abcorrect : bool
        correct maggies to AB (default False)

    responses : list of str
        names of responses to use

    templates : kcorrect.template.Template object
        templates to use

    redshift_range : list of np.float32
        minimum and maximum redshifts (default [0., 1.])

    nredshift : int or np.int32
        number of redshifts in interpolation grid

    Attributes
    ----------

    abcorrect : bool
        correct maggies to AB

    Amatrix : scipy.interpolate.interp1d object
        interpolator to produce A matrix (set to None until set_Amatrix called)

    responses : list of str
        names of responses to use

    nredshift : int or np.int32
        number of redshifts in interpolation grid

    redshift_range : list of np.float32
        minimum and maximum redshifts (default [0., 1.])

    redshifts : ndarray of np.float32
        redshifts in grid

    templates : kcorrect.template.Template object
        templates to use
"""
    def __init__(self, responses=None, templates=None, redshift_range=[0., 1.],
                 nredshift=2000, abcorrect=False):
        self.abcorrect = abcorrect
        self.Amatrix = None
        self.responses = responses
        self.templates = templates
        self.nredshift = nredshift
        self.redshift_range = np.float32(redshift_range)
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

        This method just returns maggies and/or ivar unchanged,
        as for this object we expect AB maggies on input.

        Uses the AB conversions produced by D. Eisenstein, in his
        message sdss-calib/1152

            u(AB,2.5m) = u(database, 2.5m) - 0.036
            g(AB,2.5m) = g(database, 2.5m) + 0.012
            r(AB,2.5m) = r(database, 2.5m) + 0.010
            i(AB,2.5m) = i(database, 2.5m) + 0.028
            z(AB,2.5m) = z(database, 2.5m) + 0.040

        fit_coeffs() calls this on its inputs.
"""
        if(ivar is not None):
            return(maggies, ivar)
        else:
            return(maggies)

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
            fluxes of each band in AB maggies

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

    def _process_inputs(self, redshift=None, maggies=None, ivar=None,
                        coeffs=None):
        """Returns whether input should be an array, and casts everything right

        Parameters
        ----------

        redshift : a quantity or ndarray
            input redshift defining whether we have an array or scalar

        maggies : a quantity or ndarray, or None
            input maggies

        ivar : a quantity or ndarray, or None
            input ivar

        coeffs : a quantity or ndarray, or None
            input coeffs

        Returns
        -------

        array : bool
            True if it is an array

        n : int
            number of redshifts (1 if not an array)

        redshift : np.float32 or [n] ndarray thereof
            redshift to use

        maggies : ndarray of np.float32 or None
            AB maggies to use

        ivar : ndarray of np.float32 or None
            ivar to use

        coeffs : ndarray of np.float32 or None
            coeffs to use
        
        Notes
        -----

        Uses redshift to determine whether we should treat the calculation
        as an array or scalar. Then checks the appropriate sizes (based 
        on the number of responses and the number of seds in the object).
        If maggies, ivar, or coeffs are None, then the corresponding
        output is None.

        Applies this class's to_ab() function on maggies and ivar to
        convert to return AB maggies.
"""
        if(redshift is None):
            raise ValueError("redshift must be defined")

        redshift = np.float32(redshift)

        if(len(redshift.shape) == 0):
            array = False
            n = 1
        elif(len(redshift.shape) == 1):
            array = True
            n = redshift.size
        else:
            raise TypeError("redshift must be 0- or 1-dimensional")

        if(maggies is not None):
            maggies = np.float32(maggies)
            if(array):
                if(len(maggies.shape) != 2):
                    raise ValueError("maggies must be 2-dimensional if redshift is 1-dimensional")
                if(maggies.shape[0] != n):
                    raise ValueError("maggies must have values for each redshift")
                if(maggies.shape[1] != len(self.responses)):
                    raise ValueError("maggies must have values for each band")
            else:
                if(len(maggies.shape) != 1):
                    raise ValueError("maggies must be 1-dimensional if redshift is 0-dimensional")
                if(maggies.shape[0] != len(self.responses)):
                    raise ValueError("maggies must have values for each band")

        if(ivar is not None):
            ivar = np.float32(ivar)
            if(array):
                if(len(ivar.shape) != 2):
                    raise ValueError("ivar must be 2-dimensional if redshift is 1-dimensional")
                if(ivar.shape[0] != n):
                    raise ValueError("ivar must have values for each redshift")
                if(ivar.shape[1] != len(self.responses)):
                    raise ValueError("ivar must have values for each band")
            else:
                if(len(ivar.shape) != 1):
                    raise ValueError("ivar must be 1-dimensional if redshift is 0-dimensional")
                if(ivar.shape[0] != len(self.responses)):
                    raise ValueError("ivar must have values for each band")

        if(coeffs is not None):
            coeffs = np.float32(coeffs)
            if(array):
                if(len(coeffs.shape) != 2):
                    raise ValueError("coeffs must be 2-dimensional if redshift is 1-dimensional")
                if(coeffs.shape[0] != n):
                    raise ValueError("ivar must have values for each redshift")
                if(coeffs.shape[1] != self.templates.nsed):
                    raise ValueError("ivar must have values for each template")
            else:
                if(len(coeffs.shape) != 1):
                    raise ValueError("coeffs must be 1-dimensional if redshift is 0-dimensional")
                if(coeffs.shape[0] != self.templates.nsed):
                    raise ValueError("ivar must have values for each band")

        if(self.abcorrect):
            if(maggies is not None):
                if(ivar is not None):
                    maggies, ivar = self.to_ab(maggies=maggies, ivar=ivar)
                else:
                    maggies = self.to_ab(maggies=maggies)

        return(array, n, redshift, maggies, ivar, coeffs)

    def fit_coeffs(self, redshift=None, maggies=None, ivar=None):
        """Fit coefficients

        Parameters
        ----------

        redshift : np.float32 or ndarray of np.float32
            redshift(s)

        maggies : ndarray of np.float32
            fluxes of each band in AB maggies

        ivar : ndarray of np.float32
            inverse variance of each band

        Returns
        -------

        coeffs : ndarray of np.float32
            coefficients for each template

        Notes
        -----

        maggies are assumed to be Galactic-extinction corrected already.

        Calls this class's to_ab() method on input maggies.

        If redshift is an array, even with just one element, coeffs is
        returned as an [nredshift, ntemplate] array.

        Otherwise coeffs is returned as an [ntemplate] array.
"""
        if(redshift is None):
            raise TypeError("Must specify redshift to fit coefficients")

        # Check a bunch of things about the input arrays
        (array, n, redshift, maggies, ivar,
         dumdum) = self._process_inputs(redshift=redshift, maggies=maggies,
                                        ivar=ivar)

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

    def _reconstruct(self, Amatrix=None, redshift=None, coeffs=None,
                     band_shift=0.):
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

        # Check a bunch of things about the input arrays
        (array, n, redshift, d1, d2,
         coeffs) = self._process_inputs(redshift=redshift, coeffs=coeffs)

        # Consider blueward shift of bandpass due to redshift
        # of observation and due to band_shift
        shift = (1. + redshift) * (1. + band_shift) - 1.

        # Calculate maggies
        try:
            A = self.Amatrix(shift)
        except ValueError:
            return(default_zeros)

        if(array):
            maggies = np.einsum('ijk,ki->ij', A,
                                coeffs.T.reshape(self.templates.nsed, n))
        else:
            maggies = A.dot(coeffs)

        # For band_shift !=0, require this normalization given that
        # AB source is not altered.
        maggies = maggies / (1. + band_shift)

        return(maggies)

    def reconstruct(self, redshift=None, coeffs=None, band_shift=0.):
        """Reconstruct AB maggies associated with coefficients

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
            AB maggies in each band

        Notes
        -----

        Returns AB maggies, but note that if to_ab() is non-trivial,
        these may not be directly comparable to the input maggies.
"""
        return(self._reconstruct(Amatrix=self.Amatrix, redshift=redshift,
                                 coeffs=coeffs, band_shift=band_shift))