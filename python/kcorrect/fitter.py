#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: fit.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import numpy as np
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import kcorrect.filter


class Fitter(object):
    """Nonnegative SED fitting object

    Parameters
    ----------

    filters : list of str
        names of filters to use

    nredshift : int or np.int32
        number of redshifts in interpolation grid

    redshift_range : list of np.float32
        minimum and maximum redshifts (default [0., 1.])

    templates : list of kcorrect.sed.SED
        templates to use 

    Attributes
    ----------

    filters : list of str
        names of filters to use

    nredshift : int or np.int32
        number of redshifts in interpolation grid

    redshift_range : list of np.float32
        minimum and maximum redshifts (default [0., 1.])
    
    redshifts : ndarray of np.float32
        redshifts in grid

    templates : list of kcorrect.sed.SED
        templates to use 

    rmatrix : ndarray of np.float32
        [nredshift, ntemplates, nfilters] matrix of projections
"""
    def __init__(self, filters=None, templates=None, redshift_range=[0., 1.],
                 nredshift=2000):
        self.filters = filters
        self.templates = templates
        self.nredshift = nredshift
        self.redshift_range = redshift_range
        self.redshifts = (self.redshift_range[0] +
                          (self.redshift_range[1] - self.redshift_range[0]) *
                          np.arange(nredshift, dtype=np.float32) /
                          np.float32(nredshift - 1))
        return

    def set_rmatrix(self):
        """Set rmatrix, the matrix of projections"""
        # Make sure filters are loaded
        f = kcorrect.filter.FilterDict()
        for filt in self.filters:
            f.load_filter(filt)

        # Create rmatrix
        self.rmatrix = np.zeros((self.nredshift,
                                 len(self.filters),
                                 self.templates.nsed), dtype=np.float32)
        for iz, z in enumerate(self.redshifts):
            self.templates.set_redshift(redshift=z)
            for ifilt, filt in enumerate(self.filters):
                self.rmatrix[iz, ifilt, :] = f[filt].project(sed=self.templates)

        # Now create Amatrix interpolator
        self.Amatrix = interpolate.interp1d(self.redshifts, self.rmatrix,
                                            axis=0)

        return

    def fit_coeffs(self, redshift=None, maggies=None, ivar=None):
        """Fit coefficients

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
