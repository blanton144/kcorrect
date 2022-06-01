#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: utils.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import numpy as np


def sed_ab(wave=None):
    """Calculate AB source spectrum

    Parameters
    ----------

    wave : ndarray of np.float32
        wavelength in Angstrom

    Returns
    -------

    ab : ndarray of np.float32
        AB flux in erg cm^{-2} s^{-1} Angstrom^{-1}
"""
    # This is 3631 Jy * c / wavelength^2
    # Keeping cspeed in order unity units to avoid overflows
    cspeed = 2.99792   # speed of light in 10^{10} cm / s
    ab = 3.631e-10 / wave**2 * cspeed
    return(ab)
