#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: utils.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import numpy as np
import kcorrect.template


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
    cspeed = 2.99792   # speed of light in 10^{18} A / s
    ab = 3.631e-2 / wave**2 * cspeed
    return(ab)


def read_basel(filename=None):
    """Read Basel-style spectrum file

    Parameters
    ----------

    filename : str
        file to read from

    Returns
    -------

    info : dict
        meta-data for model

    wave : ndarray of np.float32
        [1221] wavelengths in nm

    flux : ndarray of np.float32
        [N, 1221] fluxes in erg/cm^2/s/Hz

    Notes
    -----

    1221 hardcoded as number of spectral elements.

    Model parameters stored as arrays in info attribute are:
       'modelno'
       'teff'
       'logg'
       'mh'
       'vturb'
       'xh'
"""
    fp = open(filename, "r")
    fpstr = fp.read()
    fp.close

    words = fpstr.split()

    nwave = 1221

    nunits = (len(words) - nwave) // (nwave + 6)

    wavestr = words[0:nwave]
    wave = np.array([np.float32(x) for x in wavestr], dtype=np.float32)

    flux = np.zeros((nunits, nwave), dtype=np.float32)
    modelno = np.zeros(nunits, dtype=np.int32)
    teff = np.zeros(nunits, dtype=np.float32)
    logg = np.zeros(nunits, dtype=np.float32)
    mh = np.zeros(nunits, dtype=np.float32)
    vturb = np.zeros(nunits, dtype=np.float32)
    xh = np.zeros(nunits, dtype=np.float32)

    for unit in range(nunits):
        start = nwave + unit * (nwave + 6)

        modelno[unit] = np.int32(words[start])
        teff[unit] = np.float32(words[start + 1])
        logg[unit] = np.float32(words[start + 2])
        mh[unit] = np.float32(words[start + 3])
        vturb[unit] = np.float32(words[start + 4])
        xh[unit] = np.float32(words[start + 5])

        fluxstr = words[start + 6:start + 6 + nwave]
        flux[unit, :] = np.array([np.float32(x) for x in fluxstr],
                                 dtype=np.float32)

    info = dict()
    info['modelno'] = modelno
    info['teff'] = teff
    info['logg'] = logg
    info['mh'] = mh
    info['vturb'] = vturb
    info['xh'] = xh

    return(info, wave, flux)
