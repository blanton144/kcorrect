#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Filename: utils.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import numpy as np

import kcorrect.response


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


def sdss_ab_correct(maggies=None, ivar=None,
                    abfix=[-0.036, 0.012, 0.010, 0.028, 0.040]):
    """Return AB maggies based on SDSS standard maggies

    Parameters
    ----------

    maggies : ndarray of np.float32
        array of fluxes in standard SDSS system

    ivar : ndarray of np.float32
        inverse variances in standard SDSS system (optional)

    abfix : ndarray or list of np.float32
        [5]-array of magnitude offset to apply to standard values

    Returns
    -------

    ab_maggies : ndarray of np.float32
        array of fluxes converted to AB

    ab_ivar : ndarray of np.float32
        inverse variances converted to AB (if ivar input)

    Notes
    -----

    Uses the AB conversions produced by D. Eisenstein, in his
    message sdss-calib/1152

    ::

      u(AB,2.5m) = u(database, 2.5m) - 0.036
      g(AB,2.5m) = g(database, 2.5m) + 0.012
      r(AB,2.5m) = r(database, 2.5m) + 0.010
      i(AB,2.5m) = i(database, 2.5m) + 0.028
      z(AB,2.5m) = z(database, 2.5m) + 0.040

    Unless an alternative is specified with the abfix parameter.
"""
    maggies = np.float32(maggies)
    if(ivar is not None):
        ivar = np.float32(ivar)

        if(maggies.shape != ivar.shape):
            raise("maggies and ivar must be the same shape")

    abfix = np.array(abfix, dtype=np.float32)
    if(abfix.shape != (5,)):
        raise("abfix must have 5 values")
    abfac = 10.**(- 0.4 * abfix)

    if(maggies.ndim == 1):
        if(maggies.size != 5):
            raise("sdss_ab_correct expects 5 SDSS maggies values (ugriz)")
        ab_maggies = maggies * abfac
        if(ivar is not None):
            ab_ivar = ivar / abfac**2
    else:
        ab_maggies = np.zeros(maggies.shape, dtype=np.float32)
        for i, cabfac in enumerate(abfac):
            ab_maggies[..., i] = maggies[..., i] * cabfac
        if(ivar is not None):
            ab_ivar = np.zeros(maggies.shape, dtype=np.float32)
            for i, cabfac in enumerate(abfac):
                ab_ivar[..., i] = ivar[..., i] / cabfac**2

    if(ivar is not None):
        return(ab_maggies, ab_ivar)
    else:
        return(ab_maggies)


def error_floor(floor=None, maggies=None, ivar=None):
    """Calculate new inverse variance with fractional error floors

    Parameters
    ----------

    floor : ndarray of np.float32
        [Nr] fractional error floor for each response

    maggies : ndarray of np.float32
        [N, Nr] or [Nr] array of maggies

    ivar : ndarray of np.float32
        [N, Nr] or [Nr] array of inverse variance of maggies

    Returns
    -------

    ivar : ndarray of np.float32
        [N, Nr] or [Nr] array of inverse variances
"""
    floor = np.array(floor)
    maggies = np.array(maggies)
    ivar = np.array(ivar)

    if(maggies.ndim == 1):
        array = False
    else:
        array = True

    if(maggies.shape[-1] != floor.shape[0]):
        raise("floor must have same number of responses as maggies")

    if(maggies.shape != ivar.shape):
        raise("maggies and ivar must be same shape")

    if(array):
        for i in np.arange(floor.size):
            iok = np.where(ivar[:, i] > 0.)[0]
            if(len(iok) > 0):
                ferr = 1. / (np.abs(maggies[iok, i]) * np.sqrt(ivar[iok, i]))
                ioklow = iok[np.where(ferr < floor[i])[0]]
                ivar[ioklow, i] = 1. / (np.abs(maggies[ioklow, i]) *
                                        floor[i])**2
    else:
        for i in np.arange(floor.size):
            if(ivar[i] > 0.):
                ferr = 1. / (np.abs(maggies[i]) * np.sqrt(ivar[i]))
                if(ferr < floor[i]):
                    ivar[i] = 1. / (np.abs(maggies[i]) * floor[i])**2

    return(ivar)


def sdss_asinh_to_maggies(mag=None, mag_err=None, extinction=None):
    """Calculate maggies and ivar based on catalog parameters

    Parameters
    ----------

    mag : ndarray of np.float32
        [N, 5] or [5] array of asinh magnitudes from SDSS

    mag_err : ndarray of np.float32
        [N, 5] or [5] array of asinh magnitude errors from SDSS

    extinction : ndarray or list of np.float32
        [N, 5] or [5] array of Galactic extinctions from SDSS

    Returns
    -------

    maggies : ndarray of np.float32
        [N, 5] or [5] array of maggies

    ivar : ndarray of np.float32
        [N, 5] or [5] array of inverse variances

    Notes
    -----

    If extinction set, applies extinction (after converting to maggies). 

    Does not apply AB corrections.

    If mag_err is None on input, only maggies is returned.
"""
    mag = np.float32(mag)

    if(mag_err is not None):
        mag_err = np.float32(mag_err)
        if(mag.shape != mag_err.shape):
            raise("mag and mag_err must be the same shape")

    if(extinction is not None):
        extinction = np.float32(extinction)
        if(mag.shape != extinction.shape):
            raise("mag and extinction must be the same shape")
    else:
        extinction = np.zeros(mag.shape, dtype=np.float32)

    if(mag.shape[-1] != 5):
        raise("mag must have 5 values per object")

    b0 = np.array([1.4e-10, 0.9e-10, 1.2e-10, 1.8e-10, 7.4e-10],
                  dtype=np.float32)

    maggies = np.zeros(mag.shape, dtype=np.float32)
    for k in np.arange(5, dtype=int):
        maggies[..., k] = 2. * b0[k] * np.sinh(- np.log(b0[k])
                                               - (0.4 * np.log(10.)
                                                  * mag[..., k]))
        maggies[..., k] = (maggies[..., k] * 10.**(0.4 * extinction[..., k]))

    if(mag_err is not None):
        maggies_ivar = np.zeros(mag.shape, dtype=np.float32)
        for k in np.arange(5, dtype=int):
            maggies_err = (2. * b0[k] * np.cosh(- np.log(b0[k])
                                                - (0.4 * np.log(10.) *
                                                   mag[..., k])) *
                           0.4 * np.log10(10.) * mag_err[..., k])
            maggies_ivar[..., k] = 1. / maggies_err**2

    if(mag_err is None):
        return(maggies)
    else:
        return(maggies, maggies_ivar)


def maggies2flambda(maggies=None, ivar=None, responses=None):
    """Return arrays of wave, flambda for maggies

    Parameters
    ----------

    maggies : np.float32 or ndarray of np.float32
        maggies to convert

    ivar : np.float32 or ndarray of np.float32
        inverse variance of maggies to convert (or None)

    responses : list of str
        names of responses that each maggies corresponds to

    Returns
    -------

    wave : np.float32 or ndarray of np.float32
        effective wavelength of each response curve

    flux : np.float32 or ndarray of np.float32
        flux in erg/s/cm^2/A corresponding to maggies

    flux_ivar : np.float32 or ndarray of np.float32
        inverse variance of flux (if ivar input)

    Notes
    -----

    If ivar is set on input, output is (wave, flux, flux_ivar).

    If ivar is not set on input, output is (wave, flux)

    Effective wavelength is derived from response curve; curve is
    loaded into the response dictionary if it isn't there already.
"""
    
    r = kcorrect.response.ResponseDict()
    for response in responses:
        if(response not in r):
            r.load_response(response)

    wave = np.array([r[x].lambda_eff for x in responses],
                    dtype=np.float32)
    m2f = sed_ab(wave)

    if(maggies.ndim == 1):
        flambda = m2f * maggies
        if(ivar is not None):
            fivar = ivar / m2f**2
    else:
        flambda = np.outer(np.ones(maggies.shape[0], dtype=np.float32),
                           m2f) * maggies
        if(ivar is not None):
            fivar = np.outer(np.ones(maggies.shape[0], dtype=np.float32),
                             1. / m2f**2) * ivar

    if(ivar is not None):
        return(wave, flambda, fivar)
    else:
        return(wave, flambda)
