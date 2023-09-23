import os

import numpy as np
import pytest

import kcorrect.template


def test_init_sed():
    """Test initialization of SED"""
    nwave = 10000
    wave = np.exp(np.log(1000.) + (np.log(1.e+6) - np.log(1000.)) *
                  (np.arange(nwave, dtype=np.float32) + 0.5) /
                  np.float32(nwave))
    flux = np.ones(nwave, dtype=np.float32)
    s = kcorrect.template.SED(wave=wave, flux=flux)
    assert s.nwave == nwave
    assert len(s.wave) == nwave
    assert len(s.flux[0, :]) == nwave
    assert (s.wave == s.restframe_wave).min() == True
    assert (s.flux == s.restframe_flux).min() == True
    return


def test_write_and_read_sed():
    """Test reading and writing of SED"""
    nwave = 10000
    wave = np.exp(np.log(1000.) + (np.log(1.e+6) - np.log(1000.)) *
                  (np.arange(nwave, dtype=np.float32) + 0.5) /
                  np.float32(nwave))
    flux = np.ones(nwave, dtype=np.float32)
    s = kcorrect.template.SED(wave=wave, flux=flux)

    s.tofits('tmp-sed-write-and-read.fits')

    t = kcorrect.template.SED(filename='tmp-sed-write-and-read.fits')
    assert t.nwave == nwave
    assert len(t.wave) == nwave
    assert len(t.flux[0, :]) == nwave
    assert (t.wave == wave).min() == True
    assert (t.flux == flux).min() == True
    assert (t.restframe_wave == wave).min() == True
    assert (t.restframe_flux == flux).min() == True

    os.remove('tmp-sed-write-and-read.fits')
    return


def test_redshift():
    """Test redshifting SED"""
    nwave = 10000
    wave = np.exp(np.log(1000.) + (np.log(1.e+6) - np.log(1000.)) *
                  (np.arange(nwave, dtype=np.float32) + 0.5) /
                  np.float32(nwave))
    flux = np.ones(nwave, dtype=np.float32) + np.log10(wave)
    s = kcorrect.template.SED(wave=wave, flux=flux)

    z = 0.4
    s.set_redshift(redshift=z)
    assert (s.wave == s.restframe_wave * (1. + z)).min() == True
    assert (s.flux == s.restframe_flux / (1. + z)).min() == True

    return


def test_project():
    """Test projection onto filters (does not test quantitatively)"""

    filename = os.path.join(kcorrect.KCORRECT_DIR,
                            'data', 'templates',
                            'kcorrect-lrg1-v4.fits')
    s = kcorrect.template.SED(filename=filename)
    f = kcorrect.response.ResponseDict()
    responses = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0']
    for response in responses:
        f.load_response(response)

    z = 0.4
    s.set_redshift(redshift=z)
    for j, response in enumerate(responses):
        maggies = f[response].project(sed=s)
        assert maggies > 0.
