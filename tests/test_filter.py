import pytest
import numpy as np
import kcorrect.sed
import kcorrect.filter


def test_init_filter():
    """Test initialization of filter"""
    nwave = 1000
    wave = np.exp(np.log(4000.) + (np.log(6000.) - np.log(4000.)) *
                  (np.arange(nwave, dtype=np.float32) + 0.5) /
                  np.float32(nwave))
    response = np.exp(- 0.5 * (wave - 5000.)**2 / (500.)**2)
    f = kcorrect.filter.Filter(wave=wave, response=response)
    assert f.nwave == nwave
    assert len(f.wave) == nwave
    assert len(f.response) == nwave
    return


def test_load_filter():
    """Test loading of filter into FilterDict"""
    f = kcorrect.filter.FilterDict()

    # Test the load
    f.load_filter('sdss_u0')
    assert(type(f['sdss_u0']) == kcorrect.filter.Filter)

    # Test the singleton nature of the class
    f = 0
    f = kcorrect.filter.FilterDict()
    assert(type(f['sdss_u0']) == kcorrect.filter.Filter)

    return


def test_ab_mag_simple():
    """Test that AB source gets magnitude 0"""
    f = kcorrect.filter.FilterDict()
    f.load_filter('sdss_u0')

    nwave = 10000
    wave = np.exp(np.log(3000.) + (np.log(6000.) - np.log(3000.)) *
                  (np.arange(nwave, dtype=np.float32) + 0.5) /
                  np.float32(nwave))
    flux = 3631e-23 * 2.99792e+10 / wave**2
    s = kcorrect.sed.SED(wave=wave, flux=flux)

    maggies = f['sdss_u0'].project(sed=s)

    assert np.abs(maggies - 1.) < 1.e-7

    return
