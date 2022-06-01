import pytest
import numpy as np
import kcorrect.template
import kcorrect.response


def test_init_response():
    """Test initialization of response"""
    nwave = 1000
    wave = np.exp(np.log(4000.) + (np.log(6000.) - np.log(4000.)) *
                  (np.arange(nwave, dtype=np.float32) + 0.5) /
                  np.float32(nwave))
    response = np.exp(- 0.5 * (wave - 5000.)**2 / (500.)**2)
    f = kcorrect.response.Response(wave=wave, response=response)
    assert f.nwave == nwave
    assert len(f.wave) == nwave
    assert len(f.response) == nwave
    return


def test_load_response():
    """Test loading of response into ResponseDict"""
    f = kcorrect.response.ResponseDict()

    # Test the load
    f.load_response('sdss_u0')
    assert(type(f['sdss_u0']) == kcorrect.response.Response)

    # Test the singleton nature of the class
    f = 0
    f = kcorrect.response.ResponseDict()
    assert(type(f['sdss_u0']) == kcorrect.response.Response)

    return


def test_ab_mag_simple():
    """Test that AB source gets magnitude 0"""
    f = kcorrect.response.ResponseDict()
    f.load_response('sdss_u0')

    nwave = 10000
    wave = np.exp(np.log(3000.) + (np.log(6000.) - np.log(3000.)) *
                  (np.arange(nwave, dtype=np.float32) + 0.5) /
                  np.float32(nwave))
    flux = 3631e-23 * 2.99792e+10 / wave**2
    s = kcorrect.template.SED(wave=wave, flux=flux)

    maggies = f['sdss_u0'].project(sed=s)

    assert np.abs(maggies - 1.) < 1.e-7

    return
