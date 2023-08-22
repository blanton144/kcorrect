import pytest
import os
import re
import numpy as np
import kcorrect
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


def test_load_all_response():
    """Test loading of response into ResponseDict"""
    f = kcorrect.response.ResponseDict()

    rdir = os.path.join(kcorrect.KCORRECT_DIR, 'data', 'responses')
    files = os.listdir(rdir)
    for file in files:
        print(file)
        if(os.path.isfile(os.path.join(rdir, file))):
            m = re.match('^(.*)\\.par$', file)
            if(m is not None):
                response = m.group(1)
                print(response)
                f.load_response(response)

    return


def test_ab_mag_simple():
    """Test that AB source gets magnitude 0"""
    f = kcorrect.response.ResponseDict()
    f.load_response('sdss_u0')

    nwave = 10000
    wave = np.exp(np.log(3000.) + (np.log(6000.) - np.log(3000.)) *
                  (np.arange(nwave, dtype=np.float32) + 0.5) /
                  np.float32(nwave))
    flux = 3631e-23 * 2.99792e+18 / wave**2
    s = kcorrect.template.SED(wave=wave, flux=flux)

    maggies = f['sdss_u0'].project(sed=s)

    assert np.abs(maggies - 1.) < 1.e-7

    return


def test_vega2ab():
    """Test that Vega-to-AB conversion is about right"""
    f = kcorrect.response.ResponseDict()
    f.load_response('sdss_u0')
    f.load_response('twomass_Ks')

    assert hasattr(f['sdss_u0'], 'vega2ab')
    assert np.abs(f['sdss_u0'].vega2ab - 0.93196) < 1.e-4

    assert hasattr(f['twomass_Ks'], 'vega2ab')
    assert np.abs(f['twomass_Ks'].vega2ab - 1.84730) < 1.e-4

    return


def test_solar_magnitudes():
    """Test that solar absolute magnitudes is about right"""
    f = kcorrect.response.ResponseDict()
    f.load_response('sdss_u0')
    f.load_response('twomass_Ks')

    assert hasattr(f['sdss_u0'], 'solar_magnitude')
    assert np.abs(f['sdss_u0'].solar_magnitude - 6.38696) < 1.e-4

    assert hasattr(f['twomass_Ks'], 'solar_magnitude')
    assert np.abs(f['twomass_Ks'].solar_magnitude - 5.13359) < 1.e-4

    return
