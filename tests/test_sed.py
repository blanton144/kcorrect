import pytest
import numpy as np
import kcorrect.sed


def test_init_sed():
    """Test initialization of SED"""
    nwave = 10000
    wave = np.exp(np.log(1000.) + (np.log(1.e+6) - np.log(1000.)) *
                  (np.arange(nwave, dtype=np.float32) + 0.5) /
                  np.float32(nwave))
    flux = np.ones(nwave, dtype=np.float32)
    s = kcorrect.sed.SED(wave=wave, flux=flux)
    assert s.nwave == nwave
    assert len(s.wave) == nwave
    assert len(s.flux) == nwave
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
    s = kcorrect.sed.SED(wave=wave, flux=flux)

    s.tofits('tmp-sed-write-and-read.fits')

    t = kcorrect.sed.SED(filename='tmp-sed-write-and-read.fits')
    assert t.nwave == nwave
    assert len(t.wave) == nwave
    assert len(t.flux) == nwave
    assert (t.wave == wave).min() == True
    assert (t.flux == flux).min() == True
    assert (t.restframe_wave == wave).min() == True
    assert (t.restframe_flux == flux).min() == True
    return


def test_redshift():
    """Test redshifting SED"""
    nwave = 10000
    wave = np.exp(np.log(1000.) + (np.log(1.e+6) - np.log(1000.)) *
                  (np.arange(nwave, dtype=np.float32) + 0.5) /
                  np.float32(nwave))
    flux = np.ones(nwave, dtype=np.float32) + np.log10(wave)
    s = kcorrect.sed.SED(wave=wave, flux=flux)

    z = 0.4
    s.set_redshift(redshift=z)
    assert (s.wave == s.restframe_wave * (1. + z)).min() == True
    assert (s.flux == s.restframe_flux / (1. + z)).min() == True

    return
