import pytest
import os
import numpy as np
import kcorrect.kcorrect


def test_kcorrect_fit_coeffs():
    """Test fitting of SED (does not test quantatively)"""

    responses = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0',
                 'wise_w1']

    kc = kcorrect.kcorrect.Kcorrect(responses=responses,
                                    redshift_range=[0.2, 0.3],
                                    nredshift=100)

    redshift = 0.2532
    maggies = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)
    ivar = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)
    assert coeffs.size == kc.templates.nsed
    assert np.all(coeffs >= 0.)

    redshift = np.array([0.2532])
    maggies = np.array([[1., 1., 1., 1., 1., 1.]], dtype=np.float32)
    ivar = np.array([[1., 1., 1., 1., 1., 1.]], dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)
    assert coeffs.size == kc.templates.nsed
    assert coeffs.shape[0] == 1
    assert coeffs.shape[1] == kc.templates.nsed
    assert np.all(coeffs >= 0.)

    redshift = np.array([0.2532, 0.2532])
    maggies = np.array([[1., 1., 1., 1., 1., 1.],
                        [2., 2., 2., 2., 2., 2.]], dtype=np.float32)
    ivar = np.array([[1., 1., 1., 1., 1., 1.],
                     [1., 1., 1., 1., 1., 1.]], dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    assert coeffs.size == kc.templates.nsed * 2
    assert coeffs.shape[0] == 2
    assert coeffs.shape[1] == kc.templates.nsed
    assert np.all(coeffs >= 0.)

    return


def test_kcorrect():
    """Test K-correction calculation (doesn't test quantitatively"""

    responses = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0',
                 'wise_w1']

    kc = kcorrect.kcorrect.Kcorrect(responses=responses,
                                    redshift_range=[0., 0.3],
                                    nredshift=100)

    redshift = 0.2532
    maggies = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)
    ivar = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    k = kc.kcorrect(redshift=redshift, coeffs=coeffs)

    assert k.size == maggies.size

    redshift = np.array([0.1, 0.2532], dtype=np.float32)
    maggies = np.ones((2, len(responses)), dtype=np.float32)
    ivar = np.ones((2, len(responses)), dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    k = kc.kcorrect(redshift=redshift, coeffs=coeffs)

    assert k.size == maggies.size
    assert k.shape[0] == maggies.shape[0]
    assert k.shape[1] == len(responses)

    return


def test_derived():
    """Test derived quantity calculation (doesn't test quantitatively)"""

    responses = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0',
                 'wise_w1']

    kc = kcorrect.kcorrect.Kcorrect(responses=responses,
                                    redshift_range=[0., 0.3],
                                    nredshift=100)

    redshift = 0.2532
    maggies = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)
    ivar = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    d = kc.derived(redshift=redshift, coeffs=coeffs)

    assert len(list(d.keys())) == 7
    for cd in d:
        if(cd == 'mtol'):
            assert d[cd].size == len(maggies)
        else:
            assert d[cd].size == 1

    redshift = np.array([0.1, 0.2532], dtype=np.float32)
    maggies = np.ones((2, len(responses)), dtype=np.float32)
    ivar = np.ones((2, len(responses)), dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    d = kc.derived(redshift=redshift, coeffs=coeffs)

    assert len(list(d.keys())) == 7
    for cd in d:
        if(cd == 'mtol'):
            assert d[cd].size == maggies.size
            assert d[cd].shape == maggies.shape
        else:
            assert d[cd].size == redshift.size
            assert d[cd].shape == redshift.shape

    return


def test_absmag():
    """Test absolute magnitude quantity calculation (doesn't test quantitatively)"""

    responses = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0',
                 'wise_w1']

    kc = kcorrect.kcorrect.Kcorrect(responses=responses,
                                    redshift_range=[0., 0.3],
                                    nredshift=100)

    redshift = 0.2532
    maggies = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)
    ivar = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    absmag = kc.absmag(redshift=redshift, maggies=maggies,
                       ivar=ivar, coeffs=coeffs)

    assert absmag.size == maggies.size
    assert absmag.shape == maggies.shape

    redshift = np.array([0.2532, 0.11], dtype=np.float32)
    maggies = np.ones((2, len(responses)), dtype=np.float32)
    ivar = np.ones((2, len(responses)), dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    absmag = kc.absmag(redshift=redshift, maggies=maggies,
                       ivar=ivar, coeffs=coeffs)

    assert absmag.size == maggies.size
    assert absmag.shape == maggies.shape

    return


def test_absmag_reconstruct():
    """Test absolute magnitude quantity calculation (doesn't test quantitatively)"""

    responses = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0',
                 'wise_w1']

    kc = kcorrect.kcorrect.Kcorrect(responses=responses,
                                    redshift_range=[0., 0.3],
                                    nredshift=100)

    redshift = 0.2532
    maggies = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)
    ivar = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    absmag, absmag_reconstruct = kc.absmag(redshift=redshift, maggies=maggies,
                                           ivar=ivar, coeffs=coeffs, reconstruct=True)

    assert absmag.size == maggies.size
    assert absmag.shape == maggies.shape
    assert absmag_reconstruct.size == maggies.size
    assert absmag_reconstruct.shape == maggies.shape

    redshift = np.array([0.2532, 0.11], dtype=np.float32)
    maggies = np.ones((2, len(responses)), dtype=np.float32)
    ivar = np.ones((2, len(responses)), dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    absmag, absmag_reconstruct = kc.absmag(redshift=redshift, maggies=maggies,
                                           ivar=ivar, coeffs=coeffs, reconstruct=True)

    assert absmag.size == maggies.size
    assert absmag.shape == maggies.shape
    assert absmag_reconstruct.size == maggies.size
    assert absmag_reconstruct.shape == maggies.shape

    return


def test_absmag_limit():
    """Test absolute magnitude quantity calculation (doesn't test quantitatively)"""

    responses = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0',
                 'wise_w1']

    kc = kcorrect.kcorrect.Kcorrect(responses=responses,
                                    redshift_range=[0., 0.3],
                                    nredshift=100)

    redshift = 0.2532
    maggies = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)
    ivar = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    absmag, absmag_limit = kc.absmag(redshift=redshift, maggies=maggies,
                                     ivar=ivar, coeffs=coeffs, limit=True)

    assert absmag.size == maggies.size
    assert absmag.shape == maggies.shape
    assert absmag_limit.size == maggies.size
    assert absmag_limit.shape == maggies.shape

    redshift = np.array([0.2532, 0.11], dtype=np.float32)
    maggies = np.ones((2, len(responses)), dtype=np.float32)
    ivar = np.ones((2, len(responses)), dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    absmag, absmag_limit = kc.absmag(redshift=redshift, maggies=maggies,
                                     ivar=ivar, coeffs=coeffs, limit=True)

    assert absmag.size == maggies.size
    assert absmag.shape == maggies.shape
    assert absmag_limit.size == maggies.size
    assert absmag_limit.shape == maggies.shape

    return


def test_absmag_reconstruct_limit():
    """Test absolute magnitude quantity calculation (doesn't test quantitatively)"""

    responses = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0',
                 'wise_w1']

    kc = kcorrect.kcorrect.Kcorrect(responses=responses,
                                    redshift_range=[0., 0.3],
                                    nredshift=100)

    redshift = 0.2532
    maggies = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)
    ivar = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    absmag, absmag_reconstruct, absmag_limit = kc.absmag(redshift=redshift, maggies=maggies,
                                                         ivar=ivar, coeffs=coeffs, limit=True,
                                                         reconstruct=True)

    assert absmag.size == maggies.size
    assert absmag.shape == maggies.shape
    assert absmag_limit.size == maggies.size
    assert absmag_limit.shape == maggies.shape
    assert absmag_reconstruct.size == maggies.size
    assert absmag_reconstruct.shape == maggies.shape

    redshift = np.array([0.2532, 0.11], dtype=np.float32)
    maggies = np.ones((2, len(responses)), dtype=np.float32)
    ivar = np.ones((2, len(responses)), dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    absmag, absmag_reconstruct, absmag_limit = kc.absmag(redshift=redshift, maggies=maggies,
                                                         ivar=ivar, coeffs=coeffs, limit=True,
                                                         reconstruct=True)

    assert absmag.size == maggies.size
    assert absmag.shape == maggies.shape
    assert absmag_limit.size == maggies.size
    assert absmag_limit.shape == maggies.shape
    assert absmag_reconstruct.size == maggies.size
    assert absmag_reconstruct.shape == maggies.shape

    return


def test_sdss_kcorrect():
    """Test K-correction calculation in SDSS (doesn't test quantitatively"""

    kc = kcorrect.kcorrect.KcorrectSDSS(redshift_range=[0., 0.3],
                                        nredshift=100)

    redshift = 0.2532
    maggies = np.array([1., 1., 1., 1., 1.], dtype=np.float32)
    ivar = np.array([1., 1., 1., 1., 1.], dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    k = kc.kcorrect(redshift=redshift, coeffs=coeffs)

    assert k.size == maggies.size

    redshift = np.array([0.1, 0.2532], dtype=np.float32)
    maggies = np.ones((2, 5), dtype=np.float32)
    ivar = np.ones((2, 5), dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    k = kc.kcorrect(redshift=redshift, coeffs=coeffs)

    assert kc.responses == ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0']
    assert k.size == maggies.size
    assert k.shape[0] == maggies.shape[0]
    assert k.shape[1] == 5

    return


def test_sdss_kcorrect_asinh():
    """Test K-correction calculation in SDSS (doesn't test quantitatively"""

    kc = kcorrect.kcorrect.KcorrectSDSS(redshift_range=[0., 0.3],
                                        nredshift=100)

    redshift = 0.2532
    mag = np.array([1., 1., 1., 1., 1.], dtype=np.float32)
    err = np.array([1., 1., 1., 1., 1.], dtype=np.float32)
    extinction = np.array([1., 1., 1., 1., 1.], dtype=np.float32)

    coeffs = kc.fit_coeffs_asinh(redshift=redshift, mag=mag, mag_err=err,
                                 extinction=extinction)

    k = kc.kcorrect(redshift=redshift, coeffs=coeffs)

    assert k.size == mag.size

    redshift = np.array([0.1, 0.2532], dtype=np.float32)
    mag = np.ones((2, 5), dtype=np.float32)
    err = np.ones((2, 5), dtype=np.float32)
    extinction = np.ones((2, 5), dtype=np.float32)

    coeffs = kc.fit_coeffs_asinh(redshift=redshift, mag=mag, mag_err=err,
                                 extinction=extinction)

    k = kc.kcorrect(redshift=redshift, coeffs=coeffs)

    assert kc.responses == ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0',
                            'sdss_z0']
    assert k.size == mag.size
    assert k.shape[0] == mag.shape[0]
    assert k.shape[1] == 5

    return


def test_gst_kcorrect():
    """Test K-correction calculation in GALEX-SDSS-2MASS (doesn't test quantitatively)"""

    kc = kcorrect.kcorrect.KcorrectGST(redshift_range=[0., 0.3],
                                       nredshift=100)

    redshift = 0.2532
    maggies = np.ones(10, dtype=np.float32)
    ivar = np.ones(10, dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    k = kc.kcorrect(redshift=redshift, coeffs=coeffs)

    assert k.size == maggies.size

    redshift = np.array([0.1, 0.2532], dtype=np.float32)
    maggies = np.ones((2, 10), dtype=np.float32)
    ivar = np.ones((2, 10), dtype=np.float32)

    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

    k = kc.kcorrect(redshift=redshift, coeffs=coeffs)

    assert kc.responses == ['galex_FUV', 'galex_NUV', 'sdss_u0', 'sdss_g0',
                            'sdss_r0', 'sdss_i0', 'sdss_z0', 'twomass_J',
                            'twomass_H', 'twomass_Ks']
    assert k.size == maggies.size
    assert k.shape[0] == maggies.shape[0]
    assert k.shape[1] == 10

    return
