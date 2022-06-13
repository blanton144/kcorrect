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
    maggies = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)
    ivar = np.array([1., 1., 1., 1., 1., 1.], dtype=np.float32)

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

    print(coeffs)
    print(k)
    assert k.size == maggies.size

    return
