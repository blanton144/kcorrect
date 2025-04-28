import pytest
import os
import numpy as np
import kcorrect
import kcorrect.template
import kcorrect.fitter


def test_fitter():
    """Test fitting of SED (does not test quantatively)"""

    filename = os.path.join(kcorrect.KCORRECT_DIR, 'data',
                            'templates',
                            'kcorrect-default-v4.fits')
    templates = kcorrect.template.SED(filename=filename)

    responses = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0',
                 'twomass_Ks', 'wise_w4']

    fitter = kcorrect.fitter.Fitter(templates=templates,
                                    responses=responses,
                                    redshift_range=[0., 0.1],
                                    nredshift=100)
    fitter.set_Amatrix()

    redshift = 0.0532
    maggies = np.array([1., 1., 1., 1., 1., 1., 1.], dtype=np.float32)
    ivar = np.array([1., 1., 1., 1., 1., 1., 1.], dtype=np.float32)

    coeffs = fitter.fit_coeffs(maggies=maggies, ivar=ivar, redshift=redshift)

    fitspec = coeffs.dot(templates.flux)
    fitmaggies = fitter.reconstruct(redshift=redshift, coeffs=coeffs)

    assert np.all(fitspec >= 0.)
    assert np.all(fitmaggies >= 0.)
    return


def test_fitter_nointerp():
    """Test that fitting of SED is consistent with and without interpolation """

    filename = os.path.join(kcorrect.KCORRECT_DIR, 'data',
                            'templates',
                            'kcorrect-default-v4.fits')
    templates = kcorrect.template.SED(filename=filename)

    responses = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0',
                 'twomass_Ks', 'wise_w4']

    fitter = kcorrect.fitter.Fitter(templates=templates,
                                    responses=responses,
                                    redshift_range=[0., 0.1],
                                    nredshift=100)
    fitter.set_Amatrix()

    redshift = 0.0532
    maggies = np.array([1., 1., 1., 1., 1., 1., 1.], dtype=np.float32)
    ivar = np.array([1., 1., 1., 1., 1., 1., 1.], dtype=np.float32)

    coeffs = fitter.fit_coeffs(maggies=maggies, ivar=ivar, redshift=redshift)

    fitspec = coeffs.dot(templates.flux)
    fitmaggies = fitter.reconstruct(redshift=redshift, coeffs=coeffs)

    templates_nointerp = kcorrect.template.SED(filename=filename, interpolate=False)
    fitter_nointerp = kcorrect.fitter.Fitter(templates=templates_nointerp,
                                             responses=responses,
                                             redshift_range=[0., 0.1],
                                             nredshift=100)
    fitter_nointerp.set_Amatrix()

    coeffs_nointerp = fitter_nointerp.fit_coeffs(maggies=maggies, ivar=ivar, redshift=redshift)

    fitspec_nointerp = coeffs_nointerp.dot(templates_nointerp.flux)
    fitmaggies_nointerp = fitter_nointerp.reconstruct(redshift=redshift, coeffs=coeffs_nointerp)

    inz = np.where(fitspec > 0)[0]
    assert np.max(np.abs((fitspec[inz] - fitspec_nointerp[inz]) / fitspec[inz])) < 1.e-6
    assert np.max(np.abs(fitmaggies - fitmaggies_nointerp) < 1.e-6)
    return
