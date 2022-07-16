.. _basics:


K-correction Basics
=========================

The basics of K-corrections are well-described in `Hogg et al. (2002)
<https://ui.adsabs.harvard.edu/abs/2002astro.ph.10394H/abstract>`_ and
in `Blanton & Roweis (2007)
<https://ui.adsabs.harvard.edu/abs/2007AJ....133..734B/abstract>`_.
Here is a briefer description including how they are estimated by
``kcorrect``.

If you want to convert apparent magnitudes in band R to absolute
magnitudes in band Q, you need to calculate the K-correction, which is
defined by the equation:

.. math::

  m_R = M_Q + {\rm DM}(z) + K_{QR}(z),

where :math:`m_R` is the apparent magnitude, :math:`M_Q` is the
absolute magnitude, :math:`{\rm DM}(z)` is the distance modulus,
accounting for the angular diameter distance and cosmological
surface-brightness dimming, and :math:`K_{QR}(z)` is the
K-correction.

By absolute magnitude we mean: the apparent magnitude in band
:math:`Q` that the object would have if it were observed at rest, 10
pc away, using an aperture that contains its total flux. The distance
modulus accounts for the difference between an object's actual
distance and 10 pc. The K-correction accounts for the fact that you
observed a redshifted galaxy in band :math:`R` but the absolute
magnitude requires a rest-frame observation in band
:math:`Q`. Obviously the difference between the fluxes observed in
different bandpasses is fully determined by the galaxy SED and the
description of the bandpasses.

In order to get the appropriate SED for a set of galaxy fluxes,
``kcorrect`` fits an SED which is a nonnegative linear combination of
some small number of templates. The templates have been optimized to
minimize the residuals between the actual galaxy fluxes and the galaxy
fluxes reconstructed from the galaxy SED fit. The K-correction is then
calculated from this best-fit SED.

To perform the fits, the software requires broad band flux
measurements.  These are accepted as AB maggies. Maggies are the ratio
of the source to the AB standard source in each band, using the
integrals in the `kcorrect paper
<https://ui.adsabs.harvard.edu/abs/2007AJ....133..734B/abstract>`_.
They have a simple relationship to magnitudes:

.. math::
   m = − 2.5 \log 10 \mu.

where :math:`m` is the magnitude and :math:`\mu` is maggies. An
advantage of the maggie unit system relative to magnitudes is that it
is linear, and thus can when necessary accommodate negative flux
estimates.

The AB standard source is a flat spectrum object with :math:`f_\nu =
3631 {\rm ~Jy} = 3.631 \times 10^{−20} {\rm ~ergs} {\rm ~cm}^{−2} {\rm
~s}^{−1} {\rm ~Hz}^{−1}`. Such a source would have all magnitudes
equal to zero.

There are still many available catalogs that are defined on the Vega
standard source system. The :py:class:`Response
<kcorrect.response.Response>` object associated with a filter has the
``vega2ab`` attribute which defines the correction from Vega to AB in
magnitudes. 

The band :math:`Q` can in principle be anything; it does not have to
be the same as :math:`R`, and it also does not have to be an actual
bandpass at all.

``kcorrect`` supports a particular choice of a "shifted" bandpass with
its ``band_shift`` option, where a band shift of :math:`z` would be
denoted :math:`^{z}Q` and would indicate that the rest-frame band pass
:math:`Q` blue shifted by a factor :math:`(1+z)`. In this case the
K-correction for an object at redshift :math:`z` from observed
bandpass :math:`R` to a rest frame band pass :math:`^{z}R` would be
independent of the object's SED and equal to
:math:`-2.5\log_{10}(1+z)`, because the bandpasses exactly
overlap. The advantage of this choice is that by choosing the band
shift to be near the typical redshift of a sample, one can minimize
the errors due to K-corrections when comparing objects within the
sample.

Two final notes on SDSS units, which since ``kcorrect`` grew out of SDSS
work seems appropriate here.

* The official SDSS catalog numbers are not our best guess for the AB
  system, so SDSS data has a small conversion factor that needs to be
  applied using :py:func:`sdss_ab_correct()
  <kcorrect.utils.sdss_ab_correct>`.

* The "magnitudes" published by SDSS are so-called "asinh magnitudes"
  or "luptitudes", described by `the SDSS imaging documentation
  <https://www.sdss.org/dr17/algorithms/magnitudes/#asinh>`_. These
  can be converted to maggies with :py:func:`sdss_asinh_to_maggies()
  <kcorrect.utils.sdss_asinh_to_maggies>`.
