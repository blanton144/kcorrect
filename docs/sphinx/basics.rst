.. _basics:


K-correction Basics
=========================

The basics of K-corrections are well-described in `Hogg et al. (2002)
<https://ui.adsabs.harvard.edu/abs/2002astro.ph.10394H/abstract>`_.
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

The band :math:`Q` can in principle be anything; it does not have to
be the same as :math:`R`, and it also does not have to be an actual
bandpass at all. A particular choice is a "shifted" bandpass, where a
band shift of :math:`z` would be denoted :math:`^{z}Q` and would
indicate that the rest-frame band pass :math:`Q` blue shifted by a
factor :math:`(1+z)`. In this case the K-correction for an object at
redshift :math:`z` from observed bandpass :math:`R` to a rest frame
band pass :math:`^{z}R` would be independent of the object's SED and
equal to :math:`-2.5\log_{10}(1+z)`, because the bandpasses exactly
overlap. The advantage of this choice is that by choosing the band
shift to be near the typical redshift of a sample, one can minimize
the errors due to K-corrections when comparing objects within the
sample.

To perform the fits, the software requires broad band flux
measurements.  These are accepted as AB maggies. These are the ratio
of the source to the AB standard source in each band. To relate these
quantities to magnitudes, an object with total flux f in maggies has
magnitude

m = − 2.5 log 10 f .
An advantage of the maggie unit system is that it is linear, and thus can when necessary accommodate negative flux estimates. Two notes:

As discussed below, SDSS catalog numbers are NOT on our best guess for the AB system. A set of offsets must be applied to the magnitudes to achieve our best guess. Our high level IDL code deals with this automatically.
SDSS catalog magnitudes obtained from the official survey database are luptitudes, which for reasonably bright objects are equivalent to magnitudes. While maggies are simply related to magnitudes by 10−0.4m, the conversion from luptitudes to maggies is a bit more complicated; see the description accompanying the DR2 documentation). Our high level IDL code deals with these conversions automatically.
But you may still be wondering what I mean by an "AB" magnitude. The AB system is designed such that a flat spectrum object with fν = 3631 Jy = 3.631 × 10−20 ergs cm−2 s−1 Hz−1 should have every magnitude equal to zero. The beauty of the AB system is that the uniform definition makes it convenient to synthesize AB magnitudes from theory or models. The tragic flaw is that the quality of the AB-ness of a system is very dependent on precise spectrophotometry of standards and the carefulness of the calibrators, since no objects have a flat spectrum. There is a tension between these two needs --- similar to other tensions throughout astronomy between making precise measurements and making interpretable ones.

Finally, I have also included code to calculate photometric redshifts based on the templates. This procedure is as simple as running the K-correction code at each redshift and finding that redshift which provides the best fit in the χ2 sense.
