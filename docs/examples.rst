.. _examples:


Examples
=========================

Basic K-correction Calculation
------------------------------

The way ``kcorrect`` works is that the user initializes an object of
the :py:class:`Kcorrect class <kcorrect.kcorrect.Kcorrect>` for the
specific set of response functions for the photometry they are
interested in. This object takes a bit of time to calculate, because
it precalculates the projections of the K-correction templates onto
the response function on a grid of redshifts.

Given redshifts and photometry in units of maggies (and their inverse
variance), the object's :py:func:`fit_coeffs <kcorrect.kcorrect.Kcorrect.fit_coeffs>` method fits the templates to
the photometry. 

These coefficients can then be used to calculate the kcorrections with
the :py:func:`kcorrect <kcorrect.kcorrect.Kcorrect.kcorrect>` method.

Below is an example run on a single object, for which ``redshift`` is
a scalar, and ``maggies`` and ``ivar`` are 1-dimensional arrays. We
can input the arrays as ndarrays or as lists.

.. code::

   import kcorrect.kcorrect

   responses = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0']
   kc = kcorrect.kcorrect.Kcorrect(responses=responses)

   redshift = 0.0826
   maggies = [27.26e-9, 73.98e-9, 132.56e-9, 198.52e-9, 238.05e-9]
   ivar = [1.528e+18, 5.23e+18, 2.429e+18, 1.042e+18, 0.1653e+18]
   
   # "coeffs" is a [5]-array with coefficients multiplying each template
   coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

   # "k" is a [5]-array with the K-corrections in magnitude units
   k = kc.kcorrect(redshift=redshift, coeffs=coeffs)


For multiple objects (and optionally for single objects), ``redshift``
is a 1-dimensional array, whereas ``maggies`` and ``ivar`` are
2-dimensional arrays. Again, these can be input as lists too, as shown
below.

.. code::

   import kcorrect.kcorrect

   responses = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0']
   kc = kcorrect.kcorrect.Kcorrect(responses=responses)

   redshift = [0.0826, 0.111]
   maggies = [[27.26e-9, 73.98e-9, 132.56e-9, 198.52e-9],
	            [4.022e-9, 29.36e-9, 98.230e-9, 155.63e-9]]
   ivar = [[1.528e+18, 5.23e+18, 2.429e+18, 1.042e+18],
	         [2.360e+18, 9.87e+18, 3.006e+18, 1.122e+18]]  
   
   
   # "coeffs" is a [2,5]-array coefficients multiplying each template
   coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

   # "k" is a [2,4]-array with the K-corrections in magnitude units
   k = kc.kcorrect(redshift=redshift, coeffs=coeffs)


Absolute magnitudes
-------------------

The :py:func:`absmag <kcorrect.kcorrect.Kcorrect.absmag>` method
returns absolute magnitudes.

The code returns:

.. math::

  m_Q = m_R - {\rm DM}(z) - K_{QR}(z),

for those bands for which the observed maggies are positive and have a
positive inverse variance.  For other bands, the returned
absolute magnitude is determined using the reconstructed maggies (see
next section).

.. code::

   absmag = kc.absmag(redshift=redshift, maggies=maggies, ivar=ivar, coeffs=coeffs)

If one of the maggies is zero or negative, the corresponding absolute magnitude
is returned as ``-9999``. In the code below, the u-band absolute magnitude
of the first object should have this value.

.. code::

    maggies = [[-3.26e-9, 73.98e-9, 132.56e-9, 198.52e-9],
               [4.022e-9, 29.36e-9, 98.230e-9, 155.63e-9]]
    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)
    absmag = kc.absmag(redshift=redshift, maggies=maggies, ivar=ivar, coeffs=coeffs)

In these cases, the method can calculate the one-sigma absolute magnitude limit based
on the inverse variance, if it is not zero, by setting the ``limit`` keyword to ``True``.
This returns the one-sigma detection limit for all objects and bands (even the detected
ones). 
 
.. code::

    absmag, absmag_limit = kc.absmag(redshift=redshift, maggies=maggies, ivar=ivar, coeffs=coeffs, limit=True)

If one of the input inverse variances is zero or negative, the corresponding
absolute magnitude and (if it is requested) the limit are returned as
``-9999``.

.. code::

    ivar = [[0., 0., 2.429e+18, 1.042e+18],
		        [4.022e-9, 29.36e-9, 98.230e-9, 155.63e-9]]
    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)
    absmag, absmag_limit = kc.absmag(redshift=redshift, maggies=maggies, ivar=ivar, coeffs=coeffs, limit=True)

There is also an option to return the reconstructed absolute magnitude
from the full SED fit. This will differ somewhat from the absolute magnitude
generated through applying the K-correction, of course, because the SED fit
is not a perfect fit to the data. It may be a useful quantity in assessing
the goodness of fit (though the reconstructed maggies are a more direct
comparison to the input data, see below). It also may be useful as a guess of
the absolute magnitude for missing bands.

.. code::

    absmag, absmag_reconstruct = kc.absmag(redshift=redshift, maggies=maggies, ivar=ivar, coeffs=coeffs, reconstruct=True)

Reconstructing spectra and fluxes
---------------------------------

The model expressed by the coefficients is a full SED. The
:py:class:`Kcorrect <kcorrect.kcorrect.Kcorrect>` object has an
attribute ``templates`` that is an instance of the :py:class:`Template
<kcorrect.kcorrect.Template>` class.

Each template is in units of :math:`{\rm ~erg} {\rm ~s}^{-1} {\rm
~cm}^{-2} {\rm ~A}^{-1} {\rm ~}M_\odot^{-1}` as observed for a galaxy
at 10pc distance. The coefficients are in units of "the solar masses
corresponding to a galaxy at 10pc" (see next section).

It is straightforward to use the templates to reconstruct what that
SED looks like. If you have found coefficients as in the first example
above, the following code will return the SED in :math:`{\rm ~erg}
{\rm ~s}^{-1} {\rm ~cm}^{-2} {\rm ~A}^{-1}`. The templates have to be
shifted to the observed frame (conserving bolometric flux) to get the
correct observed spectrum.

.. code::

   import matplotlib.pyplot as plt
   import numpy as np

   wave = kc.templates.restframe_wave * (1. + redshift)
   spec = coeffs.dot(kc.templates.restframe_flux) / (1. + redshift)

   plt.plot(np.log10(wave), np.log10(spec))
   plt.xlabel('$\\log_{10} wavelength$')
   plt.ylabel('$\\log_{10} flux$ (erg s$^{-1}$ cm$^{-2}$ Ang$^{-1}$')
   plt.xlim([3., 4.5])
   plt.ylim([-18., -14.])
   plt.show()

We can also reconstruct the fluxes from the model (of course we can,
because the K-correction determination must do so!). This is useful to
do to compare the best fit SED to the observations. In the case here
you should find agreement within a few percent between ``maggies`` and
``rmaggies``.

.. code::

    rmaggies = kc.reconstruct(redshift=redshift, coeffs=coeffs)


Derived parameters (i.e. stellar mass)
--------------------------------------

The :py:func:`derived <kcorrect.kcorrect.Kcorrect.derived>` method
returns some derived parameters. These parameters are described in the
`kcorrect paper
<https://ui.adsabs.harvard.edu/abs/2007AJ....133..734B/abstract>`_,
but it is important to take them with a grain of salt.

The only one worth taking at all seriously is ``mremain`` and
``mtol``, the surviving stellar mass and mass-to-light ratios in the
stellar population fit. But this parameter is shown to disagree with
other estimates, with a trend of a few tenths of dex across stellar
mass (``kcorrect`` declining). While I don't know if any stellar mass
indicator from broad band photometry is great, the one in ``kcorrect``
is particularly simple (and also doesn't come with any error bar).

Importantly, the mass-to-light ratios are in the output bandpasses
(and you have to specify ``band_shift`` if you want shifted output
bandpasses).

Like the absolute magnitudes, the stellar masses use the ``cosmo``
attribute of the :py:class:`Kcorrect <kcorrect.kcorrect.Kcorrect>`
object, which by default is the ``Planck18`` cosmology from
``astropy``.

.. code::

   derived = kc.derived(redshift=redshift, coeffs=coeffs)

   # This has one entry per object
   stellar_mass = derived['mtol']

   # This has one entry per object per output bandpass
   mtol = derived['mtol']
	 

Changing the output responses
-----------------------------

As one gets to higher redshift, the K-corrections from a given
observed frame bandpass to its rest frame counterpart become more and
more dependent on the SED model being correct.

One approach to dealing with this is to define a set of output rest
frame bandpasses that are shifted versions of your input bandpasses,
where the shift is the typical redshift of your sample. This minimizes
the internal error in your sample, at the expense of calculating
quantities that are less likely to be comparable to catalogs in the
literature. This option can be utilized by just specifying
``band_shift`` for the :py:func:`kcorrect
<kcorrect.kcorrect.Kcorrect.kcorrect>`, :py:func:`absmag
<kcorrect.kcorrect.Kcorrect.absmag>`, or :py:func:`derived
<kcorrect.kcorrect.Kcorrect.derived>` methods when you use them:

.. code::

   absmag = kc.absmag(redshift=redshift, maggies=maggies, ivar=ivar, coeffs=coeffs, band_shift=0.1)

A second approach is to define a set of output rest frame bandpasses
that correspond closely to the effective rest frame wavelength of the
observed bandpass for galaxies at the typical observed redshift in
your sample. For example, at :math:`z\sim 0.7` the observed SDSS
:math:`r`, :math:`i`, and :math:`z` bands are close in wavelength to
the rest frame :math:`U`, :math:`B`, and :math:`V` bands. So if we
observed galaxies in the SDSS bands at those redshifts, we could find
the K-corrections to :math:`UBV`. To do so, we have to instantiate a
:py:class:`Kcorrect <kcorrect.kcorrect.Kcorrect>` object that can do
so using the ``response_map`` and ``response_out`` arguments:

.. code::

   import kcorrect.kcorrect

   responses_in = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0']
   responses_out = ['bessell_U', 'bessell_B', 'bessell_V']
   responses_map = ['sdss_r0', 'sdss_i0', 'sdss_z0']
   kc = kcorrect.kcorrect.Kcorrect(responses=responses_in,
                                   responses_out=responses_out,
				   responses_map=responses_map)

   # These are the ugriz observations (made up!)
   redshift = 0.72
   maggies = [27.26e-9, 73.98e-9, 132.56e-9, 198.52e-9, 238.05e-9]
   ivar = [1.528e+18, 5.23e+18, 2.429e+18, 1.042e+18, 0.1653e+18]
   
   # "coeffs" is a [5]-array with coefficients multiplying each template
   coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

   # "k" is a [3]-array with the K-corrections in magnitude units,
   # from riz to UBV.
   k = kc.kcorrect(redshift=redshift, coeffs=coeffs)

   # "absmag" is also a [3]-array resuling from applying K-corrections
   # and the distance modulus
   absmag = kc.absmag(redshift=redshift, maggies=maggies, ivar=ivar, coeffs=coeffs)

