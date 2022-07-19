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
