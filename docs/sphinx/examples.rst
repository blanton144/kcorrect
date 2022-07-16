.. _examples:


Examples
=========================

Basic K-correction Calculation
------------------------------

.. code::

   import kcorrect.kcorrect

   responses = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0']
   kc = kcorrect.kcorrect.Kcorrect(responses=responses)

   redshift = 0.0826
   maggies = [27.26e-9, 73.98e-9, 132.56e-9, 198.52e-9, 238.05e-9]
   ivar = [1.528e+18, 5.23e+18, 2.429e+18, 1.042e+18, 0.1653e+18]
   
   coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)
   kc.kcorrect(redshift=redshift, coeffs=coeffs)
