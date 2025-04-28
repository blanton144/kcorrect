.. _templates:


Templates
=========================

``kcorrect`` requires SED templates, which it fits to the given
photometry. There is a default set of templates which is limited and
works well at low redshifts. However, users can provide their own
templates. We show an example using a very broad set of templates that
is useful for observations that extend to the infrared.

If you would like to create custom templates, and have questions about
the data format for the template files, please contact me!

Default Templates
-----------------

The default set of templates are the five templates derived in the
`kcorrect paper
<https://ui.adsabs.harvard.edu/abs/2007AJ....133..734B/abstract>`_.
These were developed for ``kcorrect`` version ``v4``.

They were designed to explain a variety of low redshift photometric
data sets and some spectroscopic indices from the SDSS, and are only
appropriate in the rest-frame ultraviolet through the near-infrared.

If you would like to examine these templates, the following code will
read them into a :py:class:`Template <kcorrect.template.Template>`
object. This contains ``restframe_wave`` (in Angstroms) and
``restframe_flux`` attributes. The flux is in erg/cm^2/s/A at 10 pc.
The corresponding values for the integrated star-formation history,
the remaining stellar mass, and other quantities, are available in the
object (and are read from the FITS file referenced below).

.. code::

   import kcorrect
   import kcorrect.template

   filename = os.path.join(kcorrect.KCORRECT_DIR, 'data',
                           'templates',
                           'kcorrect-default-v4.fits')
   templates = kcorrect.template.Template(filename=filename)


Broad Set of Templates
----------------------

For the purposes of examining mid-infrared active galactic nuclei, in
`Pai and Blanton (2024)
<https://ui.adsabs.harvard.edu/abs/2024ApJ...977..102P/abstract>`_ we
created a very large template set, including 57,600 simple stellar
populations and 36 mid-IR active galactic nuclei models from CLUMPY,
described in `Nenkova et al (2008)
<https://ui.adsabs.harvard.edu/abs/2008ApJ...685..160N/abstract>`_.

Although this sounds a bit ridiculously large, for reasonably sized
samples it is tractable these days to fit non-negatively to large
numbers of templates. There are two advantages of doing so. First, we
do not have to reduce the dimensionality of the problem in any clever
way; if there is a reasonable stellar population and active galactic
nucleus fit to some set of photometry, the optimization will find
it. Second, by using a Monte Carlo technique to estimate the resulting
errors in the fits, we can obtain errors that more accurately reflect
our lack of knowledge of the underlying variations in the
templates. The disadvantage is of course that, though tractable, it is
still pretty slow to do the fits.

The templates can be found at `this Dropbox link (1.3 GB)
<https://www.dropbox.com/scl/fi/0xvt955y1h55vqa426456/templates_broad.fits?rlkey=duwuxyg8j223r6a6pc2im8s6t&st=0xzoq9zg&dl=0>`_. If
you have questions about the specific format of this file, please
contact me. The code below demonstrates how to use this file. Note
that these examples take some time to complete (of order an hour).

The first set of code pre-processes the templates in the context of a
particular redshift range of interest and a particular set of response
curves. This is an expensive operation, but it only need be performed
once. The template fitting calculation will involve an interpolation
of the resulting table over redshift, so you do not want the number of
redshifts (``nredshift``) to be too small. The :py:class:`Kcorrect
<kcorrect.kcorrect.Kcorrect>` object has a method :py:func:`tofits
<kcorrect.kcorrect.Kcorrect.tofits>` that allows you to save the
result to disk for later use.

.. code::

   import kcorrect.template
   import kcorrect.kcorrect

   responses = ['galex_FUV', 'galex_NUV',
                'decam_g', 'decam_r', 'decam_z',
                'wise_w1', 'wise_w2', 'wise_w3', 'wise_w4']

   templates = kcorrect.template.Template(filename='templates_broad.fits')

   kc = kcorrect.kcorrect.Kcorrect(responses=responses,
                                   templates=templates,
                                   redshift_range=[-0.002, 0.4],
                                   nredshift=50)

   kc.tofits('kcorrect_broad.fits')

This preprocessing step can take a lot of memory (up to about 20 GB),
because the both the templates and an interpolation object for all the
templates are stored. We can tell the ``Template`` class to not store
this object by specifying ``interpolate=False``.  In this case, the
interpolation necessary is done when the :py:class:`Kcorrect
<kcorrect.kcorrect.Kcorrect>` object is instantiated, and a full
interpolation object is never created. This is a bit slower
performance but uses much less memory (less than 4 GB).

.. code::

   import kcorrect.template
   import kcorrect.kcorrect

   responses = ['galex_FUV', 'galex_NUV',
                'decam_g', 'decam_r', 'decam_z',
                'wise_w1', 'wise_w2', 'wise_w3', 'wise_w4']

   templates = kcorrect.template.Template(filename='templates_broad.fits',
                                          interpolate=False)

   kc = kcorrect.kcorrect.Kcorrect(responses=responses,
                                   templates=templates,
                                   redshift_range=[-0.002, 0.4],
                                   nredshift=50)

   kc.tofits('kcorrect_broad.fits')

Once the ``kcorrect_broad.fits`` file exists, it can be used to
actually perform template fitting and K-correction determination. Note
that whenever you are using a new set of responses or templates, you
should always make sure that the SED-fitting is adequately
reconstructing the original set of maggies you gave it; as shown
below, you can do that with the :py:func:`reconstruct
<kcorrect.kcorrect.Kcorrect.reconstruct>` method.

.. code::

   import kcorrect.kcorrect

   redshift = 0.030317

   maggies = [5.8345693e-09, 2.3105990e-08,
              1.9847762e-06, 4.2561787e-06, 5.6498902e-06,
              5.0805811e-06, 2.7973852e-06, 1.9943982e-06, 1.2780572e-06]

   ivar = [1.11281764e+18, 1.61874720e+18,
           2.82055401e+14, 6.13363357e+13, 3.48078746e+13,
           4.30457938e+13, 1.41988431e+14, 1.22633740e+14, 4.96459265e+12]

   kc = kcorrect.kcorrect.Kcorrect(filename='kcorrect_broad.fits')

   # For this case, coeffs is large! [1, 57636]
   # If you look carefully at this case, only 6 of the coefficients are non-zero!
   coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=ivar)

   # Check the reconstructed maggies against the original
   rmaggies = kc.reconstruct(redshift=redshift, coeffs=coeffs)

   # We can then calculate the absolute magnitudes as usual
   absmag = kc.absmag(redshift=redshift, maggies=maggies, ivar=ivar, coeffs=coeffs)

