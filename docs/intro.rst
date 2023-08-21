
.. _intro:

:tocdepth: 3


Introduction
=========================

This Python module is for calculating K-corrections of galaxies in the
ultraviolet, optical, and near-infrared.  The method, as well as
results from an earlier version of the software, are described in a
`paper in the Astronomical Journal
<https://ui.adsabs.harvard.edu/abs/2007AJ....133..734B/abstract>`_.

The primary class is :py:class:`Kcorrect <kcorrect.kcorrect.Kcorrect>`, which fits very restricted spectral
energy distribution models to galaxy photometry. These fits are based
on stellar population synthesis models and may also be interpreted
physically (e.g. in terms of stellar masses or mass-to-light ratios),
with a somewhat limited accuracy.

Of potential additional use are:

* :py:class:`Fitter <kcorrect.fitter.Fitter>`: a general photometry fitter class for users who want to fit their own SEDs.

* :py:class:`Template <kcorrect.template.Template>`: a class for expressing templates.
		
* :py:class:`Response <kcorrect.response.Response>`: a class for expressing bandpass responses, and projecting SEDs onto them.

