
.. _intro:


Introduction
=========================

This Python module is for calculating K-corrections in the
ultraviolet, optical, and near-infrared.  The method, as well as
results from an earlier version of the software, are described in a
`paper in the Astronomical Journal
<https://ui.adsabs.harvard.edu/abs/2007AJ....133..734B/abstract>`_.

The primary class is :ref:`Kcorrect<kcorrect>`, which fits very restricted spectral
energy distribution models to galaxy photometry. These fits are based
on stellar population synthesis models and may also be interpreted
physically (e.g. in terms of stellar masses or mass-to-light ratios),
with a somewhat limited accuracy.

Of potential additional use are:

* :ref:`Fitter<fitter>`: a general photometry fitter class for users who want to fit their own SEDs.

* :ref:`Template<template>`: a class for expressing templates.
		
* :ref:`Response<response>`: a class for expressing bandpass responses, and projecting SEDs onto them.

The easiest interface to the code, if it is available to you, is the IDL code. Detailed documentation is available for all of the routines. For users of the IDL version of the code, we have special high-level routines to handle data from a number of commonly used surveys: SDSS, GALEX, 2MASS, DEEP2 and GOODS. However, the lower level IDL routines are written generally enough that one can use those routine to handle data from any survey in the restframe UV, optical, or near-infrared.

The lowest level code compiles into a shared object library callable by C, so that non-IDL users can incorporate the K-correction routines directly into their code. It is possible in principle to link the C libraries into code based on SM, TCL/Tk, or Python, and I would be interested in helping interested parties to do this.

One can imagine using the results of this code to calculate the evolution of the luminosity function, the distribution of galaxy colors, as well as to develop galaxy classification algorithms.
