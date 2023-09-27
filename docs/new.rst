
.. _new:


What's New in kcorrect 5
=========================

For this version the software has been migrated to pure Python. There
is no longer a C or IDL interface. There are fewer high-level
utilities as well; basically, more of the data preparation is left to
the users. The new version does not require any environmental
variables to be set.

The methodology is otherwise almost the same, using the same set of
templates as ``v4``.

As of 5.1.0, the responses have been converted to a
``fixed_width`` format readable by ``astropy.io.ascii``, and the 
dependencies on ``fitsio`` and ``pydl`` have been removed.

The only methodological difference to the K-corrections is that the
nonnegative least squares is performed with :py:func:`the SciPy NNLS function
<scipy.optimize.nnls>`
rather than the iterative method used in the past. Both methods solve
the same minimization problem so the answers are very close to one
another.

The stellar masses and absolute magnitudes are no longer referenced to
:math:`h=1`. Instead, they are given by default using the ``Planck18``
version of the :py:class:`FlatLambdaCDM
<astropy.cosmology.FlatLambdaCDM>` cosmology, with the Hubble Constant
for that cosmology.

If you are still using ``v4``, please look at its `online
documentation <http://kcorrect.org>`_ for more information.
