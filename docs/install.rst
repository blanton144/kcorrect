
.. _install:


Installation
=========================

The ``kcorrect`` software versions ``5.0.0`` and later
can be installed with ``pip``:

.. code-block:: python

   pip install kcorrect

Alternatively, the software can be found `on GitHub
<https://github.com/blanton144/kcorrect>`_, and cloned and installed
by hand from there.

The SciPy version can affect the results and performance of the
software. Changes to the ``scipy.optimize.nnls`` function lead 
to errors and excessive memory use in SciPy versions 1.12.x, 1.13.x,
and 1.14.x. Versions 1.11.x, and 1.15.x and beyond (so far) are
free of these issues, and we recommend using those versions with
``kcorrect``.
