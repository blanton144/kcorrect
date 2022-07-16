
.. _new:


What's New in kcorrect 5
=========================

For this version the software has been migrated to pure Python. There
is no longer a C or IDL interface. There are fewer high-level
utilities as well; basically, more of the data preparation is left to
the users.

The methodology is otherwise almost the same. The only methodological
difference is that the nonnegative least squares is performed with
`SciPi's
NNLS<https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html>`_
rather than the iterative method used in the past. Both methods solve
the same minimization problem so the answers are very close to one
another.
