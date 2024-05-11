.. role:: header_no_toc
  :class: class_header_no_toc

.. title:: kcorrect

:tocdepth: 3

.. rubric:: :header_no_toc:`kcorrect documentation`

<<<<<<< before updating
Overview
=========================
=======
.. code-block:: console

   >> conda create env -n <env_name> python=3.10
   >> conda activate <env_name>


Once you have created a new environment, you can install this project for local
development using the following commands:

.. code-block:: console

   >> pip install -e .'[dev]'
   >> pre-commit install
   >> conda install pandoc


Notes:

1) The single quotes around ``'[dev]'`` may not be required for your operating system.
2) ``pre-commit install`` will initialize pre-commit for this local repository, so
   that a set of tests will be run prior to completing a local commit. For more
   information, see the Python Project Template documentation on
   `pre-commit <https://lincc-ppt.readthedocs.io/en/latest/practices/precommit.html>`_.
3) Installing ``pandoc`` allows you to verify that automatic rendering of Jupyter notebooks
   into documentation for ReadTheDocs works as expected. For more information, see
   the Python Project Template documentation on
   `Sphinx and Python Notebooks <https://lincc-ppt.readthedocs.io/en/latest/practices/sphinx.html#python-notebooks>`_.
>>>>>>> after updating

This page has documentation for the kcorrect product to calculate
astronomical K-corrections for galaxies.

.. toctree::
   :maxdepth: 2

   intro
   new
   conditions
   install
   basics
   examples
   API <autoapi/index>
