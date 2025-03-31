chiq_post
=========

.. code-block:: bash

    chiq_post.py [-h] toml

Description
-----------

``chiq_post.py`` is a post script for outputting text-based data from results stored in HDF5 file. For the basic idea for how to analyze the susceptibility, see :ref:`Algorithm_Eigen`

Positional Arguments
---------------------

**toml**
    Parameter file in TOML format. See :ref:`input` for details.

Options
-------

**-h, --help**
    Show this help message and exit.

Example
-------

Basic usage:

.. code-block:: console

    $ chiq_post.py parameters.toml
