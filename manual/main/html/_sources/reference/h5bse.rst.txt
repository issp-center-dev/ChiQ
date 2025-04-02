.. _Reference_h5bse:

Access to HDF5 file
-------------------

``h5BSE`` class
~~~~~~~~~~~~~~~~

The data in HDF5 files should be accessed via class ``bse.h5bse.h5BSE``.

class methods
^^^^^^^^^^^^^

**__init__(filename, groupname="", compression="gzip", compression_opts=None)**


**open(mode='r')**
    Keep the file open.

**close()**
    Close the file if ``open`` method is called.

**save(key, data)**
    Save data to the file. See the available `key` and data structure below.

**get(key)**  `return data(Dict)`
    Get data from the file.

**get_keylist_data(input_output='input')**  `return List[(key, path), ...]`
    Get the list of keys in the file.


``key`` argument
^^^^^^^^^^^^^^^^^

``key`` argument is a tuple that contains three elements or less. The maximum case is as follows:

.. code-block:: python

    key=(type, w, q)

- **type** `(string)`: ``'beta'``, ``'X_loc'``, ``'chi_q'``, etc. See the table below.

- **w** `(int)`: Index for the bosonic frequency.

- **q** `(string)`: A string to distinguish different q-points.

If ``w`` and/or ``q`` are not necessary, one can omit them.
See :ref:`Reference_hdf5` for the list of available ``type``.

Example
~~~~~~~

This example shows how to access the data (``gamma0`` in this example) in the HDF5 file.

.. code-block:: python

    import numpy as np
    from bse.h5bse import h5BSE


    def _convert_to_matrix(block_matrix, n_block, n_inner):
        assert isinstance(block_matrix, dict)
        dim = n_inner * n_block
        mat = np.zeros((n_block, n_inner, n_block, n_inner), dtype=complex)
        for block, bmat in block_matrix.items():
            i, j = block
            mat[i, :, j, :] = bmat
        return mat.reshape((dim, dim))


    def get_gamma0(h5_file):
        BS = h5BSE(h5_file)
        BS.open('r')

        # get info
        block_name = BS.get(key=('block_name'))
        inner_name = BS.get(key=('inner_name'))
        print(f"block_name: {block_name}")
        print(f"inner_name: {inner_name}")
        n_block = len(block_name)
        n_inner = len(inner_name)
        dim = n_inner * n_block

        # get gamma0
        gamma0_dict = BS.get(key=('gamma0',))

        BS.close()

        # convert dict to ndarray
        gamma0_mat = _convert_to_matrix(gamma0_dict, n_block, n_inner)
        assert gamma0_mat.shape == (dim, dim)

        return gamma0_mat


    gamma0 = get_gamma0("dmft_bse.h5")
    print(f"\ngamma0 loaded. shape={gamma0.shape}")


Data is obtained using the ``get`` method of the ``bse.h5bse.h5BSE`` class.
The returned data is a dictionary type (block matrix), and should be converted to ndarray.

Standard output of the above script for a six-orbital system with ``j=5/2``:

.. code-block:: bash

    $ python3 get_gamma0.py
    block_name: [b'0-ud-0-ud']
    inner_name: [b'0-0', b'0-1', b'0-2', b'0-3', b'0-4', b'0-5', b'1-0', b'1-1', b'1-2', b'1-3', b'1-4', b'1-5', b'2-0', b'2-1', b'2-2', b'2-3', b'2-4', b'2-5', b'3-0', b'3-1', b'3-2', b'3-3', b'3-4', b'3-5', b'4-0', b'4-1', b'4-2', b'4-3', b'4-4', b'4-5', b'5-0', b'5-1', b'5-2', b'5-3', b'5-4', b'5-5']

    gamma0 loaded
      dtype = complex128
      shape = (36, 36)
