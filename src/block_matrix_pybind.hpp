// g++ -O2 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` wrapper.cpp -o wrapper`python3-config --extension-suffix` -I /usr/include/eigen3

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <iostream>
// #include <typeinfo>

#include "./block_matrix.hpp"

/**
 * Conversion between bse::block_matrix class and
 * python dictionary (pybind11::dict).
 */

/**
 * Convert python dictionary to bse::block_matrix
 *
 * Parameters
 * ----------
 * d : pybind11::dict
 *     Input Python dictionary containing block matrices
 * matrixInfo : std::vector<int>
 *     Vector containing information about matrix block sizes
 *
 * Returns
 * -------
 * bse::block_matrix<Scalar>
 *     Block matrix constructed from input dictionary
 *
 * Raises
 * ------
 * std::runtime_error
 *     If inconsistent data type is found in dictionary blocks
 */
template<typename Scalar>
bse::block_matrix<Scalar> dict_to_blockmatrix(const pybind11::dict &d, const std::vector<int> &matrixInfo){
    bse::block_matrix<Scalar> bm(matrixInfo);

    using matrix_t = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    for (auto item : d){
        std::tuple<int, int> key = item.first.cast<std::tuple<int, int> >();
        int i = std::get<0>(key);
        int j = std::get<1>(key);
        // std::cout << " (i, j) = (" << i << ", " << j << ")" << std::endl;

        // Check datatype (double, complex<double>, etc)
        if(!pybind11::isinstance<pybind11::array_t<Scalar> >(item.second)){
            std::stringstream ss;
            ss << "Inconsistent datatype in (" << i << ", " << j << ") block.";
            throw std::runtime_error(ss.str());
        }

        matrix_t m = item.second.cast<matrix_t>();
        bm.assign(i, j, m);
    }
    return bm;
}


/**
 * Convert bse::block_matrix to python dictionary
 *
 * Parameters
 * ----------
 * bm : bse::block_matrix<Scalar>
 *     Input block matrix to convert
 *
 * Returns
 * -------
 * pybind11::dict
 *     Python dictionary containing the block matrices
 */
template<typename Scalar>
pybind11::dict blockmatrix_to_dict(const bse::block_matrix<Scalar> &bm) {
    pybind11::dict d;
    for (int i = 0; i < bm.num_blocks(); ++i) {
      for (int j = 0; j < bm.num_blocks(); ++j) {
        if (bm.exists(i, j)) {
          pybind11::tuple key = pybind11::make_tuple(i, j);
          d[key] = bm(i, j);
        }
      }
    }
    return d;
}
