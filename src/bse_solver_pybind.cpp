#include <iostream>

#include "./block_matrix_pybind.hpp" // Include block matrix Python bindings
#include "./bse.hpp" // Include BSE solver header

namespace py = pybind11; // Namespace alias for pybind11

/**
 * @brief Python binding for BSE (Bethe-Salpeter Equation) solver
 *
 * This class provides Python bindings for solving the Bethe-Salpeter equation.
 * It wraps the C++ BSE solver implementation and exposes key functionality to Python.
 */
class BSESolver {
public:
  /**
   * @brief Constructor for BSESolver
   *
   * @param beta Inverse temperature
   * @param matrix_info_A Block matrix information for A matrices
   * @param matrix_info_B Block matrix information for B matrices
   * @param matrix_info_C Block matrix information for C matrices
   */
  BSESolver(double beta, const std::vector<int> &matrix_info_A, const std::vector<int> &matrix_info_B, const std::vector<int> &matrix_info_C)
    : matrix_info_A(matrix_info_A), // Store A matrix info
      matrix_info_B(matrix_info_B), // Store B matrix info
      matrix_info_C(matrix_info_C), // Store C matrix info
      solver(beta, matrix_info_A, matrix_info_B, matrix_info_C) {} // Initialize solver

  /**
   * @brief Set input matrices for the solver
   *
   * @param d Python dictionary containing matrix data
   * @param str Matrix type identifier string
   *
   * Supported matrix types:
   * - "X_loc": Local X matrix
   * - "X0_loc": Local X0 matrix
   * - "X0_q": X0(q) matrix
   * - "chi_loc": Local susceptibility matrix
   * - "chi0_loc": Local bare susceptibility matrix
   * - "gamma0": Bare vertex matrix
   * - "Phi": Phi matrix
   * - "Phi_sum": Summed Phi matrix
   */
  void setMatrix(const py::dict &d, const py::str &str) {
    std::string type(str); // Convert Python string to C++ string

    // Set X_loc matrix
    if (type == "X_loc") {
      X_loc = dict_to_blockmatrix<std::complex<double> >(d, matrix_info_A);
      solver.SetXLoc(X_loc);
    }
    // Set X0_loc matrix
    else if (type == "X0_loc") {
      X0_loc = dict_to_blockmatrix<std::complex<double> >(d, matrix_info_B);
      solver.SetX0Loc(X0_loc);
    }
    // Set X0(q) matrix
    else if (type == "X0_q") {
      X0_q = dict_to_blockmatrix<std::complex<double> >(d, matrix_info_B);
      solver.SetX0q(X0_q);
    }
    // Set chi_loc matrix
    else if (type == "chi_loc") {
      chi_loc = dict_to_blockmatrix<std::complex<double> >(d, matrix_info_C);
      solver.SetChiLoc(chi_loc);
    }
    // Set chi0_loc matrix
    else if (type == "chi0_loc") {
      chi0_loc = dict_to_blockmatrix<std::complex<double> >(d, matrix_info_C);
      solver.SetChi0Loc(chi0_loc);
    }
    // Set chi0_q matrix
    else if (type == "chi0_q") {
      chi0_q = dict_to_blockmatrix<std::complex<double> >(d, matrix_info_C);
      solver.SetChi0q(chi0_q);
    }
    // Set gamma0 matrix
    else if (type == "gamma0") {
      gamma0 = dict_to_blockmatrix<std::complex<double> >(d, matrix_info_C);
      solver.SetGamma0(gamma0);
    }
    // Set Phi matrix
    else if (type == "Phi") {
      Phi = dict_to_blockmatrix<std::complex<double> >(d, matrix_info_B);
      solver.SetPhi(Phi);
    }
    // Set Phi_sum matrix
    else if (type == "Phi_sum") {
      Phi_sum = dict_to_blockmatrix<std::complex<double> >(d, matrix_info_C);
      solver.SetPhi_Sum(Phi_sum);
    }
    // Handle invalid matrix type
    else{
      std::stringstream ss;
      ss << "Invalid type '" << type << "'";
      throw std::runtime_error(ss.str());
    }
  }

  /**
   * @brief Calculate various quantities based on the specified type
   *
   * @param str Calculation type identifier string
   *
   * Supported calculation types:
   * - "chi0": Calculate bare susceptibility
   * - "BSE"/"bse": Solve full BSE
   * - "RPA"/"rpa": Calculate RPA susceptibility
   * - "RRPA"/"rrpa": Calculate renormalized RPA
   * - "SCL"/"scl": Calculate using SCL approximation
   */
  void calc(const py::str &str) {
    std::string type(str); // Convert Python string to C++ string

    // Calculate bare susceptibility
    if (type == "chi0") {
      solver.CalcChi0q();
    }
    // Solve full BSE
    else if (type == "BSE" || type == "bse") {
      solver.CalcChiq_BSE();
    }
    // Calculate RPA susceptibility
    else if (type == "RPA" || type == "rpa") {
      solver.CalcChiq_RPA();
    }
    // Calculate renormalized RPA
    else if (type == "RRPA" || type == "rrpa") {
      solver.CalcChiq_RRPA();
    }
    // Calculate using SCL approximation
    else if (type == "SCL" || type == "scl") {
      solver.CalcChiq_SCL();
    }
    // Handle invalid calculation type
    else{
      std::stringstream ss;
      ss << "Invalid type '" << type << "'";
      throw std::runtime_error(ss.str());
    }
  }

  /**
   * @brief Get calculated matrices
   *
   * @param str Matrix type identifier string
   * @return Python dictionary containing the requested matrix
   *
   * Supported matrix types:
   * - "chi0_q": Bare susceptibility
   * - "chi_q": Full BSE susceptibility
   * - "chi_q_rpa": RPA susceptibility
   * - "chi_q_rrpa": Renormalized RPA susceptibility
   * - "chi_q_scl": SCL susceptibility
   * - "I_q": I(q) matrix
   * - "I_q_scl": I(q) matrix in SCL
   */
  py::dict getMatrix(const py::str &str) {
    std::string type(str); // Convert Python string to C++ string

    bse::block_matrix<std::complex<double> > bm;
    // Get bare susceptibility
    if (type == "chi0_q") {
      bm = solver.GetChi0q();
    }
    // Get full BSE susceptibility
    else if (type == "chi_q") {
      bm = solver.GetChiq_BSE();
    }
    // Get RPA susceptibility
    else if (type == "chi_q_rpa") {
      bm = solver.GetChiq_RPA();
    }
    // Get renormalized RPA susceptibility
    else if (type == "chi_q_rrpa") {
      bm = solver.GetChiq_RRPA();
    }
    // Get SCL susceptibility
    else if (type == "chi_q_scl") {
      bm = solver.GetChiq_SCL();
    }
    // Get I(q) matrix
    else if (type == "I_q") {
      bm = solver.GetIq();
    }
    // Get I(q) matrix in SCL
    else if (type == "I_q_scl") {
      bm = solver.GetIq();
    }
    // Handle invalid matrix type
    else{
      std::stringstream ss;
      ss << "Invalid type '" << type << "'";
      throw std::runtime_error(ss.str());
    }

    return blockmatrix_to_dict<std::complex<double> >(bm); // Convert block matrix to Python dict
  }

private:
  bsesol::bse_solver<std::complex<double> > solver; ///< Underlying BSE solver instance
  std::vector<int> matrix_info_A; ///< Block matrix information for A matrices
  std::vector<int> matrix_info_B; ///< Block matrix information for B matrices
  std::vector<int> matrix_info_C; ///< Block matrix information for C matrices

  // Input matrices
  bse::block_matrix<std::complex<double> > X_loc; ///< Local X matrix
  bse::block_matrix<std::complex<double> > X0_loc; ///< Local X0 matrix
  bse::block_matrix<std::complex<double> > X0_q; ///< X0(q) matrix
  bse::block_matrix<std::complex<double> > chi_loc; ///< Local susceptibility
  bse::block_matrix<std::complex<double> > chi0_loc; ///< Local bare susceptibility
  bse::block_matrix<std::complex<double> > chi0_q; ///< Lattice bare susceptibility chi0(q)
  bse::block_matrix<std::complex<double> > gamma0; ///< Bare vertex
  bse::block_matrix<std::complex<double> > Phi; ///< Phi matrix
  bse::block_matrix<std::complex<double> > Phi_sum; ///< Summed Phi matrix
};

/**
 * @brief Pybind11 module definition for BSE solver
 */
PYBIND11_MODULE(bse_solver, m) {
  // Create Python class binding
  py::class_<BSESolver>(m, "BSESolver")
	.def(py::init<double, const std::vector<int> &, const std::vector<int> &, const std::vector<int> &>()) // Constructor binding
	.def("set", &BSESolver::setMatrix) // Bind setMatrix method
	.def("calc", &BSESolver::calc) // Bind calc method
	.def("get", &BSESolver::getMatrix) // Bind getMatrix method
  ;
}
