#pragma once

#include "./block_matrix.hpp"
#include "./memory_info.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * @brief Namespace containing BSE solver implementation
 */
namespace bsesol {

  /**
   * @brief Class implementing Bethe-Salpeter equation solver
   *
   * This class provides functionality to solve the Bethe-Salpeter equation (BSE)
   * and calculate various susceptibilities including:
   * - BSE susceptibility: Full BSE solution including vertex corrections
   * - RPA susceptibility: Random Phase Approximation
   * - RRPA susceptibility: Renormalized RPA
   * - SCL susceptibility: Strong-Coupling Limit approximation
   *
   * The solver works with three types of block matrices:
   * - A-type: Blocks [block], inner structure [in, iw]
   * - B-type: Blocks [iw, block], inner structure [in]
   * - C-type: Blocks [block], inner structure [in]
   *
   * Where:
   * - block: Block index for orbital/spin channels
   * - in: Inner index within each block
   * - iw: Frequency index
   *
   * @tparam Scalar The scalar type used for calculations (e.g. double, complex<double>)
   */
  template<typename Scalar>
  class bse_solver {
  public:
    /**
     * @brief Constructor for BSE solver
     *
     * Initializes solver with matrix dimensions and temperature
     *
     * @param beta Inverse temperature (1/T)
     * @param matrixInfo_A_ Matrix info for A-type matrices (block: [block], inner: [in, iw])
     * @param matrixInfo_B_ Matrix info for B-type matrices (block: [iw, block], inner: [in])
     * @param matrixInfo_C_ Matrix info for C-type matrices (block: [block], inner: [in])
     * @param verbose Enable verbose memory info output
     */
    bse_solver(double beta,
               const std::vector<int> &matrixInfo_A_,
               const std::vector<int> &matrixInfo_B_,
               const std::vector<int> &matrixInfo_C_,
               bool verbose=false)
            : beta(beta), matrixInfo_A(matrixInfo_A_), matrixInfo_B(matrixInfo_B_), matrixInfo_C(matrixInfo_C_), memInfo(verbose) {
      // Extract dimensions from matrix info
      int n_block = matrixInfo_A.size();  // Number of blocks
      assert (matrixInfo_B.size() % n_block == 0);
      int n_iw = int(matrixInfo_B.size() / n_block);  // Number of frequencies
      int n_in = matrixInfo_C[0];  // Inner dimension

      // Validate A-type matrix dimensions
      assert(matrixInfo_A.size() == n_block);
      for (int i=0; i<n_block; i++){
        assert(matrixInfo_A[i] == n_in * n_iw);  // Each block should have size n_in * n_iw
      }

      // Validate B-type matrix dimensions
      assert(matrixInfo_B.size() == n_iw * n_block);
      for (int i=0; i<n_iw * n_block; i++){
        assert(matrixInfo_B[i] == n_in);  // Each block should have size n_in
      }

      // Validate C-type matrix dimensions
      assert(matrixInfo_C.size() == n_block);
      for (int i=0; i<n_block; i++){
        assert(matrixInfo_C[i] == n_in);  // Each block should have size n_in
      }
    }

    /**
     * @brief Set X0_loc matrix - Local bare susceptibility
     * @param other X0_loc matrix to set
     */
    void SetX0Loc(const bse::block_matrix<Scalar> &other) {
      this->X0_Loc = &other;
      memInfo.add(other, "X0_loc");
      // this->chi0_q_ptr = nullptr;  // Reset chi0_q pointer since it depends on X0_loc
    }

    /**
     * @brief Set X_loc matrix - Local full vertex
     * @param other X_loc matrix to set
     */
    void SetXLoc(const bse::block_matrix<Scalar> &other) {
      this->X_loc = &other;
      memInfo.add(other, "X_loc");
    }

    /**
     * @brief Set X0q matrix - Momentum-dependent bare susceptibility
     * @param other X0q matrix to set
     */
    void SetX0q(const bse::block_matrix<Scalar> &other) {
      this->X0_q = &other;
      memInfo.add(other, "X0_q");
      // this->chi0_q_ptr = nullptr;  // Reset chi0_q pointer since it depends on X0_q
    }

    /**
     * @brief Set chi_loc matrix - Local physical susceptibility
     * @param other chi_loc matrix to set
     */
    void SetChiLoc(const bse::block_matrix<Scalar> &other) {
      this->chi_loc = &other;
    }

    /**
     * @brief Set chi0_loc matrix - Local bare physical susceptibility
     * @param other chi0_loc matrix to set
     */
    void SetChi0Loc(const bse::block_matrix<Scalar> &other) {
      this->chi0_loc = &other;
      // this->chi0_q_ptr = nullptr;  // Reset chi0_q pointer since it depends on chi0_loc
    }

    /// Set chi0q
    /// \param other chi0q
    void SetChi0q(const bse::block_matrix<Scalar> &other) {
      this->chi0_q_in = &other;
      memInfo.add(other, "chi0_q");
    }

    /**
     * @brief Set gamma0 matrix - Bare vertex
     * @param other gamma0 matrix to set
     */
    void SetGamma0(const bse::block_matrix<Scalar> &other) {
      this->gamma0 = &other;
    }

    /**
     * @brief Set Phi matrix - Vertex function
     * @param other Phi matrix to set
     */
    void SetPhi(const bse::block_matrix<Scalar> &other) {
      this->Phi = &other;
      memInfo.add(other, "Phi");
    }

    /*
    /// Set Phi_L - Left vertex function
    /// \param other Phi_L
    void SetPhi_L(const bse::block_matrix<Scalar> &other) {
      this->Phi_L = other;
    }

    /// Set Phi_R - Right vertex function
    /// \param other Phi_R
    void SetPhi_R(const bse::block_matrix<Scalar> &other) {
      this->Phi_R = other;
    }
    */

    /**
     * @brief Set Phi_Sum matrix - Summed vertex function
     * @param other Phi_Sum matrix to set
     */
    void SetPhi_Sum(const bse::block_matrix<Scalar> &other) {
      this->Phi_Sum = &other;
      memInfo.add(other, "Phi_sum");
    }

    /// Get chi_q_bse - BSE susceptibility result
    /// @property
    bse::block_matrix<Scalar>& GetChiq_BSE() {
      return this->chi_q_bse;
    }

    /// Get chi0_q - Bare susceptibility result
    /// @property
    bse::block_matrix<Scalar>& GetChi0q() {
      return this->chi0_q;
    }

    /// Get chi_q_rpa - RPA susceptibility result
    /// @property
    bse::block_matrix<Scalar>& GetChiq_RPA() {
      return this->chi_q_rpa;
    }

    /// Get chi_q_rrpa - RRPA susceptibility result
    /// @property
    bse::block_matrix<Scalar>& GetChiq_RRPA() {
      return this->chi_q_rrpa;
    }

    /// Get chi_q_scl - SCL susceptibility result
    /// @property
    bse::block_matrix<Scalar>& GetChiq_SCL() {
      return this->chi_q_scl;
    }

    /// Get Iq - Irreducible vertex result
    /// @property
    bse::block_matrix<Scalar>& GetIq() {
      return this->Iq;
    }

    /**
     * @brief Print memory usage summary
     */
    void Summary() {
      memInfo.summary();
    }

    /**
     * @brief Calculate chi0_q susceptibility
     *
     * Computes the momentum-dependent bare susceptibility:
     * chi0_q = chi0_loc + (1/beta) * sum_w (X0_q - X0_loc)
     */
    void CalcChi0q() {
      require(X0_q, "X0_q");
      require(X0_Loc, "X0_Loc");
      require(chi0_loc, "chi0_loc");

      // Calculate chi0_q by summing over frequencies
      auto X0_q_sum = Sumfreq_B(*X0_q - *X0_Loc);
      this->chi0_q = *chi0_loc + (1.0 / std::complex<double>(beta)) * X0_q_sum;

      // for RPA, RRPA
      // this->chi0_q_ptr = &this->chi0_q;
    }

    /**
     * @brief Calculate BSE susceptibility
     *
     * Computes the full BSE susceptibility including vertex corrections:
     * chi_q = chi_loc + (1/beta) * sum_w Z_q
     * where Z_q = X_loc * [1 - P_q * X_loc]^-1 * P_q * X_loc
     * and P_q = X0_loc^-1 - X0_q^-1
     */
    void CalcChiq_BSE() {
      require(X0_q, "X0_q");
      require(X0_Loc, "X0_Loc");
      require(X_loc, "X_loc");
      require(chi_loc, "chi_loc");

      // Calculate P_q = X0_loc^-1 - X0_q^-1
      bse::block_matrix<Scalar> Pq(matrixInfo_B);
      Pq = bse::inverse(*X0_Loc) - bse::inverse(*X0_q);
      memInfo.add(Pq, "Pq");

      // Calculate Z_q = X_loc * [1 - P_q * X_loc]^-1 * P_q * X_loc
      bse::block_matrix<Scalar> Zq(matrixInfo_A);
      bse::block_matrix<Scalar> Identity(matrixInfo_A);
      Identity = bse::identity_like(Identity);

      // First compute X_loc * P_q
      bse::block_matrix<Scalar> Xloc_Pq = Prod_A_B(*X_loc, Pq);
      memInfo.add(Xloc_Pq, "Xloc_Pq");

      // Then compute Z_q
      Zq = (bse::inverse(Identity - Xloc_Pq) - Identity) * *X_loc;
      memInfo.add(Zq, "Zq");

      // Calculate final chi_q_bse by summing over frequencies
      this->chi_q_bse = *chi_loc + (1.0 / std::complex<double>(beta)) * Sumfreq_A(Zq);

      // Calculate irreducible vertex I_q = chi_q_bse^-1 - chi_loc^-1
      this->Iq = bse::inverse(*chi_loc) - bse::inverse(chi_q_bse);
    }

    /**
     * @brief Calculate RPA susceptibility
     *
     * Computes the RPA susceptibility:
     * chi_q_rpa = [1 - chi0_q * gamma0]^-1 * chi0_q
     */
    void CalcChiq_RPA() {
      require(gamma0, "gamma0");
      require(chi0_q_in, "chi0_q");
      // require(chi0_q_ptr, "chi0_q_ptr");
      // if (!chi0_q_ptr){
      //   CalcChi0q();
      // }

      bse::block_matrix<std::complex<double> > Identity(matrixInfo_C);
      Identity = bse::identity_like(Identity);
      this->chi_q_rpa = bse::inverse(Identity - *chi0_q_in * *gamma0) * *chi0_q_in;
    }

    /**
     * @brief Calculate RRPA susceptibility
     *
     * Computes the renormalized RPA susceptibility:
     * chi_q_rrpa = [1 - chi0_q * U_eff]^-1 * chi0_q
     * where U_eff = chi0_loc^-1 - chi_loc^-1
     */
    void CalcChiq_RRPA() {
      require(chi0_loc, "chi0_loc");
      require(chi_loc, "chi_loc");
      require(chi0_q_in, "chi0_q");
      // require(chi0_q_ptr, "chi0_q_ptr");
      // if (!chi0_q_ptr){
      //   CalcChi0q();
      // }

      bse::block_matrix<std::complex<double> > Identity(matrixInfo_C);
      Identity = bse::identity_like(Identity);

      // Calculate effective interaction U_eff
      auto Ueff = bse::inverse(*chi0_loc) - bse::inverse(*chi_loc);
      this->chi_q_rrpa = bse::inverse(Identity - *chi0_q_in * Ueff) * *chi0_q_in;
    }

    /**
     * @brief Calculate SCL susceptibility
     *
     * Computes the self-consistent ladder susceptibility:
     * chi_q_scl = Phi_sum * [1 - lambda_q]^-1 * Phi_sum
     * where lambda_q = sum_w (Phi * Lambda_q * Phi)
     * and Lambda_q = X0_loc^-1 - X0_q^-1
     */
    void CalcChiq_SCL() {
      require(X0_q, "X0_q");
      require(X0_Loc, "X0_Loc");
      require(Phi, "Phi");
      require(Phi_Sum, "Phi_Sum");

      // Calculate Lambda_q = X0_loc^-1 - X0_q^-1
      bse::block_matrix<Scalar> Lambdaq(matrixInfo_B);
      Lambdaq = bse::inverse(*X0_Loc) - bse::inverse(*X0_q);
      memInfo.add(Lambdaq, "Lambdaq");

      // Calculate lambda_q by summing over frequencies
      bse::block_matrix<Scalar> lambda_q;
      lambda_q = Sumfreq_B(*Phi * Lambdaq * *Phi);  // matrix C
      memInfo.add(lambda_q, "lambda_q");

      // Setup identity matrix
      bse::block_matrix<Scalar> Identity;
      Identity = bse::identity_like(lambda_q);

      // Calculate final chi_q_scl
      this->chi_q_scl = *Phi_Sum * bse::inverse(Identity - lambda_q) * *Phi_Sum;
      memInfo.add(chi_q_scl, "chi_q_scl");

      // Calculate irreducible vertex I_q
      bse::block_matrix<Scalar> chi_minv_sqrt;
      chi_minv_sqrt = bse::inverse(*Phi_Sum);
      memInfo.add(chi_minv_sqrt, "chi_minv_sqrt");
      this->Iq = chi_minv_sqrt * lambda_q * chi_minv_sqrt;
      memInfo.add(this->Iq, "Iq");
    }


    // bool DiagonalizeChiqSum(bse::block_matrix<Scalar> &chiq_sum, std::vector<double> &eigen_values) {
    //   eigen_values.resize(matrixInfo.size(), 0);
    //   return true;
    // }

    /**
     * @brief Sum over frequencies for A-type matrix
     *
     * Converts A-type matrix (block: [block], inner: [in, iw])
     * to C-type matrix (block: [block], inner: [in])
     * by summing over frequency indices
     *
     * @param bm Input A-type block matrix
     * @return C-type block matrix result
     */
    bse::block_matrix<Scalar> Sumfreq_A(const bse::block_matrix<Scalar> &bm) {
      // matrix-A
      //   block: [block], inner: [in, iw]
      // matrix-C
      //   block: [block], inner: [in]
      assert (isBlockMatrix_A(bm));

      int nb = matrixInfo_A.size();
      int nw = int(matrixInfo_A[0] / matrixInfo_C[0]);

      using matrix_t = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

      bse::block_matrix<Scalar> chiq_sum(matrixInfo_C);
      chiq_sum = zeros_like(chiq_sum);

      // Parallelize over block pairs
#pragma omp parallel for default(none) firstprivate(nb, nw) shared(chiq_sum, bm)
      for (int num = 0; num < nb*nb; ++num) {
        int i = int(num / nb);
        int j = num % nb;
        if (!bm.exists(i, j)) {
          continue;
        }
        matrix_t mat_A = bm(i, j);
        int dim1 = matrixInfo_C[i];
        int dim2 = matrixInfo_C[j];
        matrix_t m(dim1, dim2);
        m.setZero();

        // Sum over frequency indices
        for (int in1=0; in1<dim1; in1++){
          for (int in2=0; in2<dim2; in2++){
            for (int iw1=0; iw1<nw; iw1++){
              for (int iw2=0; iw2<nw; iw2++){
                m(in1, in2) += mat_A(in1 * nw + iw1, in2 * nw + iw2);
              }
            }
          }
        }
        chiq_sum.assign(i, j, m);
      }
      return chiq_sum;
    }

    /**
     * @brief Sum over frequencies for B-type matrix
     *
     * Converts B-type matrix (block: [iw, block], inner: [in])
     * to C-type matrix (block: [block], inner: [in])
     * by summing over frequency indices
     *
     * @param bm Input B-type block matrix
     * @return C-type block matrix result
     */
    bse::block_matrix<Scalar> Sumfreq_B(const bse::block_matrix<Scalar> &bm) {
      // matrix-B
      //   block: [iw, block], inner: [in]
      // matrix-C
      //   block: [block], inner: [in]
      assert (isBlockMatrix_B(bm));

      int nb = matrixInfo_C.size();
      int nw = int(matrixInfo_B.size() / nb);

      using matrix_t = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

      bse::block_matrix<Scalar> chiq_sum(matrixInfo_C);  // return
      chiq_sum = zeros_like(chiq_sum);

      // Parallelize over block pairs
#pragma omp parallel for default(none) firstprivate(nb, nw) shared(chiq_sum, bm)
      for (int num = 0; num < nb*nb; ++num) {
        int i = int(num / nb);
        int j = num % nb;
        matrix_t m(matrixInfo_C[i], matrixInfo_C[j]);
        m.setZero();

        // Sum over frequency indices
        for (int iw=0; iw<nw; iw++){
          int i2 = i + nb * iw;
          int j2 = j + nb * iw;
          if (!bm.exists(i2, j2)) {
            continue;
          }
          m += bm(i2, j2);
        }
        chiq_sum.assign(i, j, m);
      }
      return chiq_sum;
    }

  private:
    using block_type = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    std::vector<int> matrixInfo_A, matrixInfo_B, matrixInfo_C;  // Matrix dimension info
    bse::memory_info memInfo;  // Memory usage tracker

    double beta;  // Inverse temperature

    // Input matrices
    bse::block_matrix<Scalar> const* X0_Loc = nullptr;  // Local bare susceptibility
    bse::block_matrix<Scalar> const* X_loc = nullptr;   // Local full vertex
    bse::block_matrix<Scalar> const* X0_q = nullptr;    // q-dependent bare susceptibility
    bse::block_matrix<Scalar> const* Phi = nullptr;     // Vertex function
    //bse::block_matrix<Scalar> Phi_L;                  // Left vertex function
    //bse::block_matrix<Scalar> Phi_R;                  // Right vertex function
    bse::block_matrix<Scalar> const* Phi_Sum = nullptr; // Summed vertex function
    bse::block_matrix<Scalar> const* chi_loc = nullptr; // Local physical susceptibility
    bse::block_matrix<Scalar> const* chi0_loc = nullptr;// Local bare physical susceptibility
    bse::block_matrix<Scalar> const* gamma0 = nullptr;  // Bare vertex

    bse::block_matrix<Scalar> const* chi0_q_in = nullptr;  // for RPA, RRPA
    // bse::block_matrix<Scalar> const* chi0_q_ptr = nullptr;  // Pointer to chi0_q (set in CalcChi0q)

    // Output matrices
    bse::block_matrix<Scalar> chi0_q;      // Bare susceptibility
    bse::block_matrix<Scalar> chi_q_bse;   // BSE susceptibility
    bse::block_matrix<Scalar> chi_q_rpa;   // RPA susceptibility
    bse::block_matrix<Scalar> chi_q_rrpa;  // RRPA susceptibility
    bse::block_matrix<Scalar> chi_q_scl;   // SCL susceptibility
    bse::block_matrix<Scalar> Iq;          // Irreducible vertex (BSE, SCL)

    /**
     * @brief Check if required input has been set
     * @param bm Block matrix pointer to check
     * @param s Name of matrix for error message
     * @throws runtime_error if matrix is not set
     */
    void require(bse::block_matrix<Scalar> const* bm, const std::string &s) {
      if (!bm){
        std::stringstream ss;
        ss << "'" << s << "' must be set before calling calc.";
        throw std::runtime_error(ss.str());
      }
    }

    /**
     * @brief Convert B-type matrix to A-type matrix
     *
     * Converts between matrix formats:
     * B-type: block [iw, block], inner [in]
     * A-type: block [block], inner [in, iw]
     *
     * @param bm_B Input B-type block matrix
     * @return A-type block matrix
     */
    bse::block_matrix<Scalar> Convert_B2A(const bse::block_matrix<Scalar> &bm_B) {
      // matrix-A
      //   block: [block], inner: [in, iw]
      // matrix-B
      //   block: [iw, block], inner: [in]
      assert (isBlockMatrix_B(bm_B));

      using matrix_t = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
      int nb = matrixInfo_C.size();
      int nw = int(matrixInfo_B.size() / nb);

      bse::block_matrix<Scalar> bm_A(matrixInfo_A);  // return
      bm_A = zeros_like(bm_A);

      // Loop over blocks and frequencies
      for (int b1=0; b1<nb; b1++){
        for (int b2=0; b2<nb; b2++){
          int dim1 = matrixInfo_C[b1];
          int dim2 = matrixInfo_C[b2];
          matrix_t mat_A(dim1*nw, dim2*nw);
          mat_A.setZero();
          for (int iw1=0; iw1<nw; iw1++){
            for (int iw2=0; iw2<nw; iw2++){
              if (!bm_B.exists(iw1 * nb + b1, iw2 * nb + b2)){
                continue;
              }
              matrix_t mat_B = bm_B(iw1 * nb + b1, iw2 * nb + b2);
              // Copy elements with reordered indices
              for (int in1=0; in1<dim1; in1++){
                for (int in2=0; in2<dim2; in2++){
                  mat_A(in1 * nw + iw1, in2 * nw + iw2) = mat_B(in1, in2);
                }
              }
            }
          }
          if (!mat_A.isZero(0)){
            bm_A.assign(b1, b2, mat_A);
          }
        }
      }
      return bm_A;
    }

    /**
     * @brief Convert A-type matrix to B-type matrix
     *
     * Converts between matrix formats:
     * A-type: block [block], inner [in, iw]
     * B-type: block [iw, block], inner [in]
     *
     * @param bm_A Input A-type block matrix
     * @return B-type block matrix
     */
    bse::block_matrix<Scalar> Convert_A2B(const bse::block_matrix<Scalar> &bm_A) {
      // matrix-A
      //   block: [block], inner: [in, iw]
      // matrix-B
      //   block: [iw, block], inner: [in]
      assert (isBlockMatrix_A(bm_A));

      using matrix_t = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
      int nb = matrixInfo_C.size();
      int nw = int(matrixInfo_B.size() / nb);

      bse::block_matrix<Scalar> bm_B(matrixInfo_B);  // return
      bm_B = zeros_like(bm_B);

      // Loop over blocks and frequencies
      for (int b1=0; b1<nb; b1++){
        for (int b2=0; b2<nb; b2++){
          if (!bm_A.exists(b1, b2)){
            continue;
          }
          matrix_t mat_A = bm_A(b1, b2);
          int dim1 = matrixInfo_C[b1];
          int dim2 = matrixInfo_C[b2];
          for (int iw1=0; iw1<nw; iw1++){
            for (int iw2=0; iw2<nw; iw2++){
              matrix_t mat_B(dim1, dim2);
              mat_B.setZero();
              // Copy elements with reordered indices
              for (int in1=0; in1<dim1; in1++){
                for (int in2=0; in2<dim2; in2++){
                  mat_B(in1, in2) = mat_A(in1 * nw + iw1, in2 * nw + iw2);
                }
              }
              if (!mat_B.isZero(0)){
                bm_B.assign(iw1 * nb + b1, iw2 * nb + b2, mat_B);
              }
            }
          }
        }
      }
      return bm_B;
    }

    /**
     * @brief Calculate product of A-type and B-type matrices
     *
     * Optimized matrix multiplication that takes advantage of
     * the block structure of B-type matrices
     *
     * @param bm_A Input A-type block matrix
     * @param bm_B Input B-type block matrix
     * @return A-type block matrix result
     */
    bse::block_matrix<Scalar> Prod_A_B(const bse::block_matrix<Scalar> &bm_A, const bse::block_matrix<Scalar> &bm_B) {
      assert (isBlockMatrix_A(bm_A));
      assert (isBlockMatrix_B(bm_B));

      // Convert A to B-type, multiply, then convert back to A-type
      // This is more efficient since B-type matrices are block diagonal in frequency
      bse::block_matrix<Scalar> bm_A_inB = Convert_A2B(bm_A);
      bse::block_matrix<Scalar> bm_AB_inB = bm_A_inB * bm_B;
      bse::block_matrix<Scalar> bm_AB_inA = Convert_B2A(bm_AB_inB);

      return bm_AB_inA;
    }

    /**
     * @brief Check if matrix has A-type structure
     * @param bm Block matrix to check
     * @return true if matrix has A-type structure
     */
    bool isBlockMatrix_A(const bse::block_matrix<Scalar> &bm) {
      return CheckBlockMatrix(bm, matrixInfo_A);
    }

    /**
     * @brief Check if matrix has B-type structure
     * @param bm Block matrix to check
     * @return true if matrix has B-type structure
     */
    bool isBlockMatrix_B(const bse::block_matrix<Scalar> &bm) {
      return CheckBlockMatrix(bm, matrixInfo_B);
    }

    /**
     * @brief Check if matrix has C-type structure
     * @param bm Block matrix to check
     * @return true if matrix has C-type structure
     */
    bool isBlockMatrix_C(const bse::block_matrix<Scalar> &bm) {
      return CheckBlockMatrix(bm, matrixInfo_C);
    }

    /**
     * @brief Check if block matrix has expected structure
     * @param bm_A Block matrix to check
     * @param _matrix_info Expected matrix info structure
     * @return true if matrix has expected structure
     */
    bool CheckBlockMatrix(const bse::block_matrix<Scalar> &bm_A, const std::vector<int> &_matrix_info) {
      bse::block_matrix<Scalar> bm_empty(_matrix_info);
      return bm_A.has_same_block_structure(bm_empty);
    }

  };
}
