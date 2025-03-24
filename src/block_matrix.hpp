#pragma once

// Standard library includes
#include <vector>
#include <numeric>
#include <memory>
#include <set>

// OpenMP support if available
#ifdef _OPENMP
#include <omp.h>
#endif

// Eigen library includes
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

/**
 * @brief Block matrix library
 * 
 * This namespace contains classes and functions for working with block matrices,
 * which are matrices divided into submatrices (blocks) for efficient operations.
 */
namespace bse {

  /**
   * @brief Block matrix class
   * 
   * This class implements a block matrix structure where the matrix is divided into
   * blocks that can be manipulated independently. The blocks are stored sparsely,
   * meaning only non-zero blocks are stored in memory.
   * 
   * @tparam Scalar Scalar type used for matrix elements (e.g. double, float)
   */
  template<typename Scalar>
  class block_matrix {
  public:
    /**
     * @brief Default constructor
     * 
     * If an object is constructed via default constructor,
     * it must be initialized via copy before any use.
     * Creates an empty block matrix with no blocks.
     */
    block_matrix()
            : block_matrix(std::vector<int>()) {
    }

    /**
     * @brief Construct a block_matrix object with the given block structure
     * 
     * All the blocks are initialized to 0. The matrix is divided into blocks
     * according to the dimensions specified in dim_blocks.
     * 
     * @param dim_blocks Vector containing dimensions of each block
     */
    block_matrix(const std::vector<int> &dim_blocks)
            : dim_blocks_(dim_blocks), left_top_corner_blocks_(),
              dim_(std::accumulate(dim_blocks.begin(), dim_blocks.end(), 0)),
              p_blocks_() {
      // Reserve space for all possible blocks
      p_blocks_.reserve(dim_blocks.size() * dim_blocks.size());
      
      int left_top_corner = 0;
      int offset = 0;
      int bl = 0;
      
      // Initialize block structure
      for (int bl_row = 0; bl_row < dim_blocks.size(); ++bl_row) {
        for (int bl_col = 0; bl_col < dim_blocks.size(); ++bl_col) {
          p_blocks_.push_back(std::unique_ptr<block_type>());
          left_top_corner_blocks_.push_back(offset);
          offset += dim_blocks[bl_col];
          bl++;
        }
        offset = dim_blocks[bl_row] * dim_blocks.size();
      }
    }

    /**
     * @brief Copy constructor
     * 
     * Creates a deep copy of another block matrix.
     * 
     * @param other Block matrix to be copied
     */
    block_matrix(const block_matrix<Scalar> &other) : block_matrix() {
      *this = other;
    }

    /**
     * @brief Check if a block exists (is non-zero)
     * 
     * @param block_row Row index of block
     * @param block_col Column index of block
     * @return true if block exists, false otherwise
     */
    bool exists(int block_row, int block_col) const {
      range_check(block_row);
      range_check(block_col);
      return bool(p_blocks_[flatten_index(block_row, block_col)]);
    }

    /**
     * @brief Get linear dimension of the whole matrix
     * 
     * Returns the total size of the matrix (sum of block dimensions)
     * 
     * @return Total dimension
     */
    int dim() const {
      return dim_;
    }

    /**
     * @brief Get number of blocks
     * 
     * Returns the number of blocks in each row/column
     * 
     * @return Number of blocks
     */
    int num_blocks() const {
      return dim_blocks_.size();
    }

    /**
     * @brief Get dimension of a block
     * 
     * Returns the size of the specified block
     * 
     * @param block Index of block
     * @return Dimension of the block
     */
    int dim_block(int block) const {
      range_check(block);
      return dim_blocks_[block];
    }

    /**
     * @brief Comparison operator
     * 
     * Return true if the positions of the existing blocks match between the two objects
     * and the existing blocks contain the same values.
     * 
     * @param other Object to compare with
     * @return true if equal, false otherwise
     */
    bool operator==(const block_matrix<Scalar> &other) const {
      // Check if block structures match
      if (this->dim_blocks_ != other.dim_blocks_ || this->left_top_corner_blocks_ != other.left_top_corner_blocks_) {
        return false;
      }
      
      // Compare each block
      for (int i = 0; i < num_blocks(); ++i) {
        for (int j = 0; j < num_blocks(); ++j) {
          if (exists(i, j) != other.exists(i, j)) {
            return false;
          }
          //compare values if block exists
          if (exists(i, j) && this->operator()(i, j) != other(i, j)) {
            return false;
          }
        }
      }
      return true;
    }

    /**
     * @brief Get size of object in bytes
     * 
     * Calculates total memory usage of non-zero blocks
     * 
     * @return Size in bytes
     */
    std::size_t get_bytes() const {
      std::size_t nelem = 0;
      for (const auto &b : p_blocks_) {
        if (b) {
          nelem += b->rows() * b->cols();
        }
      }
      return sizeof(Scalar) * nelem;
    }

    /**
     * @brief Check if two objects have same block structure
     * 
     * Compares only the block dimensions, not the actual contents
     * 
     * @param other Object to compare with
     * @return true if same structure, false otherwise
     */
    bool has_same_block_structure(const block_matrix<Scalar> &other) const {
      return dim_blocks_ == other.dim_blocks_;
    }

    /**
     * @brief Copy assignment operator
     * 
     * Performs deep copy of blocks. All existing blocks in target are released
     * and replaced with copies from source.
     * 
     * @param other Object to copy from
     * @return Reference to this object
     */    
    block_matrix<Scalar> &operator=(const block_matrix<Scalar> &other) {
      dim_blocks_ = other.dim_blocks_;
      left_top_corner_blocks_ = other.left_top_corner_blocks_;
      dim_ = other.dim_;

      //release all existing blocks
      p_blocks_.resize(0);

      //perform deep copy
      p_blocks_.reserve(dim_blocks_.size() * dim_blocks_.size());
      for (int bl = 0; bl < dim_blocks_.size() * dim_blocks_.size(); ++bl) {
        if (other.p_blocks_[bl]) {
          p_blocks_.push_back(
                  std::unique_ptr<block_type>(new block_type(*other.p_blocks_[bl]))
          );
        } else {
          p_blocks_.push_back(std::unique_ptr<block_type>());
        }
      }

      return *this;
    }

    /**
     * @brief Assign new values to a block
     * 
     * Updates or creates a block with new values
     * 
     * @param block_row Row index of block
     * @param block_col Column index of block
     * @param new_block Matrix of new values
     */
    template<typename Derived>
    void assign(int block_row, int block_col, const Eigen::MatrixBase<Derived> &new_block) {
      range_check(block_row);
      range_check(block_col);

      // Check dimensions match
      if (new_block.rows() != dim_blocks_[block_row] || new_block.cols() != dim_blocks_[block_col]) {
        throw std::runtime_error("Size of new block is not consistent with the block structure.");
      }

      // Update existing block or create new one
      if (exists(block_row, block_col)) {
        *p_blocks_[flatten_index(block_row, block_col)] = new_block;
      } else {
        p_blocks_[flatten_index(block_row, block_col)].reset(new block_type(new_block));
      }
    }

    /**
     * @brief Get read-only view of a block
     * 
     * @param block_row Row index of block
     * @param block_col Column index of block
     * @return Const reference to block matrix
     * @throws std::runtime_error if block does not exist
     */
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &
    operator()(int block_row, int block_col) const {
      if (!exists(block_row, block_col)) {
        throw std::runtime_error("Block does not exists");
      }
      return *p_blocks_[flatten_index(block_row, block_col)];
    };

    /**
     * @brief Convert block matrix to regular matrix
     * 
     * Creates a dense matrix representation of the block matrix
     * 
     * @return Dense matrix representation
     */
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    convert_to_matrix() const {
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> m(dim(), dim());
      m.setZero();
      int num_block = num_blocks();
      
      // Parallel loop over all blocks
#pragma omp parallel for default(none) firstprivate(num_block, left_top_corner_blocks_, p_blocks_) shared(m)
      for (int num = 0; num < num_block*num_block; ++num) {
        int i = int(num / num_block);
        int j = num % num_block;
        if (!exists(i, j)) {
          continue;
        }
        // Copy block to appropriate location in dense matrix
        m.block(left_top_corner_blocks_[i], left_top_corner_blocks_[j], dim_block(i),
                dim_block(j)) = *p_blocks_[flatten_index(i, j)];
      }

      return m;
    };

  private:
    /** @brief Type alias for block matrix */
    using block_type = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    /**
     * @brief Check if block index is in valid range
     * 
     * @param block_row_or_col Block index to check
     * @throws std::runtime_error if index is out of range
     */
    void range_check(int block_row_or_col) const {
      if (block_row_or_col < 0 || num_blocks() <= block_row_or_col) {
        throw std::runtime_error(
            "block_row_or_col is out of range " +
            std::to_string(block_row_or_col) +
            " num_blocks = " + std::to_string(num_blocks())
          );
      }
    }

    /**
     * @brief Convert 2D block indices to flattened 1D index
     * 
     * Maps (row,col) coordinates to linear index in storage array
     * 
     * @param block_row Row index
     * @param block_col Column index
     * @return Flattened index
     */
    int flatten_index(int block_row, int block_col) const {
      return block_row * num_blocks() + block_col;
    }

    /** @brief Dimensions of blocks */
    std::vector<int> dim_blocks_;

    /** @brief Left-top corner coordinates of blocks */
    std::vector<int> left_top_corner_blocks_;

    /** @brief Linear dimension of the matrix */
    int dim_;

    /** @brief Pointers to blocks in row major order */
    std::vector<std::unique_ptr<block_type>> p_blocks_;
  };

  namespace detail {

    /**
     * @brief Create a list of blocks to be combined
     * 
     * Groups blocks that are connected according to the provided functor
     * 
     * @tparam F Type of functor
     * @param connected Functor that returns whether two blocks should be in same group
     * @param num_blocks Number of blocks
     * @return Vector of block groups, each containing indices of blocks in same group
     */
    template<typename F>
    std::vector<std::vector<int>>
    gather_blocks(const F &connected, int num_blocks) {
      // Track unprocessed blocks
      std::set<int> remains;
      for (int bl = 0; bl < num_blocks; ++bl) {
        remains.insert(bl);
      }

      std::vector<std::vector<int>> bl_groups;

      int gr = 0;
      // Process blocks until all are assigned to groups
      while (remains.size() > 0) {
        // Start new group with first remaining block
        bl_groups.push_back(std::vector<int>{*remains.begin()});
        remains.erase(*remains.begin());
        
        //find blocks that belong to the current group
        while (true) {
          bool updated = false;
          for (auto bl : bl_groups[gr]) {
            for (auto it = remains.begin(); it != remains.end(); ++it) {
              assert (bl != *it);
              if (connected(bl, *it)) {
                updated = true;
                bl_groups.back().push_back(*it);
                remains.erase(*it);
                break;
              }
            }
            if (updated) {
              //We start over a loop each time bl_groups and remains have been updated.
              break;
            }
          }
          //if remains has not been updated, exit the loop for a new group.
          if (!updated) {
            break;
          }
        }
        ++gr;
      }

      return bl_groups;
    }
  }

  /**
   * @brief Multiply scalar by block matrix
   * 
   * Performs element-wise multiplication of scalar with matrix
   * 
   * @tparam Scalar Scalar type
   * @param s Scalar value
   * @param m Block matrix
   * @return Result of multiplication
   */
  template<typename Scalar>
  block_matrix<Scalar> operator*(Scalar s, block_matrix<Scalar> m) {
    int nb = m.num_blocks();
    using matrix_t = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    auto bm_r = zeros_like(m);
    
    // Parallel loop over all blocks
#pragma omp parallel for default(none) firstprivate(nb, s) shared(m, bm_r)
    for (int num = 0; num < nb*nb; num++) {
        int i = int(num / nb);
        int j = num % nb;
        if (!m.exists(i, j)) {
          continue;
        }
        matrix_t bm(m.dim_block(i), m.dim_block(j));
        bm.setZero();
        bm = m(i, j);
        bm *= s;
        bm_r.assign(i, j, bm);
    }

    return bm_r;
  }

  /**
   * @brief Add two block matrices
   * 
   * Performs element-wise addition of matrices
   * 
   * @tparam Scalar Scalar type
   * @param bm First matrix
   * @param bm2 Second matrix
   * @return Sum of matrices
   * @throws std::runtime_error if matrices have incompatible block structure
   */
  template<typename Scalar>
  block_matrix<Scalar> operator+(const block_matrix<Scalar> &bm, const block_matrix<Scalar> &bm2) {
    if (!bm.has_same_block_structure(bm2)) {
      throw std::runtime_error("Error in operator+: block structure is not compatible");
    }

    int nb = bm.num_blocks();
    using matrix_t = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    auto bm_r = zeros_like(bm);
    
    // Parallel loop over all blocks
#pragma omp parallel for default(none) firstprivate(nb) shared(bm, bm2, bm_r)
    for (int num = 0; num < nb*nb; num++) {
      int i = int(num / nb);
      int j = num % nb;
      if (!bm.exists(i, j) && !bm2.exists(i, j)) {
        continue;
      }
      matrix_t m(bm.dim_block(i), bm.dim_block(j));
      m.setZero();
      if (bm.exists(i, j)) {
        m += bm(i, j);
      }
      if (bm2.exists(i, j)) {
        m += bm2(i, j);
      }
      bm_r.assign(i, j, m);
    }

    return bm_r;
  }

  /**
   * @brief Subtract two block matrices
   * 
   * Performs element-wise subtraction of matrices
   * 
   * @tparam Scalar Scalar type
   * @param bm First matrix
   * @param bm2 Second matrix
   * @return Difference of matrices
   * @throws std::runtime_error if matrices have incompatible block structure
   */
  template<typename Scalar>
  block_matrix<Scalar> operator-(const block_matrix<Scalar> &bm, const block_matrix<Scalar> &bm2) {
    if (!bm.has_same_block_structure(bm2)) {
      throw std::runtime_error("Error in operator+: block structure is not compatible");
    }

    int nb = bm.num_blocks();
    using matrix_t = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    auto bm_r = zeros_like(bm);

    // Parallel loop over all blocks
#pragma omp parallel for default(none) firstprivate(nb) shared(bm, bm2, bm_r)
    for (int num = 0; num < nb * nb; num++) {
      int i = int(num / nb);
      int j = num % nb;
      if (!bm.exists(i, j) && !bm2.exists(i, j)) {
        continue;
      }
      matrix_t m(bm.dim_block(i), bm.dim_block(j));
      m.setZero();
      if (bm.exists(i, j)) {
        m += bm(i, j);
      }
      if (bm2.exists(i, j)) {
        m -= bm2(i, j);
      }
      bm_r.assign(i, j, m);
    }

    return bm_r;
  }

  /**
   * @brief Compute inverse of block matrix
   * 
   * Inverts the matrix by first grouping connected blocks and inverting each group
   * 
   * @tparam Scalar Scalar type
   * @param bm Input matrix
   * @return Inverse matrix
   */
  template<typename Scalar>
  block_matrix<Scalar> inverse(const block_matrix<Scalar> &bm) {
    using matrix_type = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    auto inv_bm = zeros_like(bm);

    // Function to determine if blocks are connected
    auto connected = [&bm](int i, int j) {
      return bm.exists(i, j) || bm.exists(j, i);
    };
    
    // Group connected blocks
    auto bl_groups = detail::gather_blocks(connected, bm.num_blocks());

    int nb = bl_groups.size();
    
    // Parallel loop over block groups
#pragma omp parallel for default(none) firstprivate(nb, bl_groups) shared(bm, inv_bm)
    for(int num = 0; num <nb; num++){
      auto gr = bl_groups[num];
      int matrix_size = 0;
      std::vector<int> left_upper_corner(gr.size());
      
      // Calculate dimensions and offsets
      for (int bl = 0; bl < gr.size(); ++bl) {
        left_upper_corner[bl] = matrix_size;
        matrix_size += bm.dim_block(gr[bl]);
      }

      // Create super-block from group
      matrix_type super_block(matrix_size, matrix_size);
      super_block.setZero();

      for (int bl = 0; bl < gr.size(); ++bl) {
        for (int bl2 = 0; bl2 < gr.size(); ++bl2) {
          if (bm.exists(gr[bl], gr[bl2])) {
            super_block.block(
                    left_upper_corner[bl],
                    left_upper_corner[bl2],
                    bm.dim_block(gr[bl]),
                    bm.dim_block(gr[bl2])
            ) = bm(gr[bl], gr[bl2]);
          }
        }
      }

      // Invert super-block
      matrix_type inv_super_block = super_block.inverse();

      // Extract blocks from inverted super-block
      for (int bl = 0; bl < gr.size(); ++bl) {
        for (int bl2 = 0; bl2 < gr.size(); ++bl2) {
          inv_bm.assign(
                  gr[bl],
                  gr[bl2],
                  inv_super_block.block(
                          left_upper_corner[bl],
                          left_upper_corner[bl2],
                          bm.dim_block(gr[bl]),
                          bm.dim_block(gr[bl2])
                  )
          );
        }
      }
    }

    return inv_bm;
  }

  /**
   * @brief Multiply two block matrices
   * 
   * Performs matrix multiplication of block matrices
   * 
   * @tparam Scalar Scalar type
   * @param bm First matrix
   * @param bm2 Second matrix
   * @return Product of matrices
   * @throws std::runtime_error if matrices have incompatible block structure
   */
  template<typename Scalar>
  block_matrix<Scalar> operator*(const block_matrix<Scalar> &bm, const block_matrix<Scalar> &bm2) {
    if (!bm.has_same_block_structure(bm2)) {
      throw std::runtime_error("Error in operator*: block structure is not compatible");
    }

    int nb = bm.num_blocks();
    using matrix_t = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    auto bm_r = zeros_like(bm);

    // Parallel loop over result blocks
#pragma omp parallel for default(none) firstprivate(nb) shared(bm, bm2, bm_r)
    for (int num = 0; num < nb * nb; ++num) {
      int i = int(num / nb);
      int j = num % nb;
      matrix_t m(bm.dim_block(i), bm.dim_block(j));
      m.setZero();

      bool non_zero = false;
      // Multiply and sum intermediate blocks
      for (int k = 0; k < nb; ++k) {
        if (bm.exists(i, k) && bm2.exists(k, j)) {
          non_zero = true;
          m += bm(i, k) * bm2(k, j);
        }
      }

      if (non_zero) {
        bm_r.assign(i, j, m);
      }
    }
    return bm_r;
  }

  /**
   * @brief Create zero block matrix with same structure as input
   * 
   * @tparam Scalar Scalar type
   * @param m Input matrix
   * @return Zero matrix with same block structure
   */
  template<typename Scalar>
  block_matrix<Scalar> zeros_like(const block_matrix<Scalar> &m) {
    std::vector<int> dim_blocks;
    dim_blocks.reserve(m.num_blocks());
    for (int b = 0; b < m.num_blocks(); ++b) {
      dim_blocks.push_back(m.dim_block(b));
    }
    return block_matrix<Scalar>(dim_blocks);
  }

  /**
   * @brief Create identity block matrix with same structure as input
   * 
   * Creates a block diagonal matrix with identity matrices on diagonal
   * 
   * @tparam Scalar Scalar type
   * @param m Input matrix
   * @return Identity matrix with same block structure
   */
  template<typename Scalar>
  block_matrix<Scalar> identity_like(const block_matrix<Scalar> &m) {
    auto m_new = zeros_like(m);
    for (int b = 0; b < m.num_blocks(); ++b) {
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> new_block(m.dim_block(b), m.dim_block(b));
      new_block.setIdentity();
      m_new.assign(b, b, new_block);
    }
    return m_new;
  }

}
