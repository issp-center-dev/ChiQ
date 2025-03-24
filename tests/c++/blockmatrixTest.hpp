#pragma once

#include "gtest.h"
#include "../src/block_matrix.hpp"

/**
 * Test function for block matrix assignment operations
 * Tests assignment, copying, and element access for block matrices
 * @tparam Scalar Type of matrix elements
 */
template<typename Scalar>
void test_assignment() {
  // Create block matrix with dimensions {1,2,3}
  bse::block_matrix<Scalar> m(std::vector<int>{1,2,3});

  // Create 3x3 test matrix
  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> bm(3,3);

  // Fill test matrix with values i + 10*j
  for (int j=0; j<3; ++j) {
    for (int i=0; i<3; ++i) {
      bm(i,j) = i+10*j;
    }
  }

  // Assign test matrix to block (2,2)
  m.assign(2, 2, bm);

  // Verify block exists and has correct size
  ASSERT_TRUE(m.exists(2, 2));
  ASSERT_EQ(sizeof(Scalar)*3*3, m.get_bytes());

  // Get view of block (2,2)
  auto bm_view = m(2,2);

  // Verify block elements match original matrix
  for (int j=0; j<2; ++j) {
    for (int i = 0; i < 2; ++i) {
      ASSERT_TRUE(bm_view(i,j) == bm(i,j));
    }
  }

  // Test copy constructor and equality
  auto m_copy = m;
  ASSERT_TRUE(m_copy == m);
  ASSERT_TRUE(m_copy.convert_to_matrix() == m.convert_to_matrix());
}

/**
 * Test function for memory usage calculation
 * Verifies that get_bytes() returns correct memory size as blocks are added
 * @tparam Scalar Type of matrix elements
 */
template<typename Scalar>
void test_getbytes() {
  // Create block matrix with dimensions {1,2,3}
  bse::block_matrix<Scalar> m(std::vector<int>{1,2,3});
  using block_type = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;

  // Track total bytes as blocks are added
  std::size_t bytes = 0;
  for (int i=0; i<m.num_blocks(); ++i) {
    for (int j=0; j<m.num_blocks(); ++j) {
      // Verify current byte count
      ASSERT_EQ(bytes, m.get_bytes());
      // Add new block and update byte count
      m.assign(i, j, block_type(m.dim_block(i), m.dim_block(j)));
      bytes += sizeof(Scalar) * m.dim_block(i) * m.dim_block(j);
    }
  }
}
