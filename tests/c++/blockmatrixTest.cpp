#include "blockmatrixTest.hpp"

// Test constructor of BlockMatrix class
TEST(BlockMatrix, Constructor){
  // Initialize dimensions vector
  std::vector<int> dims = {1, 2, 3};
  // Create block matrix with given dimensions
  bse::block_matrix<double> m(dims);

  // Verify all blocks are initially zero/empty
  for (int b=0; b < m.num_blocks(); ++b) {
    for (int b2=0; b2 < m.num_blocks(); ++b2) {
      ASSERT_TRUE(!m.exists(b, b2));
    }
  }

  // Check total dimension equals sum of block dimensions
  assert(m.dim() == std::accumulate(dims.begin(), dims.end(), 0));
}

// Test assignment operations for different data types
TEST(BlockMatrix, Assignment){
  test_assignment<double>();
  test_assignment<std::complex<double>>();
}

// Test memory usage calculation
TEST(BlockMatrix, GetBytes){
  test_getbytes<double>();
  test_getbytes<std::vector<double>>();
}

// Test scalar multiplication of block matrix
TEST(BlockMatrix, BlockMatrix_Multiply_Scalar_Test) {
  // Create block matrix with dimensions {1,2,3}
  bse::block_matrix<double> m(std::vector<int>{1,2,3});
  // Create 3x3 Eigen matrix
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> bm(3,3);

  // Fill matrix with test values
  for (int j=0; j<3; ++j) {
    for (int i=0; i<3; ++i) {
      bm(i,j) = i+10*j;
    }
  }
  // Assign block and multiply by scalar
  m.assign(2, 2, bm);
  m = 10.0*m;

  // Get view of block and verify multiplication
  auto bm_view = m(2,2);
  for (int j=0; j<2; ++j) {
    for (int i = 0; i < 2; ++i) {
      ASSERT_TRUE(bm_view(i,j) == 10.0*bm(i,j));
    }
  }
}

// Test addition of block matrices
TEST(BlockMatrix, BlockMatrix_Add_Test) {
  // Create two block matrices
  bse::block_matrix<double> m(std::vector<int>{1,2,3});
  bse::block_matrix<double> m1(std::vector<int>{1,2,3});

  // Create test matrices
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> bm(3,3);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> bm13(1,3);

  // Fill 3x3 matrix
  for (int j=0; j<3; ++j) {
    for (int i=0; i<3; ++i) {
      bm(i,j) = i+10*j;
    }
  }

  // Fill 1x3 matrix
  for (int i=0; i<3; ++i) {
    bm13(0, i) = i;
  }

  // Assign blocks and perform addition
  m.assign(2, 2, bm);
  m1.assign(2, 2, bm);
  m1.assign(0, 2, bm13);

  m = m+m1;
  m1 = 2.0*m1;
  
  auto bm_view = m(2,2);
  auto bm1_view = m1(2,2);
  
  // Test addition of non-zero blocks
  for (int j=0; j<3; ++j) {
    for (int i = 0; i < 3; ++i) {
      ASSERT_TRUE(bm_view(i,j) == bm1_view(i,j));
    }
  }

  // Test addition with zero block
  auto bm2_view=m1(0,2);
  for (int i=0; i<3; ++i) {
    ASSERT_TRUE(2.0*bm13(0, i) == bm2_view(0, i) );
  }
}

// Test subtraction of block matrices
TEST(BlockMatrix, BlockMatrix_Minus_Test) {
  // Create two block matrices
  bse::block_matrix<double> m(std::vector<int>{1,2,3});
  bse::block_matrix<double> m1(std::vector<int>{1,2,3});

  // Create test matrix
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> bm(3,3);

  // Fill matrix with test values
  for (int j=0; j<3; ++j) {
    for (int i=0; i<3; ++i) {
      bm(i,j) = i+10*j;
    }
  }

  // Assign blocks and subtract
  m.assign(2, 2, bm);
  m1.assign(2, 2, bm);

  m = m-m1;
  
  // Verify result is zero
  auto bm_view = m(2,2);
  for (int j=0; j<3; ++j) {
    for (int i = 0; i < 3; ++i) {
      ASSERT_TRUE(fabs(bm_view(i,j))<1e-10);
    }
  }
}

// Test multiplication of block matrices
TEST(BlockMatrix, BlockMatrix_Multiply_Test) {
  // Create block matrices
  bse::block_matrix<double> m(std::vector<int>{2,2});
  bse::block_matrix<double> m1(std::vector<int>{2,2});
  
  // Create test matrices
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> bm(2,2);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> bm1(2,2);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Largebm(4,4);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Largebm1(4,4);
  Largebm.setZero();
  Largebm1.setZero();

  // Fill 00 block
  for (int j=0; j<2; ++j) {
    for (int i=0; i<2; ++i) {
      bm(i,j) = i+10*j;
      bm1(i, j) = 10*i+j;
      Largebm(i, j)= i+10*j;
      Largebm1(i, j)= 10*i+j;
    }
  }
  m.assign(0,0, bm);
  m1.assign(0, 0, bm1);

  // Fill 01 block
  for (int j=0; j<2; ++j) {
    for (int i=0; i<2; ++i) {
      bm(i,j) = i-10*j;
      bm1(i, j) = 10*i-j;
      Largebm(i, j+2)= i-10*j;
      Largebm1(i, j+2)= 10*i-j;
    }
  }
  m.assign(0,1, bm);
  m1.assign(0,1, bm1);

  // Fill 10 block
  for (int j=0; j<2; ++j) {
    for (int i=0; i<2; ++i) {
      bm(i,j) = 2*i-10*j;
      bm1(i, j) = 10*i-2*j;
      Largebm(i+2, j)= 2*i-10*j;
      Largebm1(i+2, j)= 10*i-2*j;
    }
  }
  m.assign(1,0, bm);
  m1.assign(1,0, bm1);

  // Multiply matrices
  m=m*m1;
  Largebm = Largebm*Largebm1;
  
  // Verify results
  for (int j=0; j<4; ++j) {
    for (int i = 0; i<4; ++i) {
      auto bm_view=m(int(i/2), int(j/2));
      ASSERT_TRUE(bm_view(i%2, j%2)==Largebm(i, j));
    }
  }
}

// Test block gathering functionality
TEST(BlockMatrix, Gather_Blocks_Test) {
  int super_block_size = 5;
  int num_blocks = 3 * super_block_size;
  // Define connection rule
  auto connected = [](int i, int j) {return i%3== j%3;};

  // Gather connected blocks
  auto bl_groups = bse::detail::gather_blocks(connected, num_blocks);

  // Verify results
  ASSERT_TRUE(bl_groups.size() == 3);
  for (auto gr : bl_groups) {
    ASSERT_TRUE(gr.size() == super_block_size);
  }
}

// Test matrix inversion
TEST(BlockMatrix, Inverse_Test) {
  using matrix_type = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;
  bse::block_matrix<double> m1(std::vector<int>{1,2,3});

  // Create test matrices
  matrix_type bm33(3,3);
  matrix_type bm11(1,1);
  matrix_type bm22(2,2);
  matrix_type bm12(1,2);

  // Initialize identity matrices with different scales
  bm11.setIdentity();
  bm22.setIdentity();
  bm33.setIdentity();
  bm11 *= 1;
  bm22 *= 2;
  bm33 *= 3;

  // Fill 1x2 block
  for (int i=0; i<2; ++i) {
    bm12(0, i) = 0.1*i;
  }

  // Assign blocks
  m1.assign(0, 0, bm11);
  m1.assign(1, 1, bm22);
  m1.assign(2, 2, bm33);
  m1.assign(0, 1, bm12);

  // Calculate inverse
  auto inv_m1 = bse::inverse(m1);

  // Verify dimensions and block structure
  ASSERT_TRUE(inv_m1.dim() == m1.dim());
  for (int bl=0; bl < m1.num_blocks(); ++bl) {
    for (int bl2=0; bl2 < m1.num_blocks(); ++bl2) {
      ASSERT_TRUE(
          (m1.exists(bl,bl2) || m1.exists(bl2,bl)) == (inv_m1.exists(bl,bl2) || inv_m1.exists(bl2,bl))
      );
    }
  }
  
  // Check accuracy of inverse
  matrix_type diff = inv_m1.convert_to_matrix() - m1.convert_to_matrix().inverse();
  ASSERT_TRUE(diff.cwiseAbs().maxCoeff() < 1e-10);

  // Test inverse difference
  bm11.setIdentity();
  bm22.setIdentity();
  bm33.setIdentity();
  bm11 *= 3;
  bm22 *= 1;
  bm33 *= 2;

  for (int i=0; i<2; ++i) {
    bm12(0, i) = 0.3*i;
  }
  
  // Create second test matrix
  bse::block_matrix<double> m2(std::vector<int>{1,2,3});
  m2.assign(0, 0, bm11);
  m2.assign(1, 1, bm22);
  m2.assign(2, 2, bm33);
  m2.assign(0, 1, bm12);
  
  // Test difference of inverses
  auto inv_m1_minus_m2 = bse::inverse(m1)-bse::inverse(m2);
  matrix_type diff2 = inv_m1_minus_m2.convert_to_matrix() - (m1.convert_to_matrix().inverse()- m2.convert_to_matrix().inverse());
  ASSERT_TRUE(diff2.cwiseAbs().maxCoeff() < 1e-10);
}

/*
// Commented out test for diagonal extraction
TEST(BlockMatrix, asDiagonal_TEST) {
  using matrix_type = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;
  bse::block_matrix<double> m1(std::vector<int>{1,2,3});

  matrix_type bm33(3,3);
  matrix_type bm11(1,1);
  matrix_type bm22(2,2);
  matrix_type bm12(1,2);
  matrix_type bm21(2,1);

  bm11.setIdentity();
  bm22.setIdentity();
  bm33.setIdentity();
  bm11 *= 1;
  bm22 *= 2;
  bm33 *= 3;

  for (int i=0; i<2; ++i) {
    bm12(0, i) = 0.1*i;
    bm21(i, 0) = 0.1*i;
  }

  m1.assign(0, 0, bm11);
  m1.assign(1, 1, bm22);
  m1.assign(2, 2, bm33);
  m1.assign(0, 1, bm12);
  m1.assign(1, 0, bm21);

  auto diagonal_m1 = bse::asDiagonal(m1);

  ASSERT_TRUE(diagonal_m1.dim() == m1.dim());
  for (int bl=0; bl < m1.num_blocks(); ++bl) {
    for (int bl2=0; bl2 < m1.num_blocks(); ++bl2) {
      ASSERT_TRUE(
              (m1.exists(bl,bl2) || m1.exists(bl2,bl)) == (diagonal_m1.exists(bl,bl2) || diagonal_m1.exists(bl2,bl))
      );
    }
  }

  std::cout << diagonal_m1.convert_to_matrix() << std::endl;
}
*/