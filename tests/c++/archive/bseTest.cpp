#include "gtest.h"

#include "../src/bse.hpp"

TEST(BSE, SumFreq){
  using matrix_type = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;
  bse::block_matrix<double> m1( std::vector<int> {1,2,3} );

  matrix_type bm11(1,1);
  matrix_type bm22(2,2);
  matrix_type bm12(1,2);
  matrix_type bm21(2,1);
  matrix_type bm33(3,3);

  bm11.setIdentity();
  bm22.setIdentity();
  bm33.setIdentity();
  bm11 *= 1;
  bm22 *= 2;
  bm33 *= 3;

  for (int i=0; i<2; ++i) {
    bm12(0, i) = 0.1*(i+1);
    bm21(i, 0) = -0.1*(i+1);
  }

  m1.assign(0, 0, bm11);
  m1.assign(1, 1, bm22);
  m1.assign(0, 1, bm12);
  m1.assign(1, 0, bm21);
  m1.assign(2, 2, bm33);

  // bsesol::bse_solver<double>  bse( std::vector<int> {1,2,3});
  bsesol::bse_solver<double>  bse( std::vector<int> {1,2,3}, std::vector<int> {0}, std::vector<int> {0});
  auto m1_sum = bse.Sumfreq_A(m1);
  ASSERT_TRUE(m1_sum(0,0)(0,0)-bm11.sum() < 1e-10);
  ASSERT_TRUE(m1_sum(1,1)(0,0)-bm22.sum() < 1e-10);
  ASSERT_TRUE(m1_sum(0,1)(0,0)-bm12.sum() < 1e-10);
  ASSERT_TRUE(m1_sum(1,0)(0,0)-bm21.sum() < 1e-10);
  ASSERT_TRUE(m1_sum(1,1)(0,0)-bm22.sum() < 1e-10);
}
