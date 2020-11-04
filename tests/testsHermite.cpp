// calc_test.cpp
#include <gtest/gtest.h>
#include "../src/Hermite.h"
#include "../src/constants.h"

using namespace arma;

TEST(Hermite, ComputeMatrix)
{
  mat mat_zero = mat(1, 1, fill::zeros);
  mat M;
  /** Wrong order
  M = Hermite::ComputeMatrix(-2, 1.0, 2.0, 0.01);
  ASSERT_TRUE(approx_equal(mat_zero, M, "absdiff", EPSILON)) << "Fail with non null matrix while n_max is neg";

  M = Hermite::ComputeMatrix(2, 1.0, -1.0, 0.01);
  ASSERT_TRUE(approx_equal(mat_zero, M, "absdiff", EPSILON)) << "Fail with non null matrix while z_max is smaller than z_min";

  M = Hermite::ComputeMatrix(10, 1.0, 2.0, 0);
  ASSERT_TRUE(approx_equal(mat_zero, M, "absdiff", EPSILON)) << "Fail with non null matrix while step is null";
*/
  mat M_3_0_10_1_computed = Hermite::computeMatrix(3, arma::regspace(0, 1, 10).as_row());
  mat M_3_0_10_1_ref = {{1,  1,  1,  1,   1,   1,   1,    1,    1,    1,    1},
                        {0,  2,  4,  6,   8,   10,  12,   14,   16,   18,   20},
                        {-2, 2,  14, 34,  62,  98,  142,  194,  254,  322,  398},
                        {0,  -4, 40, 180, 464, 940, 1656, 2660, 4000, 5724, 7880}};
  ASSERT_TRUE(approx_equal(M_3_0_10_1_ref, M_3_0_10_1_computed, "absdiff", EPSILON)) << "Fail with wrong matrix params nmax=3, zmin=0, zmax=10, step=1";
}