/**
 * @file tests_derivator.cpp
 *
 * This file contains unit test for the class Derivator
 */

#include <gtest/gtest.h>

#include "../src/Derivator.h"
#include "../src/constants.h"

using namespace arma;

mat test1()
{
  mat M = zeros(5, 6);
  for (int i = 1; i <= 5; ++i)
    {
      for (int j = 1; j <= 6; ++j)
        {
          M(i - 1, j - 1) = i + j - 1;
        }
    }
  return M;
}

mat test2()
{
  mat M = zeros(8, 7);
  for (int i = 1; i <= 8; ++i)
    {
      for (int j = 1; j <= 7; ++j)
        {
          M(i - 1, j - 1) = pow((i * STEP + j * STEP + 0.3), 2);
        }
    }
  return M;
}

TEST(Derivator, differentiate)
{
  mat m0 = eye(3, 2);
  mat m1 = eye(1, 1);
  mat m2 = test1();
  mat m3 = test2();
  mat m4 = {{1, 2,    7,   2.5, 15,   1.45, 8},
            {3, 8,    1.4, -12, 2.87, 4,    9},
            {6, -1.2, 5,   7,   -22,  5.3,  1}};

  mat m4_res = {{4,     -9.5, 17,    -26.05, 20.1},
                {-11.6, -6.8, 28.27, -13.74, 3.87},
                {13.4,  -4.2, -31,   56.3,   -31.6}};

  m4_res = m4_res / (STEP * STEP);

  ASSERT_TRUE(approx_equal(Derivator::differentiate(m0), zeros(3, 2), "absdiff", EPSILON));
  ASSERT_TRUE(approx_equal(Derivator::differentiate(m1), zeros(1, 1), "absdiff", EPSILON));
  ASSERT_TRUE(approx_equal(Derivator::differentiate(m2), zeros(5, 4), "absdiff", EPSILON));
  ASSERT_TRUE(approx_equal(Derivator::differentiate(m3), 2 * ones(8, 5), "absdiff", EPSILON));
  ASSERT_TRUE(approx_equal(Derivator::differentiate(m4), m4_res, "absdiff", EPSILON));
}





