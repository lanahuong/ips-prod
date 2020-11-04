/**
 * @file tests_solverschrodinger.cpp
 *
 * This file contains unit test for the class SolverSchrodinger
 */

#include <gtest/gtest.h>
#include <armadillo>
#include <vector>

#include "../src/SolverSchrodinger.h"
#include "../src/constants.h"

/**
 * @struct define the parameter of SolutionBoundedTest
 */
struct solution_bounded_state {
    double zmin;
    double zmax;
    uint n;
    // expected number of rows in the resulting matrix
    uint length;

    /**
     * @brief format when printing a solution_bounded_state object
     */
    friend std::ostream &operator<<(std::ostream &os, const solution_bounded_state &state)
    {
      os << "zmin=" << state.zmin << "; zmax=" << state.zmax;
      os << "; n=" << state.n << "; length=" << state.length << "}\n";
      return os;
    }
};

/**
 * @class SolutionBoundedTest
 * This class defines parameterized tests using instances of solution_bounded_state as parameter
 */
class SolutionBoundedTest : public testing::TestWithParam<solution_bounded_state> {
 public:
  SolutionBoundedTest()
  = default;;
};

/**
 * @brief Creates an instance of SolutionBoundedTest that will run the given test
 */
TEST_P(SolutionBoundedTest, solutionSize)
{
  solution_bounded_state state = GetParam();
  arma::mat result = SolverSchrodinger::solve1D(state.zmin, state.zmax, state.n);

  EXPECT_EQ(state.n + 1, result.n_rows);
  EXPECT_EQ(state.length, result.n_cols);
}

/**
 * Run all SolutionBoundedTest tests with given values
 */
INSTANTIATE_TEST_SUITE_P(BasicCases, SolutionBoundedTest, testing::Values(
    solution_bounded_state{-1, 1, 1, uint((2 / STEP) + 1)},
    solution_bounded_state{-10, 10, 1, uint((20 / STEP) + 1)},
    solution_bounded_state{-1, 1, 4, uint((2 / STEP) + 1)}
));

/**
 * @struct solution_from_vec_state define the parameter of SolutionFromVecTest
 */
struct solution_from_vec_state {
    arma::rowvec z;
    uint n;
    // expected number of rows in the resulting matrix
    uint length;

    /**
     * @brief format when printing a solution_from_vec_state object
     */
    friend std::ostream &operator<<(std::ostream &os, const solution_from_vec_state &state)
    {
      os << "z=" << state.z;
      os << "; n=" << state.n << "; length=" << state.length << "}\n";
      return os;
    }
};

/**
 * @class SolutionFromVecTest
 * This class defines parameterized tests using instances of solution_from_vec_state as parameter
 */
class SolutionFromVecTest : public testing::TestWithParam<solution_from_vec_state> {
 public:
  SolutionFromVecTest()
  = default;;
};

/**
 * @brief Creates an instance of SolutionFromVecTest that will run the given test
 */
TEST_P(SolutionFromVecTest, solutionSize)
{
  solution_from_vec_state state = GetParam();
  arma::mat result = SolverSchrodinger::solve1D(state.z, state.n);

  EXPECT_EQ(state.n + 1, result.n_rows);
  EXPECT_EQ(state.length, result.n_cols);
}

arma::rowvec z1 = arma::rowvec({-1, 0, 1});
arma::rowvec z2 = arma::rowvec({-1, -0.5, 0, 0.5, 1});

/**
 * Run all SolutionFromVecTest tests with given values
 */
INSTANTIATE_TEST_SUITE_P(BasicCases, SolutionFromVecTest, testing::Values(
    solution_from_vec_state{z1, 1, 3},
    solution_from_vec_state{z2, 1, 5},
    solution_from_vec_state{z1, 4, 3}
));

/**
 * @struct check_solution_state define the parameter of CheckSolutionTest
 */
struct check_solution_state {
    arma::rowvec z;
    arma::mat solution;
    // expected check result
    bool isSolution;

    /**
     * @brief format when printing a check_solution_state object
     */
    friend std::ostream &operator<<(std::ostream &os, const check_solution_state &state)
    {
      os << "z=" << state.z << "; solution=";
      os << "; isSolution=" << state.isSolution << "\n";
      return os;
    }
};

/**
 * @class CheckSolutionTest
 * This class defines parameterized tests using instances of check_solution_state as parameter
 */
class CheckSolutionTest : public testing::TestWithParam<check_solution_state> {
 public:
  CheckSolutionTest()
  = default;;
};

/**
 * @brief Creates an instance of CheckSolutionTest that will run the given test
 */
TEST_P(CheckSolutionTest, checkSolution)
{
  check_solution_state state = GetParam();
  bool result = SolverSchrodinger::test1DSolution(state.z, state.solution);

  EXPECT_EQ(state.isSolution, result);
}

/**
 * Run all CheckSolutionTest tests with given values
 */
INSTANTIATE_TEST_SUITE_P(BasicCases, CheckSolutionTest, testing::Values(
    check_solution_state{z2, z2, false}
));

