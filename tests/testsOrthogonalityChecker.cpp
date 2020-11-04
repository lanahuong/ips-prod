#include <gtest/gtest.h>
#include "../src/OrthogonalityChecker.h"

struct orthoSolutionCheckAllPairs {
    int nMaxIndex;
    int nQuadra;

    /**
     * @brief format when printing a solution_check object
     */
    friend std::ostream &operator<<(std::ostream &os, const orthoSolutionCheckAllPairs &state)
    {
      os << "{nQuadra=" << state.nQuadra << "; nMaxIndex=" << state.nMaxIndex << "}\n";
      return os;
    }
};

/**
 * @class SolutionTest
 * This class defines parameterized tests using instances of solution_check as parameter
 */
class OrthogonalityCheckerTest : public testing::TestWithParam<orthoSolutionCheckAllPairs> {
 public:
  OrthogonalityCheckerTest() = default;;
};




/**
 * @brief Creates an instance of SolutionTest that will run the given test
 */
TEST_P(OrthogonalityCheckerTest, allPairs)
{
  orthoSolutionCheckAllPairs state = GetParam();
  auto *checker = new OrthogonalityChecker(GetParam().nMaxIndex, GetParam().nQuadra);
  for (int i = 0; i <= state.nMaxIndex; i++)
    {
      for (int j = 0; j <= state.nMaxIndex; j++)
        {
          double res = checker->checkFor(i, j);
          if (i + j > (2 * state.nQuadra) - 1)
            {
              ASSERT_TRUE(std::isnan(res));
            }
          else
            {
              if (i != j)
                {
                  ASSERT_TRUE(abs(res) < EPSILON);
                }
              else
                {
                  ASSERT_TRUE(abs(res - 1) < EPSILON);
                }
            }
        }
    }
}


/**
 * Run all SolutionTest tests with given values
 */
INSTANTIATE_TEST_SUITE_P(AllPairs, OrthogonalityCheckerTest, testing::Values(
    orthoSolutionCheckAllPairs{200, HERM_QUADRA_N_MAX},
    orthoSolutionCheckAllPairs{100, HERM_QUADRA_N_MAX},
    orthoSolutionCheckAllPairs{50, HERM_QUADRA_N_MAX / 2},
    orthoSolutionCheckAllPairs{20, HERM_QUADRA_N_MAX}
));
