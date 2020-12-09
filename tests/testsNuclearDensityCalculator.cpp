/**
 * @file testsNuclearDensityCalculator.cpp
 *
 * This file contains unit test for the class NuclearDensityCalculator
 */

#include <gtest/gtest.h>
#include <armadillo>
#include <vector>

#include "../src/NuclearDensityCalculator.h"


double xyBound = 10;
double zBound = 20;
double xyPoints = 32;
double zPoints = 64;

/**
 * Test suit for all tests of the NuclearDensityCalculator class
*/ 
class NuclearDensityTest : public testing::Test {
    public:
        static NuclearDensityCalculator* ndc;
        static arma::mat* res;
        static arma::vec* rVals;
        static arma::vec* zVals;
    
        static void SetUpTestSuite() {
            ndc = new NuclearDensityCalculator();
            rVals = new arma::vec(arma::linspace(-xyBound, xyBound, xyPoints));
            zVals = new arma::vec(arma::linspace(-zBound, zBound, zPoints));
            res = new arma::mat(ndc->naive_method(*rVals, *zVals));
        }

        static void TearDownTestSuite() {
            delete rVals;
            rVals = nullptr;
            delete zVals;
            zVals = nullptr;
            delete res;
            res = nullptr;
        }
};

NuclearDensityCalculator* NuclearDensityTest::ndc = nullptr;
arma::mat* NuclearDensityTest::res = nullptr;
arma::vec* NuclearDensityTest::rVals = nullptr;
arma::vec* NuclearDensityTest::zVals = nullptr;

/**
 * structure of points to test the density calculation
*/
struct densityPoint {
    double r;
    double z;
    double v;
};

/**
 * @class DensityPointTest
 * This class defines parameterized tests using instances of densityPoint as parameter
 */
class DensityPointTest : public NuclearDensityTest,
                         public testing::WithParamInterface<densityPoint> {
    public:
        arma::cube cube;
        DensityPointTest() {
            cube = ndc->density_cartesian(xyPoints, zPoints, *rVals, *res);
        }
};

struct densityPoint p1 = {2.903226, 5.396825, 0.149483};
struct densityPoint p2 = {-4.838710, -6.666667, 0.042846};
struct densityPoint p3 = {-2.258065, 0.952381 , 0.146820};
struct densityPoint p4 = {6.774194 , 7.301587 , 0.001511};

/**
 * @brief Test the conversion to cartesian coordinates
 */
TEST_P(DensityPointTest, cartesianConversion) {
    struct densityPoint p = GetParam();
    arma::vec r_diff = abs((*rVals - p.r));
    int r = r_diff.index_min();
    arma::vec z_diff = abs((*zVals - p.z));
    int z = z_diff.index_min();

    ASSERT_NEAR(abs(p.v - cube.at(r, xyPoints/2, z)), 0.0, 1e-06);
    ASSERT_NEAR(abs(p.v - cube.at(xyPoints/2, r, z)), 0.0, 1e-06);
}

INSTANTIATE_TEST_SUITE_P(PointsWithYZero, DensityPointTest, testing::Values(
    p1, p2, p3, p4
));
