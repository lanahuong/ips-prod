#include "../src/NuclearDensityCalculator.h"
#include <gtest/gtest.h>




TEST(NuclearDensityTest, OptiVSRef) {
    NuclearDensityCalculator calc;
    int r_bound = 10;
    int r_step = 1;
    int z_bound = 20;
    int z_step = 2;
    arma::mat ref = calc.naive_method(arma::regspace(-r_bound, r_step, r_bound), arma::regspace(-z_bound, z_step, z_bound));
    arma::mat opti = calc.optimized_method1(arma::regspace(-r_bound, r_step, r_bound), arma::regspace(-z_bound, z_step, z_bound));
    ASSERT_NEAR(arma::norm(opti - ref), 0.0, 1e-08);
}
