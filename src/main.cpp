#include "main.h"
#include "NuclearDensityCalculator.h"
#include "savers/3Dsavers.h"
#include "Chrono.hpp"

using namespace std;


int main()
{

    NuclearDensityCalculator nuclearDensityCalculator;
    int r_bound = 10;
    double r_step = 01;
    int r_len = 2*r_bound/r_step+1;
    int z_bound = 20;
    double z_step = 02;
    int z_len = 2*z_bound/z_step+1;


    //arma::mat res = nuclearDensityCalculator.naive_method(arma::regspace(-r_bound, r_step, r_bound), arma::regspace(-z_bound, z_step, z_bound));
    Chrono c("main");
    for (int i = 0; i<1000; i++) arma::mat res = nuclearDensityCalculator.optimized_method3(arma::regspace(-r_bound, r_step, r_bound), arma::regspace(-z_bound, z_step, z_bound));

    return 0;
}
