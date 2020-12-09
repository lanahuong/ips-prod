#include "main.h"
#include "NuclearDensityCalculator.h"
#include "savers/3Dsavers.h"

using namespace std;

int main()
{
     NuclearDensityCalculator nuclearDensityCalculator;

    double xyBound = 10;
    double zBound = 20;
    double xyPoints = 32;
    double zPoints = 64;
     
    arma::mat rVals = arma::linspace(-xyBound, xyBound, xyPoints);
    arma::mat zVals = arma::linspace(-zBound, zBound, zPoints);
    arma::mat res = nuclearDensityCalculator.optimized_method3(rVals, zVals);
    
    arma::cube cube = nuclearDensityCalculator.density_cartesian(xyPoints, zPoints, rVals, res);

    std::cout << cubeToDf3(cube);

    return 0;
}
