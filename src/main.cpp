#include "main.h"
#include "NuclearDensityCalculator.h"
#include "Saver.h"

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

    Saver::saveToCSV(res, "tmp/density-r-z.csv");
    
    arma::cube cube = NuclearDensityCalculator::density_cartesian(xyPoints, zPoints, rVals, res);

    Saver::cubeToDf3(cube, "tmp/density-r-z.df3");

    return 0;
}
