/**
 * @file NuclearDensityCalculator.h
 */

#ifndef PROJET_IPS1_NUCLEARDENSITYCALCULATOR_H
#define PROJET_IPS1_NUCLEARDENSITYCALCULATOR_H

#include "Basis.h"
#include "constants.h"

/**
 * @class NuclearDensityCalculator
 */
class NuclearDensityCalculator {
private:
    const int N = 14;
    const double Q = 1.3;
    const double br = 1.935801664793151;
    const double bz = 2.829683956491218;
    arma::mat imported_rho_values;
    Basis basis;

    /**
     * Computes the value of rho for the given
     * nuclear parameters and the basis used in the constructor
     * or hardcoded.
     *
     * @warning For now just returns hard coded values because
     * we dont know how to compute them.
     * @return
     */
    inline double rho(int m, int n, int n_z, int mp, int np, int n_zp) const;

public:

    /**
     * Default constructor that uses hard coded valued for rho and basis truncation
     */
    NuclearDensityCalculator();

    /**
     * @see https://dubrayn.github.io/IPS-PROD/project.html#19
     */
    void printRhoDefs();

    /**
     * Cube to hold the correspondance between the 2D matrix of rho values provided by the teacher
     * and the 6D space we're addressing it from
     */
    arma::icube ind;

    /**
     * Naive method seen in the class
     * @param rVals
     * @param zVals
     * @return a matrix
     */
    arma::mat naive_method(const arma::vec& rVals, const arma::vec& zVals);

    /**
     * Optimized method 1
     */
    arma::mat optimized_method1(const arma::vec& rVals, const arma::vec& zVals);

    /**
     * Manually trying to factor out some values, but doesnt really work in a clean manner
     * Plus multithreading is a useless as the first loop has very unequal branches in term
     * of computing needed. Thus the threads are unbalanced
     * @param rVals
     * @param zVals
     * @return
     */
    arma::mat optimized_method2(const arma::vec& rVals, const arma::vec& zVals);

    /**
     * Most optimized method.
     * It uses the memoised version of the @class Basis.
     * It uses multithreading with openMP but can be ported to use native threads
     * It uses the factorisation helper to extract the four sub_sums from the
     * naive one.
     * @param rVals
     * @param zVals
     * @return
     */
    arma::mat optimized_method3(const arma::vec& rVals, const arma::vec& zVals);

    /**
    * @brief Convert the density form cylindric to cartesian coordinates
    * @param xyPoints the number of points on x and y axis
    * @param zPoints the number of points on the z axis
    * @param rVals the r values for which the density was calculated in \a res
    * @param res the matrix of density values in cylindric coordinates
    * @return a cube containing the density in cartesian coordinates
    */
    arma::cube density_cartesian(const int xyPoints, const int zPoints, const arma::vec rVals, const arma::mat res);
};

#endif //PROJET_IPS1_NUCLEARDENSITYCALCULATOR_H
