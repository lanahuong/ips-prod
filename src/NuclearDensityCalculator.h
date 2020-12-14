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
    const int N = 14; /** truncation parameter */
    const double Q = 1.3; /** truncation parameter */
    const double br = 1.935801664793151; /** radius deformation factor */
    const double bz = 2.829683956491218; /** z deformation factor */
    arma::mat imported_rho_values; /** rho values from file */
    Basis basis; /** basis of functions */

    /**
     * Computes the value of rho for the given
     * nuclear parameters and the basis used in the constructor
     * or hardcoded.
     * @param n quantum number of a
     * @param m quantum number of a
     * @param n_z quantum number of a
     * @param np quantum number of b
     * @param mp quantum number of b
     * @param n_zp quantum number of b
     *
     * @warning For now just returns hard coded values because
     * we dont know how to compute them.
     * @return rho(n, m, n_z, np, mp, n_zp)
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
     * @param rVals vector of r values (radius)
     * @param zVals vector of z values
     * @return a matrix of density values for rVals x zVals (cartesian products giving coordinates)
     */
    arma::mat naive_method(const arma::vec& rVals, const arma::vec& zVals) ;

    /**
     * Optimized method 1
     * @param rVals vector of r values (radius)
     * @param zVals vector of z values
     * @return a matrix of density values for rVals x zVals (cartesian products giving coordinates)
     */
    arma::mat optimized_method1(const arma::vec& rVals, const arma::vec& zVals) ;

    /**
     * Manually trying to factor out some values, but doesnt really work in a clean manner
     * Plus multithreading is a useless as the first loop has very unequal branches in term
     * of computing needed. Thus the threads are unbalanced
     * @param rVals vector of r values (radius)
     * @param zVals vector of z values
     * @return a matrix of density values for rVals x zVals (cartesian products giving coordinates)
     */
    arma::mat optimized_method2(const arma::vec& rVals, const arma::vec& zVals);

    /**
     * Most optimized method.
     * It uses the memoised version of the class Basis.
     * It uses multithreading with openMP but can be ported to use native threads
     * It uses the factorisation helper to extract the four sub_sums from the
     * naive one.
     * @param rVals vector of r values (radius)
     * @param zVals vector of z values
     * @return a matrix of density values for rVals x zVals (cartesian products giving coordinates)
     */
    arma::mat optimized_method3(const arma::vec& rVals, const arma::vec& zVals) const;

    /**
    * @brief Convert the density form cylindric to cartesian coordinates
    * @param xyPoints the number of points on x and y axis
    * @param zPoints the number of points on the z axis
    * @param rVals the r values for which the density was calculated in \a res
    * @param res the matrix of density values in cylindric coordinates
    * @return a cube containing the density in cartesian coordinates
    */
    arma::cube static density_cartesian(int xyPoints, int zPoints, const arma::vec& rVals, const arma::mat& res) ;
};

#endif //PROJET_IPS1_NUCLEARDENSITYCALCULATOR_H
