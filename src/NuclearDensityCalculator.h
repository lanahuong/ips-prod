/**
 * @file NuclearDensityCalculator.h
 */

#ifndef PROJET_IPS1_NUCLEARDENSITYCALCULATOR_H
#define PROJET_IPS1_NUCLEARDENSITYCALCULATOR_H

#include "Basis.h"

/**
 * @class NuclearDensityCalculator
 */
class NuclearDensityCalculator {
private:
    int N = 14;
    double Q = 1.3;
    double br = 1.935801664793151;
    double bz = 2.829683956491218;
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
    inline double rho(int m, int n, int n_z, int mp, int np, int n_zp);

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
     */
    arma::mat optimized_method2(const arma::vec& rVals, const arma::vec& zVals);

    /**
     * Most optimized method.
     * It uses the memoised version of the @class Basis.
     * It uses multithreading with openMP but can be ported to use native threads
     * It uses the factorisation helper to extract the four sub_sums from the
     * naive one.
     */
    arma::mat optimized_method3(const arma::vec& rVals, const arma::vec& zVals);
};

#endif //PROJET_IPS1_NUCLEARDENSITYCALCULATOR_H
