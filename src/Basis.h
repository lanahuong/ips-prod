/**
 * @file Basis.h
 *
 */

#include <armadillo>

/**
 * @class Basis
 */
class Basis {
public:

    int mMax;
    arma::ivec nMax;
    arma::imat n_zMax;

    Basis(double br, double bz, int N, double Q);

    arma::vec rPart(arma::vec &rVec, int m, int n);

    arma::vec zPart(arma::vec &zVec, int nz);

    arma::mat basisFunc(int m, int n, int nz, arma::vec &rVec, arma::vec zVec);

private:

    double br;
    double bz;
    arma::cube rPartVals;
    arma::vec lastRVals;
    arma::mat zPartVals;
    arma::vec lastZVals;

    /**
     * Given the definition and nMax being >=0 , if Q is null the sup is not defined
     * if Q non null, nMax is equal to
     * \f$ \lfloor  (N+2)Q^{-1/3} -0.5 Q^{-1}  \rfloor  \f$
     * @param N
     * @param Q
     */
    void findMMax(int N, double Q);

};