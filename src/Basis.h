/**
 * @file Basis.h
 * Basis functions
 */

#include <armadillo>

/**
 * @class Basis
 */
class Basis {
public:

    int mMax; /**< Holds the maximum quantum number m for the basis */
    arma::ivec nMax; /**< Holds a vector of principal quantum numbers considered for the basis */
    arma::imat n_zMax; /**< nz quantum numbers */

    /**
     * Basis constructor
     * @param br Basis deformation along radius
     * @param bz Basis deformation along z axis
     * @param N Basis truncation parameter
     * @param Q Basis truncation parameter
     */
    Basis(double br, double bz, int N, double Q);

    /**
     *
     * @param rVec
     * @param m Angular momentum ?
     * @param n Principal quantum number
     * @return
     */
    arma::vec rPart(arma::vec &rVec, int m, int n) const;

    /**
     *
     * @param zVec
     * @param nz
     * @return
     */
    arma::vec zPart(arma::vec &zVec, int nz) const;

    /**
     *
     * @param m
     * @param n
     * @param nz
     * @param rVec
     * @param zVec
     * @return
     */
    arma::mat basisFunc(int m, int n, int nz, arma::vec &rVec, arma::vec zVec);

private:

    long double br; /**< */
    long double bz; /**< */
    arma::cube rPartVals; /**< */
    arma::vec lastRVals; /**< */
    arma::mat zPartVals; /**< */
    arma::vec lastZVals; /**< */

    /**
     * Given the definition and nMax being >=0 , if Q is null the sup is not defined
     * if Q non null, nMax is equal to
     * \f$ \lfloor  (N+2)Q^{-1/3} -0.5 Q^{-1}  \rfloor  \f$
     * @param N
     * @param Q
     */
    static int calcMMax(int N, double Q);

    /**
     * Computes the nMax vector after setNMax was called.
     */
    arma::ivec calcNMax() const;

    /**
     * Initializes the n_zmax matrix
     * @param N
     * @param Q
     */
    arma::imat calcN_zMax(int N, double Q);

};