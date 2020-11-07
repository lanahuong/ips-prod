/**
 * @file Basis.h
 * Basis functions
 */

#include <armadillo>
#include "Poly.h"

/**
 * @class Basis
 */
class Basis {
public:

    int mMax{}; /**< Holds the maximum quantum number m for the basis */
    arma::ivec nMax; /**< Holds a vector of principal quantum numbers considered for the basis */
    arma::imat n_zMax; /**< nz quantum numbers */

    Basis();

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
    arma::vec rPart(const arma::vec &rVec, int m, int n) const;

    /**
     *
     * @param zVec
     * @param nz
     * @return
     */
    arma::vec zPart(const arma::vec &zVec, int nz) const;

    /**
     *
     * @param m Quantum number
     * @param n Quantum number
     * @param nz Quantum number
     * @param rVec Vector of r[i]
     * @param zVec Vector of r[j]
     * @return a matrix where Mat(i,j) corresponds to \f$ \Psi_{m_a, n_a, n_{za}} (r_i, z_j) \f$
     */
    arma::mat basisFunc(int m, int n, int nz, const arma::vec &rVec, const arma::vec &zVec) const;

private:

    long double br{}; /**< Orthogonal deformation parameter */
    long double bz{}; /**< Z deformation parameter */
    arma::cube rPartVals; /**< */
    arma::vec lastRVals; /**< */
    arma::mat zPartVals; /**< */
    arma::vec lastZVals; /**< */
    Poly poly;

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