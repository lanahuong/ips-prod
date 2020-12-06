#ifndef BASIS_H
#define BASIS_H

/**
 * @file Basis.h
 * Basis functions
 */

#include <armadillo>
#include "Poly.h"

/**
 * @class Basis
 *
 * @warning methods and variables suffied with _mem require the use of the constructor
 * which receives the values of zVals and rVals
 */
class Basis {
public:

    const int mMax{}; /**< Holds the maximum quantum number m for the basis */
    const arma::ivec nMax; /**< Holds a vector of principal quantum numbers considered for the basis */
    const arma::imat n_zMax; /**< nz quantum numbers */


    Basis() = default;
    /**
     * Basis constructor
     * @param BR Basis deformation along radius
     * @param bz Basis deformation along z axis
     * @param N Basis truncation parameter
     * @param Q Basis truncation parameter
     */
    Basis(double BR, double bz, int N, double Q);

    /**
     * Special constructor to use if we want to accelerate computing by
     * pre-computing some values and storing the computed ones.
     * @param rVals
     * @param zVals
     */
    Basis(double BR, double bz, int N, double Q, const arma::vec& rVals, const arma::vec& zVals);

    /**
     * Compute the r part of the function
     * @param rVec
     * @param m Angular momentum ?
     * @param n Principal quantum number
     * @param use_mem Is set to true if the function is called from a memoized context
     * Protects the user if the user calls the function directly while the class was
     * set to be memoised.
     * @return the r part of the function
     */
    arma::vec rPart(const arma::vec& rVec, int m, int n, bool use_mem = false);

    /**
     * Function used to access memoised values.
     */
    arma::vec rPart_mem(int m, int n);

    /**
     * Compute the z part of the function
     * @param zVec
     * @param nz
     * @param use_mem Is set to true if the function is called from a memoized context
     * Protects the user if the user calls the function directly while the class was
     * set to be memoised.
     * @return the z part of the function
     */
    arma::vec zPart(const arma::vec& zVec, int nz, bool use_mem = false);

    /**
     * Function used to access memoised values.
     */
    arma::vec zPart_mem(int nz);

    /**
     *
     * @param m Quantum number
     * @param n Quantum number
     * @param nz Quantum number
     * @param rVec Vector of r[i]
     * @param zVec Vector of r[j]
     * @return a matrix where Mat(i,j) corresponds to \f$ \Psi_{m_a, n_a, n_{za}} (r_i, z_j) \f$
     */
    arma::mat basisFunc(int m, int n, int nz, const arma::vec& rVec, const arma::vec& zVec);
    arma::mat basisFunc_mem(int m, int n, int nz);

private:
    double br{}; /**< Orthogonal deformation parameter */
    double bz{}; /**< Z deformation parameter */
    Poly poly{}; /**< Polynomial class evaluator */


    bool is_mem = false; /**< Is set to true if the rVec and zVec were given at construction -> memoisation  */
    arma::vec rvec_mem;/**< rVec given in the constructor */
    arma::vec zvec_mem;/**< zVec given in the constructor */
    arma::vec squared_rarg_mem{};/**< pre-computed vector if is_mem */
    arma::vec squared_zarg_mem{};/**< pre-computed vector if is_mem */
    arma::vec zexp_mem{};/**< pre-computed vector if is_mem */
    arma::vec rexp_mem{};/**< pre-computed vector if is_mem */
    std::vector<bool> computed_z_indices;/**< computed indices if is_mem */
    std::vector<arma::vec> computed_z_vals;/**< stored zVals if is_mem */
    std::vector<bool> computed_r_indices;/**< computed indices if is_mem */
    std::vector<arma::vec> computed_r_vals;/**< stored rVals if is_mem */

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

#endif