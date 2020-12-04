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
     * Compute the r part of the function
     * @param rVec
     * @param m Angular momentum ?
     * @param n Principal quantum number
     * @return the r part of the function
     */
    arma::vec rPart(const arma::vec &rVec, int m, int n);
    arma::vec rPart_mem(const arma::vec &rVec, int m, int n);

    /**
     * Compute the z part of the function
     * @param zVec
     * @param nz
     * @return the z part of the function
     */
    arma::vec zPart(const arma::vec &zVec, int nz);
    arma::vec zPart_mem(const arma::vec &zVec, int nz);

    /**
     *
     * @param m Quantum number
     * @param n Quantum number
     * @param nz Quantum number
     * @param rVec Vector of r[i]
     * @param zVec Vector of r[j]
     * @return a matrix where Mat(i,j) corresponds to \f$ \Psi_{m_a, n_a, n_{za}} (r_i, z_j) \f$
     */
    arma::mat basisFunc(int m, int n, int nz, const arma::vec &rVec, const arma::vec &zVec);

    /**
     * Uses memoisation to faster compute the values.
     *
     * @warning The mem depends on the rVec and zVec, if they change the result
     * will be false. To "prevent" it we check the memptr of the arma vecs,
     * but it's only a VERY WEAK PROTECTION. Computing a hash is too slow,
     * the only way would be to keep the vects IN the Basis class. (TODO ?)
     *
     * @param m Quantum number
     * @param n Quantum number
     * @param nz Quantum number
     * @param rVec Vector of r[i]
     * @param zVec Vector of r[j]
     * @return a matrix where Mat(i,j) corresponds to \f$ \Psi_{m_a, n_a, n_{za}} (r_i, z_j) \f$
     */
    arma::mat basisFunc_mem(int m, int n, int nz, const arma::vec &rVec, const arma::vec &zVec);

private:

    double br; /**< Orthogonal deformation parameter */
    double bz; /**< Z deformation parameter */
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


    arma::ivec computed_z_indices;
    std::vector<arma::vec> computed_z_vals;

    arma::imat computed_r_indices;
    std::vector<arma::vec> computed_r_vals;

};