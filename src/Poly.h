/**
 * @file Poly.h
 * This file implements a function to compute all the polynomials on a
 * vector up to a fixed order
 */

#ifndef POLY_H
#define POLY_H

/**
 * Tells the branch predictor that it is very likely that the condition is true
 * Wether it does something or not is implementation dependant.
 */
#define likely(x) __builtin_expect((x), 1)

/**
 * Tells the branch predictor that the condition is very unlikely false.
 */
#define unlikely(x) __builtin_expect((x), 0)

#include <armadillo>
#include <vector>

/**
 * @class Poly
 * Used to compute Hermite and Laguerre polynomials
 */
class Poly {
private:
  /** matrix of hermite polynomials (n,z) */
  arma::mat hermitePolynomial;
  /** cube of laguerre polynomials (m,n,z)*/
  arma::cube laguerrePolynomial;
public:
  /**
   * @brief Iteratively evaluate the Hermite polynomial on a vector
   * @param nMax max degree to compute
   * @param vec input vector
   */
  void calcHermite(uint nMax, const arma::vec &vec);

    /**
     * @brief Get the Hermite polynomial previously computed of rank n-1
     * @param n the rank of the polynomial to get
     */
    arma::vec hermite(int n)const ;

    /**
    * @brief Iteratively evaluate the Laguerre polynomial on a vector
    * @param mMax max m parameter
    * @param nMax max n parameter
    * @param z input vector
    */
    void calcLaguerre(int mMax, int nMax, const arma::vec &z);

    /**
     * @brief Get the Laguerre polynomial previously computed with parameters m and n
     * @param m the m parameter of the polynomial to get
     * @param n the n parameter of the polynomial to get
     */
    arma::vec laguerre(int m, int n)const;
};

#endif // POLY_H!
