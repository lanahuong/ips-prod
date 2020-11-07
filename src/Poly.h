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
  arma::mat hermitePolynomial;
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
     * @param n the rank + 1 of the polynomial
     */
    arma::vec hermite(int n);

    /**
    * @brief Iteratively evaluate the Hermite polynomial on a vector
    * @param nMax max degree to compute
    * @param mMax
    * @param z input vector
    */
    void calcLaguerre(int mMax, int nMax, const arma::vec &z);

    /**
     * @brief Get the Hermite polynomial previously computed of rank n-1
     * @param n the rank + 1 of the polynomial
     * @param m
     */
    arma::vec laguerre(int m, int n);
};

#endif // POLY_H!
