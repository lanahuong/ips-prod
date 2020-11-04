/**
 * @file Hermite.h
 * This file implements a function to compute all the hermine polynomials on a vector up to a fixed order
 */

#ifndef PROJET_IPS1_HERMITE_H
#define PROJET_IPS1_HERMITE_H

/**
 * Tells the branch predictor that it is very likely that the condition is true
 * Wether it does something or not is implementation dependant.
 */
#define likely(x)       __builtin_expect((x),1)

/**
 * Tells the branch predictor that the condition is very unlikely false.
 */
#define unlikely(x)     __builtin_expect((x),0)

#include <armadillo>
#include <vector>

/**
 * @class Hermite
 * Used to hold a static function to evaluate Hermite polynomials
 */
class Hermite {
 public:
  /**
   * @brief Iteratively evaluate the Hermite polynomial on a vector
   * @param nMax max degree to compute
   * @param zRowvec input vector
   * @return an matrix where each row corresponds to a degree
   */
  static arma::mat computeMatrix(uint nMax, const arma::rowvec &zRowvec);
};

#endif //PROJET_IPS1_HERMITE_H
