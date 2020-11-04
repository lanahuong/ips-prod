/**
 * @file OrthogonalityChecker.h
 * Contains the implementation of the OrthogonalityChecker
 */

#ifndef PROJET_IPS1_ORTHOGONALITYCHECKER_H
#define PROJET_IPS1_ORTHOGONALITYCHECKER_H

#include <armadillo>
#include "constants.h"
#include "hermiteCoefs.h"
#include "Hermite.h"
#include "SolverSchrodinger.h"

/**
 * @class OrthogonalityChecker
 * Works with a memory thus needs to be initialised
 * @remark could be implemented more efficiently by dynamically
 * computing the values as needed.
 */
class OrthogonalityChecker {
 public:
  /**
   * Initializes a the functions to compute the Scalar product faster
   */
  explicit OrthogonalityChecker(uint indexMax = HERM_QUADRA_N_MAX, uint nQuadra = HERM_QUADRA_N_MAX);

  /**
   * @brief Check wether two functions are orthogonal
   */
  double checkFor(uint n, uint m);

 private:
  /**
   * @brief To hold the maximum index we're going to test
   */
  uint indexMax;

  /**
   * @brief To hold the index of the choosen quadrature
   */
  uint nQuadra;

  /**
   * @brief Get the weight from the values previously computed in python and included in the header
   * hermite_coefs.h
   */
  arma::Row<double> getWeightVector();

  /**
   * @brief Simple function to compute the pseudo factorials factor given a cache stored in pseudoFactorials
   * so the factorials are only computed once
   */
  void initPseudoFactorial();

  /**
   * @brief Holds the pre computed values
   */
  arma::rowvec pseudoFactorials;

  /**
   * @brief Holds H_n evaluated on a vector composed of the Gauss-Hermite polynomial's roots
   */
  arma::mat hermiteMatrix;

  /**
   * @brief An array of arma::mat. In each matrix the first row corresponds o the
   */
  arma::mat hermiteQuadra[HERM_QUADRA_N_MAX + 1] = HERMITE_QUADRA_COEFS;

  /**
   * @brief Finds and computes the vector of Z used in Hermite::ComputeMatrix
   */
  arma::rowvec getZvector();
};

#endif //PROJET_IPS1_ORTHOGONALITYCHECKER_H
