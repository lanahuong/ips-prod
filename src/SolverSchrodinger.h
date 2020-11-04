/**
 * @file SolverSchrodinger.h
 *
 * This file contains the SolverSchrodinger class that computes the final
 * solution of the Schrodinger equation.
 */
#ifndef SOLVERSCHRODINGER_H
#define SOLVERSCHRODINGER_H

#include <armadillo>
#include <vector>
#include <cmath>
#include "Hermite.h"
#include "Derivator.h"

/**
 * @class SolverSchrodinger
 *
 * This class solves the Schrodinger equation
 */
class SolverSchrodinger {
 public:
  /**
   * @brief Solve the equation between \a zmin and \a zmax with the step \a
   * step for energie levels from 0 to \a n in 1 dimension
   */
  static arma::mat solve1D(double zmin, double zmax, uint n);

  /**
   * @brief Solve the equation between \a zmin and \a zmax with the step \a
   * step for energie levels from 0 to \a n in 1 dimension
   */
  static arma::mat solve1D(const arma::rowvec &z, uint n);

  /**
   * @brief Test the solution
   */
  static bool test1DSolution(const arma::rowvec &z, arma::mat phi);
};

#endif // !SOLVERSCHRODINGER_H
