#include "SolverSchrodinger.h"
#include "Derivator.h"
#include "constants.h"

/**
 * @param zmin minimum z where to compute the function
 * @param zmax maximum z where to compute the function
 * @param n maximum energy level
 * @return a matrix of solutions where enery vary with the row and z with the column
 */
arma::mat SolverSchrodinger::solve1D(double zmin, double zmax, uint n)
{
  // Compute the factor of the solution with n constant
  arma::rowvec z = arma::regspace(zmin, STEP, zmax).as_row();
  return SolverSchrodinger::solve1D(z, n);
}

/**
 * @param z the vector of z values where function were evaluated
 * @param n maximum energy level
 * @return a matrix of solutions where enery vary with the row and z with the column
 */
arma::mat SolverSchrodinger::solve1D(const arma::rowvec &z, uint n)
{
  // Compute the factor of the solution with n constant
  arma::vec nwiseconst = arma::regspace(0, n);
  nwiseconst[0] = 1;
  nwiseconst = arma::pow(arma::sqrt(nwiseconst), -1);

  for (uint i = 1; i <= n; i++)
    {
      nwiseconst[i] = nwiseconst[i] * nwiseconst[i - 1] / 2;
    }

  nwiseconst = nwiseconst * sqrt(sqrt(MASS * OMEGA / (arma::datum::pi * H_BAR)));

  // Compute the factor of the solution with n constant
  arma::rowvec zwiseconst = arma::square(z);
  zwiseconst = zwiseconst * -1 * MASS * OMEGA / (2 * H_BAR);
  zwiseconst = arma::exp(zwiseconst);

  // Compute the final solution
  arma::mat result = nwiseconst * zwiseconst;
  arma::mat hermitePolynomes = Hermite::computeMatrix(n, z);
  result = result % hermitePolynomes;
  return result;
}

/**
 *
 * @param z the vector of z values where function were evaluated
 * @param phi a matrix of functions
 * @return true only if the functions in the matrix \a phi verify the equation
 * @attention need a matrix of at least 3 columns (for the differentiation)
 */
bool SolverSchrodinger::test1DSolution(const arma::rowvec &z, arma::mat phi)
{
  // Compute second derivative of each function
  arma::mat dzsecond = Derivator::differentiate(phi);

  // z² in a matrix
  arma::mat ztrunc = arma::vec(phi.n_rows, arma::fill::ones) * arma::square(z.cols(1, z.n_cols - 2));

  // Truncate phi
  arma::mat phitrunc = phi.cols(1, phi.n_cols - 2);

  // Compute left member
  arma::mat left = -H_BAR * H_BAR * dzsecond / (2. * MASS);
  left = left + (MASS * OMEGA * OMEGA / 2.) * (ztrunc % phitrunc);

  // Compute right member
  arma::mat E = (arma::regspace(0, (double) phi.n_rows - 1) + 1. / 2.) * arma::rowvec(phitrunc.n_cols, arma::fill::ones) * H_BAR * OMEGA;
  arma::mat right = E % phitrunc;

  return approx_equal(left, right, "absdiff", EPSILON);
}
