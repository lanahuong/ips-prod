#include "OrthogonalityChecker.h"

/**
 * The maximum value for n and m is bound by the max degree of the polynomial we have
 * @param maxIndex the maximum index (n or m) we're going to check
 * @param nquadra the index of the quadrature we want to use. It defaults to the maximum
 */
OrthogonalityChecker::OrthogonalityChecker(uint maxIndex, uint nquadra)
{
  this->indexMax = maxIndex;
  this->nQuadra = nquadra < HERM_QUADRA_N_MAX ? nquadra : HERM_QUADRA_N_MAX;
  this->hermiteMatrix = Hermite::computeMatrix(this->indexMax + 1, getZvector());
  this->initPseudoFactorial();
}

/**
 * The Hermite polynomial has a degree of n, thus here its n+m
 * If we force n and m to be smalled than the degree of the quadrature
 *  we can ensure that the result is not an approximation
 * @param n an unsigned int refering to the first function in the scalar product
 * @param m the same for the second
 * @return nan if the values are out of bound
 */
double OrthogonalityChecker::checkFor(uint n, uint m)
{
  if (unlikely(n > indexMax || m > indexMax || n + m > (2 * nQuadra) - 1))
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
  else
    {
      double constFactor = this->pseudoFactorials.at(n) * this->pseudoFactorials.at(m);
      return accu(this->getWeightVector() % this->hermiteMatrix.row(m) % this->hermiteMatrix.row(n)) * constFactor;
    }
}

/**
 * Computes the vector that is the Hermite::ComputeMatrix vector
 * The numerator corresponds to the factor that apprears when we change the
 * integration variable.
 * The z are the roots of the hermite polynomial
 * @return the vector corrected because we change variables in the integral
 */
arma::rowvec OrthogonalityChecker::getZvector()
{
  return (this->hermiteQuadra[this->nQuadra].row(0).as_row()) / (sqrt(MASS * OMEGA / H_BAR));
}

/**
 * Get the weight from the values previously computed in python and included in the header
 * hermite_coefs.h
 * The degree of the quadrature is already stored on the instance.
 * @return a row of doubles
 */
arma::Row<double> inline OrthogonalityChecker::getWeightVector()
{
  return this->hermiteQuadra[this->nQuadra].row(1);
}

/**
 * Corresponds to the factor 1/sqrt(n!2**n), each time we multiply the previous
 * number by 1/(sqrt 2n) so there's no overflow
 */
void OrthogonalityChecker::initPseudoFactorial()
{
  this->pseudoFactorials = arma::rowvec(this->indexMax + 1, arma::fill::ones);
  for (u_long i = 1; i <= this->indexMax; i++)
    {
      this->pseudoFactorials.at(i) = this->pseudoFactorials.at(i - 1) * pow(static_cast<double>(2 * i), -0.5);
    }
  this->pseudoFactorials *= pow((MASS * OMEGA / (PI * H_BAR)), 0.25);
}