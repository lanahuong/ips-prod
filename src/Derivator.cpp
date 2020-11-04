#include "Derivator.h"
#include "Hermite.h"
#include <armadillo>
#include "constants.h"

void Derivator::differentiateTwice(arma::mat &m)
{
  arma::uword n = m.n_cols;
  arma::mat M1 = m;
  arma::mat M2 = m;
  /**
   * If the matrix is too small for the differentiation to make sense,
   * a matrix of zeros of the same size is returned
   */
  if (likely(n > 2))
    {
      M1.shed_col(0);
      M1.insert_cols(n - 1, 1);
      M2.shed_col(n - 1);
      M2.insert_cols(0, 1);

      m = (M1 + M2 - 2 * m) / (STEP * STEP);
    }
  else
    {
      m = arma::zeros(size(m));
    }
}

void Derivator::correctBounds(arma::mat &m)
{
  arma::uword n = m.n_cols;
  if (likely(n > 2))
    {
      m.shed_col(n - 1);
      m.shed_col(0);
    }
}

arma::mat Derivator::differentiate(arma::mat m)
{
  differentiateTwice(m);
  correctBounds(m);
  return m;
}
