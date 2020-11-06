#include "Poly.h"

void Poly::calcHermite(uint nMax, const arma::vec &zRowvec) {
  /**
   * If the parameters are nonsense the matrix (0) is returned
   */
  arma::uword rowLen = zRowvec.size();
  hermitePolynomial = arma::mat(nMax + 1, rowLen, arma::fill::zeros);
  hermitePolynomial.row(0) = arma::vec(rowLen, arma::fill::ones).t();
  if (likely(nMax > 0)) {
    hermitePolynomial.row(1) = 2 * zRowvec.t();
    for (arma::uword i = 2; i <= nMax; i++) {
      hermitePolynomial.row(i) =
          hermitePolynomial.row(1) % hermitePolynomial.row(i - 1) -
          2. * ((double)i - 1.) * hermitePolynomial.row(i - 2);
    }
  }
}

arma::vec Poly::hermite(int n) { return hermitePolynomial.row(n).t(); }
