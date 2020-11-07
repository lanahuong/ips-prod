#include "Poly.h"

void Poly::calcHermite(uint nMax, const arma::vec &vec) {
    /**
     * If the parameters are nonsense the matrix (0) is returned
     */
    uint rowLen = vec.size();
    hermitePolynomial = arma::mat(nMax + 1, rowLen, arma::fill::zeros);
    hermitePolynomial.row(0) = arma::vec(rowLen, arma::fill::ones).t();
    if (nMax > 0) {
        hermitePolynomial.row(1) = 2 * vec.t();
        for (uint i = 2; i <= nMax; i++) {
            hermitePolynomial.row(i) =
                    hermitePolynomial.row(1) % hermitePolynomial.row(i - 1) -
                    2. * ((double) i - 1.) * hermitePolynomial.row(i - 2);
        }
    }
}

arma::vec Poly::hermite(int n) { return hermitePolynomial.row(n).as_col(); }

arma::vec Poly::laguerre(int m, int n) {
    return laguerrePolynomial.slice(n).row(m).as_col();
}

void Poly::calcLaguerre(int mMax, int nMax, const arma::vec &z) {
    /**
     * Init of the cube and helper vectors. Initializes also the first two slices (if needed) for the recursion
     */
    laguerrePolynomial = arma::cube(mMax, z.n_elem, nMax, arma::fill::ones);
    arma::vec m = arma::regspace(0, mMax - 1).as_col();
    arma::rowvec row_ones = arma::vec(z.n_elem, arma::fill::ones).as_row();
    arma::colvec col_ones = arma::vec(mMax, arma::fill::ones).as_col();
    if (nMax > 1) laguerrePolynomial.slice(1) = 1 + (m * row_ones) - col_ones * z.as_row();

    arma::mat coef1;
    arma::mat coef2;
    /**
     * Now we build our cube, each slice is a matrix and within such slice N is constant, so we
     * will operate on slices, its faster.
     * The slice 0 is already filled with ones
     */
    for (int depth = 2; depth < nMax; depth++) {
        coef1 = 2 + (m * row_ones - col_ones * z.as_row() - 1) / (double) depth;
        coef2 = 1 + (m - 1) / (double) depth * row_ones;
        laguerrePolynomial.slice(depth) = coef1 % laguerrePolynomial.slice(depth - 1) - coef2 % laguerrePolynomial.slice(depth - 2);
    }
}
