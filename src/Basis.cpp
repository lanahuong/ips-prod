#include "Basis.h"
#include "Poly.h"
#include "constants.h"


arma::mat Basis::basisFunc(int m, int n, int nz, arma::vec &rVec, arma::vec zVec) {
    return nullptr; //TODO
}

arma::vec Basis::zPart(arma::vec &zVec, int nz) {
    return arma::vec(); //TODO
}

arma::vec Basis::rPart(arma::vec &rVec, int m, int n) const {
    int abs_m = floor(m);
    long double const_factor = pow(br, -1) * pow(PI, -0.5);
    for (int i = n + 1; i <= abs_m + n; i++) {
        const_factor *= pow(i, -0.5);
    }

    arma::vec squared_arg = arma::square(rVec / br);

    Poly poly;
    poly.calcLaguerre(abs_m + 1, n + 1, squared_arg);

    arma::vec exp = arma::exp(-squared_arg / 2.0);
    arma::vec pow = arma::pow(rVec / br, abs_m);
    return const_factor * exp % pow % poly.laguerre(abs_m, n);
}

Basis::Basis(double br, double bz, int N, double Q) : br(br), bz(bz) {
    this->mMax = calcMMax(N, Q);
    this->nMax = calcNMax();
    this->n_zMax = calc(N, Q);
    //TODO
}


int Basis::calcMMax(int N, double Q) {
    if (Q != 0) {
        return floor((N + 2) * pow(Q, -1.0 / 3.0) - 0.5 * pow(Q, -1));
    } else {
        return 0;
    }
}

arma::ivec Basis::calcNMax() const {
    return arma::regspace<arma::ivec>(mMax - 1, -1, 0) / 2 + 1;
}

arma::imat Basis::calc(int N, double Q) {
    arma::imat tmp = arma::Mat<arma::sword>(mMax, nMax.at(0), arma::fill::zeros);
    for (int m = 0; m < mMax; m++) {
        for (int n = 0; n < nMax.at(m); n++) {
            tmp.at(m, n) = floor((N + 2) * pow(Q, 2.0 / 3.0) + 0.5 - (m + 2 * n + 1) * Q);
        }
    }
    return tmp;
}


