#include "Basis.h"

arma::mat Basis::basisFunc(int m, int n, int nz, arma::vec &rVec, arma::vec zVec) {
    return nullptr; //TODO
}

arma::vec Basis::zPart(arma::vec &zVec, int nz) {
    return arma::vec(); //TODO
}

arma::vec Basis::rPart(arma::vec &rVec, int m, int n) {
    return arma::vec(); //TODO
}

Basis::Basis(double br, double bz, int N, double Q) {
    initMMax(N, Q);
    initNMax();
    initN_ZMax(N, Q);
    //TODO
}


void Basis::initMMax(int N, double Q) {
    if (Q != 0) {
        this->mMax = std::floor((N + 2) * std::pow(Q, -1.0 / 3.0) - 0.5 * std::pow(Q, -1));
    } else {
        this->mMax = -1;
    }
}

void Basis::initNMax() {
    if (mMax == -1) return;
    this->nMax = arma::regspace<arma::ivec>(mMax - 1, -1, 0) / 2 + 1;
}

void Basis::initN_ZMax(int N, double Q) {
    arma::imat tmp = arma::Mat<arma::sword>(mMax, nMax.at(0), arma::fill::zeros);
    for (int m = 0; m < mMax; m++) {
        for (int n = 0; n < nMax.at(m); n++) {
            tmp.at(m, n) = std::floor((N + 2) * std::pow(Q, 2.0 / 3.0) + 0.5 - (m + 2 * n + 1) * Q);
        }
    }
    tmp.print();
    this->n_zMax = tmp;
}


