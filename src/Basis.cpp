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
    findMMax(N, Q);

    //TODO
}


void Basis::findMMax(int N, double Q) {
    if (Q != 0) {
        this->mMax = std::floor((N + 2) * std::pow(Q, -1.0 / 3.0) - 0.5 * std::pow(Q, -1));
    } else {
        this->mMax = -1;
    }
}


