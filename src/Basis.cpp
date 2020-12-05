#include "Basis.h"
#include "Poly.h"
#include "constants.h"


arma::mat Basis::basisFunc(int m, int n, int nz, const arma::vec &rVec, const arma::vec &zVec) {
    return rPart(rVec, m, n).as_col() * zPart(zVec, nz).as_row();
}

arma::mat Basis::basisFunc_mem(int m, int n, int nz, const arma::vec &rVec, const arma::vec &zVec) {
    return rPart_mem(rVec, m, n).as_col() * zPart_mem(zVec, nz).as_row();
}


arma::vec Basis::zPart(const arma::vec &zVec, int nz) {
    double const_factor = pow(bz, -0.5) * pow(PI, -0.25);

    for (int i = 1; i <= nz; i++) {
        const_factor *= pow(2 * i, -0.5);
    }

    arma::vec squared_arg = arma::square(zVec / bz);
    arma::vec exp = arma::exp(-squared_arg / 2.0);
    poly.calcHermite(nz + 1, zVec / bz);

    return const_factor * exp % poly.hermite(nz);
}

arma::vec Basis::rPart(const arma::vec &rVec, int m, int n) {
    double const_factor = pow(br, -1) * pow(PI, -0.5);

    for (int i = n + 1; i <= m + n; i++) {
        const_factor *= pow(i, -0.5);
    }

    arma::vec squared_arg = arma::square(rVec / br);
    arma::vec exp = arma::exp(-squared_arg / 2.0);
    arma::vec pow = arma::pow(rVec / br, m);
    poly.calcLaguerre(m + 1, n + 1, squared_arg);

    return const_factor * exp % pow % poly.laguerre(m, n);
}

Basis::Basis(double BR, double BZ, int N, double Q) :br(BR), bz(BZ) {
    this->mMax = calcMMax(N, Q);
    this->nMax = calcNMax();
    this->n_zMax = calcN_zMax(N, Q);

    computed_z_indices = arma::ivec(n_zMax.max(), arma::fill::zeros);
    computed_z_vals = std::vector<arma::vec>(n_zMax.max());

    computed_r_indices = arma::imat(mMax + 1, nMax.max() + 1, arma::fill::zeros);
    computed_r_vals = std::vector<arma::vec>((mMax + 1) * (nMax.max() + 1));
}


int Basis::calcMMax(int N, double Q) {
    if (Q != 0) {
        return static_cast<int>(floor((N + 2) * pow(Q, -1.0 / 3.0) - 0.5 * pow(Q, -1)));
    }
    return 0;
}

arma::ivec Basis::calcNMax() const {
    return arma::regspace<arma::ivec>(mMax - 1, -1, 0) / 2 + 1;
}

arma::imat Basis::calcN_zMax(int N, double Q) {
    arma::imat tmp = arma::Mat<arma::sword>(mMax, nMax.at(0), arma::fill::zeros);
    for (int m = 0; m < mMax; m++) {
        for (int n = 0; n < nMax.at(m); n++) {
            tmp.at(m, n) = floor((N + 2) * pow(Q, 2.0 / 3.0) + 0.5 - (m + 2 * n + 1) * Q);
        }
    }
    return tmp;
}

arma::vec Basis::zPart_mem(const arma::vec &zVec, int nz) {
    if (computed_z_indices.at(nz) == reinterpret_cast<long int> ( zVec.memptr())){
        return computed_z_vals.at(nz);
    } else {
        arma::vec tmp = zPart(zVec, nz);
        computed_z_indices.at(nz) = reinterpret_cast<long int> ( zVec.memptr());
        computed_z_vals.at(nz) = tmp;
        return tmp;
    }
}

arma::vec Basis::rPart_mem(const arma::vec &rVec, int m, int n) {
    if (computed_r_indices.at(m, n) == reinterpret_cast<long int> ( rVec.memptr())) {
        return computed_r_vals.at(n * (mMax + 1) + m);
    } else {
        arma::vec tmp = rPart(rVec, m, n);
        computed_r_indices.at(m, n) = reinterpret_cast<long int> ( rVec.memptr());
        computed_r_vals.at(n * (mMax + 1) + m) = tmp;
        return tmp;
    }
}


Basis::Basis() = default;
