#include "Basis.h"
#include "Poly.h"
#include "constants.h"

Basis::Basis(double BR, double BZ, int N, double Q)
        :br(BR), bz(BZ)
{
    this->mMax = calcMMax(N, Q);
    this->nMax = calcNMax();
    this->n_zMax = calcN_zMax(N, Q);
}

/**
 * Special constructor to use if we want to accelerate computing by
 * pre-computing some values and storing the computed ones.
 * @param rVals
 * @param zVals
 */
Basis::Basis(double BR, double BZ, int N, double Q, const arma::vec& rVals, const arma::vec& zVals)
        :Basis(BR, BZ, N, Q)
{
    rvec_mem = rVals;
    squared_rarg_mem = arma::square(rVals/br);
    rexp_mem = arma::exp(-squared_rarg_mem/2.0);
    computed_r_vals = std::vector<arma::vec>((mMax+1)*(nMax.max()+1));
    computed_r_indices = std::vector<bool>((mMax+1)*(nMax.max()+1), false);

    zvec_mem = zVals;
    squared_zarg_mem = arma::square(zVals/bz);
    zexp_mem = arma::exp(-squared_zarg_mem/2.0);
    computed_z_vals = std::vector<arma::vec>(n_zMax.max());
    computed_z_indices = std::vector<bool>(n_zMax.max(), false);

    is_mem = true;
}

/**
 * The basis function from the lesson as is
 * @param m
 * @param n
 * @param nz
 * @param rVec
 * @param zVec
 * @return
 */
arma::mat Basis::basisFunc(int m, int n, int nz, const arma::vec& rVec, const arma::vec& zVec)
{
    return rPart(rVec, m, n).as_col()*zPart(zVec, nz).as_row();
}

arma::vec Basis::zPart(const arma::vec& zVec, int nz, bool use_mem)
{
    double const_factor = pow(bz, -0.5)*pow(PI, -0.25);

    for (int i = 1; i<=nz; i++) {
        const_factor *= pow(2*i, -0.5);
    }

    arma::vec squared_arg = use_mem ? squared_zarg_mem : arma::square(zVec/bz);
    arma::vec exp = use_mem ? zexp_mem : arma::exp(-squared_arg/2.0);
    poly.calcHermite(nz+1, zVec/bz);

    return const_factor*exp%poly.hermite(nz);
}

arma::vec Basis::rPart(const arma::vec& rVec, int m, int n, bool use_mem)
{
    double const_factor = pow(br, -1)*pow(PI, -0.5);

    for (int i = n+1; i<=m+n; i++) {
        const_factor *= pow(i, -0.5);
    }

    arma::vec squared_arg = use_mem ? squared_rarg_mem : arma::square(rVec/br);
    arma::vec exp = use_mem ? rexp_mem : arma::exp(-squared_arg/2.0);
    arma::vec pow = arma::pow(rVec/br, m);
    poly.calcLaguerre(m+1, n+1, squared_arg);

    return const_factor*exp%pow%poly.laguerre(m, n);
}

int Basis::calcMMax(int N, double Q)
{
    if (Q!=0) {
        return static_cast<int>(floor((N+2)*pow(Q, -1.0/3.0)-0.5*pow(Q, -1)));
    }
    return 0;
}

arma::ivec Basis::calcNMax() const
{
    return arma::regspace<arma::ivec>(mMax-1, -1, 0)/2+1;
}

arma::imat Basis::calcN_zMax(int N, double Q)
{
    arma::imat tmp(mMax, nMax.at(0), arma::fill::zeros);
    for (int m = 0; m<mMax; m++) {
        for (int n = 0; n<nMax.at(m); n++) {
            tmp.at(m, n) = floor((N+2)*pow(Q, 2.0/3.0)+0.5-(m+2*n+1)*Q);
        }
    }
    return tmp;
}

arma::mat Basis::basisFunc_mem(int m, int n, int nz)
{
    return rPart_mem(m, n).as_col()*zPart_mem(nz).as_row();
}

/**
 * Returns the value if it was alreadu computed else computes it and
 * saves it.
 */
arma::vec Basis::zPart_mem(int nz)
{
    if (computed_z_indices.at(nz)) {
        return computed_z_vals.at(nz);
    }
    else {
        arma::vec tmp(zPart(zvec_mem, nz, true));
        computed_z_indices.at(nz) = true;
        computed_z_vals.at(nz) = tmp;
        return tmp;
    }
}

/**
 * Returns the already computed value else computes it
 */
arma::vec Basis::rPart_mem(int m, int n)
{
    if (computed_r_indices.at(n*(mMax+1)+m)) {
        return computed_r_vals.at(n*(mMax+1)+m);
    }
    else {
        arma::vec tmp(rPart(rvec_mem, m, n, true));
        computed_r_indices.at(n*(mMax+1)+m) = true;
        computed_r_vals.at(n*(mMax+1)+m) = tmp;
        return tmp;
    }
}



