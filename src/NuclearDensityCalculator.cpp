#include "NuclearDensityCalculator.h"
#include "Chrono.h"
#include <list>

arma::mat NuclearDensityCalculator::naive_method(const arma::vec &rVals, const arma::vec &zVals) {
    Chrono local("naive_method");
    arma::mat result = arma::zeros(rVals.size(), zVals.size()); // number of points on r- and z- axes
    for (int m = 0; m < basis.mMax; m++) {
        for (int n = 0; n < basis.nMax(m); n++) {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++) {
                for (int mp = 0; mp < basis.mMax; mp++) {
                    for (int np = 0; np < basis.nMax(mp); np++) {
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++) {
                            arma::mat funcA = basis.basisFunc(m, n, n_z, rVals, zVals);
                            arma::mat funcB = basis.basisFunc(mp, np, n_zp, rVals, zVals);
                            result += funcA % funcB * rho(m, n, n_z, mp, np, n_zp);
                        }
                    }
                }
            }
        }
    }
    return result;
}

/* Since results only increases when ma = mb */
arma::mat NuclearDensityCalculator::optimized_method1(const arma::vec &rVals, const arma::vec &zVals) {
    Chrono local("optimized_method1");
    arma::mat result = arma::zeros(rVals.size(), zVals.size()); // number of points on r- and z- axes
    for (int m_a = 0; m_a < basis.mMax; m_a++) {
        for (int n_a = 0; n_a < basis.nMax(m_a); n_a++) {
            for (int n_z_a = 0; n_z_a < basis.n_zMax(m_a, n_a); n_z_a++) {
                for (int n_b = 0; n_b < basis.nMax(m_a); n_b++) {
                    for (int n_z_b = 0; n_z_b < basis.n_zMax(m_a, n_b); n_z_b++) {
                        arma::mat funcA = basis.basisFunc(m_a, n_a, n_z_a, rVals, zVals);
                        arma::mat funcB = basis.basisFunc(m_a, n_b, n_z_b, rVals, zVals);
                        result += funcA % funcB * rho(m_a, n_a, n_z_a, m_a, n_b, n_z_b);
                    }
                }
            }
        }
    }
    return result;
}


arma::mat NuclearDensityCalculator::optimized_method2(const arma::vec &rVals, const arma::vec &zVals) {
    Chrono local("optimized_method2");
    arma::mat result = arma::zeros(rVals.size(), zVals.size()); // number of points on r- and z- axes
    struct opt2_pair {
        int n, nz;
    };

    for (int m_a = 0; m_a < basis.mMax; m_a++) {
        std::list<opt2_pair> list;
        for (int n_a = 0; n_a < basis.nMax(m_a); n_a++) {
            for (int n_z_a = 0; n_z_a < basis.n_zMax(m_a, n_a); n_z_a++) {
                list.push_back({n_a, n_z_a});
            }
        }
        for (auto a : list) {
            arma::mat tmp = arma::zeros(rVals.size(), zVals.size());
            for (auto b : list) {
                tmp += basis.basisFunc_mem(m_a, b.n, b.nz, rVals, zVals) * rho(m_a, a.n, a.nz, m_a, b.n, b.nz);
            }
            result += basis.basisFunc_mem(m_a, a.n, a.nz, rVals, zVals) % tmp;
        }

    }
    return result;
}


NuclearDensityCalculator::NuclearDensityCalculator() {
    imported_rho_values.load("src/rho.arma", arma::arma_ascii);

#ifdef DEBUG
    std::cout << "[src/rho.arma defs imported]" << std::endl;
#endif
    basis = Basis(br, bz, N, Q);

    ind = arma::Cube<arma::sword>(basis.n_zMax.max(), basis.nMax.max(), basis.mMax);
    int i = 0;
    for (int m = 0; m < basis.mMax; m++) {
        for (int n = 0; n < basis.nMax.at(m); n++) {
            for (int n_z = 0; n_z < basis.n_zMax.at(m, n); n_z++) {
                ind.at(n_z, n, m) = i;
                i++;
            }
        }
    }
}


double NuclearDensityCalculator::rho(int m, int n, int n_z, int mp, int np, int n_zp) {
    int a = ind.at(n_z, n, m);
    int b = ind.at(n_zp, np, mp);
    return imported_rho_values.at(a, b);
}

void NuclearDensityCalculator::printRhoDefs() {
    uint i = 0;
    for (int m = 0; m < basis.mMax; m++) {
        for (int n = 0; n < basis.nMax(m); n++) {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++) {
                std::cout << "Basis vector " << i << ": m=" << m << " n=" << n << " n_z=" << n_z << std::endl;
                i++;
            }
        }
    }
}
