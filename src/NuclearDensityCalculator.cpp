#include "NuclearDensityCalculator.h"

arma::mat NuclearDensityCalculator::naive_method(const arma::vec &rVals, const arma::vec &zVals) {
    arma::mat result = arma::zeros(rVals.size(), zVals.size()); // number of points on r- and z- axes
    int i = 0, j = 0;
    for (int m = 0; m < basis.mMax; m++) {
        for (int n = 0; n < basis.nMax(m); n++) {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++) {
                i++; //to compute the pos in the rho matrix
                for (int mp = 0; mp < basis.mMax; mp++) {
                    for (int np = 0; np < basis.nMax(mp); np++) {
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++) {
                            j++;//to compute the pos in the rho matrix
                            arma::mat funcA = basis.basisFunc(m, n, n_z, zVals, rVals);
                            arma::mat funcB = basis.basisFunc(mp, np, n_zp, zVals, rVals);
                            //   result += funcA % funcB * rho(m, n, n_z, mp, np, n_zp); // mat += mat % mat * double
                            result += funcA % funcB * imported_rho_values.at(i, j); // Makes sense ?
                        }
                    }
                }
                j = 0; //to compute the pos in the rho matrix
            }
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
}

double NuclearDensityCalculator::rho(int m, int n, int n_z, int mp, int np, int n_zp) {
    return 0;
    //TODO with the imported values;
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
