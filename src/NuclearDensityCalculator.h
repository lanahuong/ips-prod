#ifndef PROJET_IPS1_NUCLEARDENSITYCALCULATOR_H
#define PROJET_IPS1_NUCLEARDENSITYCALCULATOR_H

#include "Basis.h"


class NuclearDensityCalculator {

private:
    int N = 14;
    double Q = 1.3;
    double br = 1.935801664793151;
    double bz = 1.829683956491218;
    arma::mat imported_rho_values;
    Basis basis;

    static double rho(int m, int n, int n_z, int mp, int np, int n_zp);

public:
    NuclearDensityCalculator();

    /**
     * @see https://dubrayn.github.io/IPS-PROD/project.html#19
     */
    void printRhoDefs();

    arma::icube ind;

    arma::mat naive_method(const arma::vec &rVals, const arma::vec &zVals);
    arma::mat naive_method2(const arma::vec &rVals, const arma::vec &zVals);

};

#endif //PROJET_IPS1_NUCLEARDENSITYCALCULATOR_H
