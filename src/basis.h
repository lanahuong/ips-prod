#include <armadillo>

class Basis
{
    public:

        int mMax;
        arma::ivec nMax;
        arma::imat n_zMax;

        Basis(double br, double bz, int N, double Q);
        arma::vec rPart(arma::vec &rVec, int m, int n);
        arma::vec zPart(arma::vec &zVec, int nz);
        arma::mat basisFunc(int m, int n, int nz, arma::vec &rVec, arma::vec zVec);

    private:
        
        double br;
        double bz;
        arma::cube rPartVals;
        arma::vec lastRVals;
        arma::mat zPartVals;
        arma::vec lastZVals;
}