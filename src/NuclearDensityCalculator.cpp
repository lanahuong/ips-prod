#include <vector>
#include <memory>

#include "NuclearDensityCalculator.h"
#include "Chrono.hpp"
#include "ThreadSafeAccumulator.hpp"
#include "FactorisationHelper.hpp"


arma::mat NuclearDensityCalculator::naive_method(const arma::vec& rVals, const arma::vec& zVals)
{
    Chrono local("naive_method");
    arma::mat result = arma::zeros(rVals.size(), zVals.size()); // number of points on r- and z- axes
    for (int m = 0; m<basis.mMax; m++) {
        for (int n = 0; n<basis.nMax(m); n++) {
            for (int n_z = 0; n_z<basis.n_zMax(m, n); n_z++) {
                for (int mp = 0; mp<basis.mMax; mp++) {
                    for (int np = 0; np<basis.nMax(mp); np++) {
                        for (int n_zp = 0; n_zp<basis.n_zMax(mp, np); n_zp++) {
                            arma::mat funcA = basis.basisFunc(m, n, n_z, rVals, zVals);
                            arma::mat funcB = basis.basisFunc(mp, np, n_zp, rVals, zVals);
                            result += funcA%funcB*rho(m, n, n_z, mp, np, n_zp);
                        }
                    }
                }
            }
        }
    }
    return result;
}

/* Since results only increases when ma = mb */
arma::mat NuclearDensityCalculator::optimized_method1(const arma::vec& rVals, const arma::vec& zVals)
{
    Chrono local("optimized_method1");
    arma::mat result = arma::zeros(rVals.size(), zVals.size()); // number of points on r- and z- axes
    for (int m_a = 0; m_a<basis.mMax; m_a++) {
        for (int n_a = 0; n_a<basis.nMax(m_a); n_a++) {
            for (int n_z_a = 0; n_z_a<basis.n_zMax(m_a, n_a); n_z_a++) {
                for (int n_b = 0; n_b<basis.nMax(m_a); n_b++) {
                    for (int n_z_b = 0; n_z_b<basis.n_zMax(m_a, n_b); n_z_b++) {
                        arma::mat funcA = basis.basisFunc(m_a, n_a, n_z_a, rVals, zVals);
                        arma::mat funcB = basis.basisFunc(m_a, n_b, n_z_b, rVals, zVals);
                        result += funcA%funcB*rho(m_a, n_a, n_z_a, m_a, n_b, n_z_b);
                    }
                }
            }
        }
    }
    return result;
}


arma::mat NuclearDensityCalculator::optimized_method2(const arma::vec& rVals, const arma::vec& zVals)
{
    Chrono local("optimized_method2");
    struct opt2_pair {
      int n, nz;
    };
    arma::mat builder = arma::zeros(rVals.size(), zVals.size());
    for (int m_a = 0; m_a<basis.mMax; m_a++) {
        std::list<opt2_pair> list;
        for (int n_a = 0; n_a<basis.nMax(m_a); n_a++) {
            for (int n_z_a = 0; n_z_a<basis.n_zMax(m_a, n_a); n_z_a++) {
                list.push_back({n_a, n_z_a});
            }
        }
        for (auto a = list.begin(); a!=list.end(); a++) {
            arma::mat tmp = arma::zeros(rVals.size(), zVals.size());
            for (auto b : list) {
                tmp += basis.basisFunc(m_a, b.n, b.nz, rVals, zVals)*rho(m_a, a->n, a->nz, m_a, b.n, b.nz);
            }
            builder += basis.basisFunc(m_a, a->n, a->nz, rVals, zVals)%tmp;
        }
    }
    return builder;
}

/**
 * One can prove the following property:
 * colA * rowA %  colB * rowB == (colA%colB) * (rowA%rowB) which allows
 * us to do even mode factorisations !
 */
arma::mat NuclearDensityCalculator::optimized_method3(const arma::vec& rVals, const arma::vec& zVals) const
{
    Chrono local("optimized_method3");
    FactorisationHelper<struct quantum_numbers, int> nza_factored_sum(select_nza, symmetry_filter);
    for (int m_a(0), m_amax(basis.mMax); m_a<m_amax; m_a++) {  /* Rather than making things hard, let's just use the most naive method */
        for (int n_a(0), n_amax(basis.nMax(m_a)); n_a<n_amax; n_a++) {
            for (int nz_a(0), nz_amax(basis.n_zMax(m_a, n_a)); nz_a<nz_amax; nz_a++) {
                for (int n_b(0), m_bmax(basis.nMax(m_a)); n_b<m_bmax; n_b++) {
                    for (int nz_b(0), nz_bmax(basis.n_zMax(m_a, n_b)); nz_b<nz_bmax; nz_b++) {
                        nza_factored_sum.add({m_a, n_a, nz_a, m_a, n_b, nz_b, 1});
                    }
                }
            }
        }
    }

    const int zSize(zVals.size()), rSize(rVals.size());
    std::shared_ptr<ThreadSafeAccumulator<arma::mat>> builder(std::make_shared<ThreadSafeAccumulator<arma::mat >>(arma::zeros(rSize, zSize), operation_type::Add));
    const arma::colvec unit(rSize, arma::fill::ones);
    const Basis basis_mem(br, bz, N, Q, rVals, zVals);
    /* nza_zpart is the loop constant */
#pragma omp parallel for default(shared)
    for (size_t i = 0; i<nza_factored_sum.size(); i++) {
        Basis basis_local(basis_mem); /* We copy the basis in each thread, strangely the private openMP dont work */
        const arma::rowvec nza_zpart(basis_local.zPart_mem(nza_factored_sum[i].factor).as_row());
        arma::mat tmp(arma::zeros(rSize, zSize));
        /* We factor out nzb of the sum left to compute */
        FactorisationHelper<quantum_numbers, int> nzb_factored(nza_factored_sum[i].factored_out, select_nzb);
        /* nzb_zpart is the loop constant */
        for (factored<quantum_numbers, int>& nzb_term : nzb_factored) {
            const arma::rowvec nzb_zpart(basis_local.zPart_mem(nzb_term.factor).as_row());
            arma::colvec all_rpart(arma::zeros(rSize));
            /* We factor out the pair ma na of the sum that is left to compute */
            FactorisationHelper<quantum_numbers, m_n_pair> mana_factored(nzb_term.factored_out, select_mana);
            for (const factored<quantum_numbers, m_n_pair>& mana_term : mana_factored) {
                /* mana_rpart is the loop constant */
                const arma::colvec mana_rpart(basis_local.rPart_mem(mana_term.factor.m_a, mana_term.factor.n_a));
                arma::colvec mbnb_rpart(arma::zeros(rSize));
                /* We could factor out the pair mb and nb but its useless */
                for (const quantum_numbers e : mana_term.factored_out) {
                    mbnb_rpart += basis_local.rPart_mem(e.m_b, e.n_b)*(rho(e.m_a, e.n_a, e.nz_a, e.m_b, e.n_b, e.nz_b)*e.count);
                }
                all_rpart += mbnb_rpart%mana_rpart;
            }
            tmp += all_rpart*nzb_zpart;
        }
        builder->push(tmp%(unit*nza_zpart));
    }
    return builder->GetResult(); /* Computes pending operations and returns the result of the accumulator */
}

/**
 *
 */
NuclearDensityCalculator::NuclearDensityCalculator()
        :basis(br, bz, N, Q)
{
    imported_rho_values.load("src/rho.arma", arma::arma_ascii);
#ifdef DEBUG
    std::cout << "[src/rho.arma defs imported]" << std::endl;
#endif
    ind = arma::Cube<arma::sword>(basis.n_zMax.max(), basis.nMax.max(), basis.mMax);
    int i = 0;
    for (int m = 0; m<basis.mMax; m++) {
        for (int n = 0; n<basis.nMax.at(m); n++) {
            for (int n_z = 0; n_z<basis.n_zMax.at(m, n); n_z++) {
                ind.at(n_z, n, m) = i++;
            }
        }
    }
}

/**
 *
 * @param m
 * @param n
 * @param n_z
 * @param mp
 * @param np
 * @param n_zp
 * @return
 */
inline double NuclearDensityCalculator::rho(int m, int n, int n_z, int mp, int np, int n_zp) const
{
    return imported_rho_values.at(ind.at(n_z, n, m), ind.at(n_zp, np, mp));
}

void NuclearDensityCalculator::printRhoDefs()
{
    uint i = 0;
    for (int m = 0; m<basis.mMax; m++) {
        for (int n = 0; n<basis.nMax(m); n++) {
            for (int n_z = 0; n_z<basis.n_zMax(m, n); n_z++) {
                std::cout << "Basis vector " << i << ": m=" << m << " n=" << n << " n_z=" << n_z << std::endl;
                i++;
            }
        }
    }
}


arma::cube  NuclearDensityCalculator::density_cartesian(const int xyPoints, const int zPoints, const arma::vec& rVals, const arma::mat& res) {
    // Create an empty cube of the correct size
    arma::cube cube(xyPoints, xyPoints, zPoints, arma::fill::zeros);
    // For each point (x,y) compute the corresponding radius, find it in the rVals list
    // and put all the values with variating z in the cube
    for (int x = 0 ; x < xyPoints ; x++) {
        double x_real = rVals(x);
        for (int y = 0 ; y < xyPoints ; y++) {
            double y_real = rVals(y);
            double r = sqrt(x_real * x_real + y_real * y_real);
            arma::vec r_diff = arma::abs(rVals - r);
            arma::uword k = r_diff.index_min();
            
            cube.tube(x,y) = res.row(k);
        }
    }
    return cube;
}
