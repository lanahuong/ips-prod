#include "NuclearDensityCalculator.h"
#include "Chrono.h"
#include "ThreadSafeAccumulator.hpp"
#include <list>
#include <memory>

#include "FactorisationHelper.hpp"

/**
 *
 * @param rVals
 * @param zVals
 * @return
 */
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

/**
 *
 * @param rVals
 * @param zVals
 * @return
 */
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
 * @param rVals
 * @param zVals
 * @return
 */
arma::mat NuclearDensityCalculator::optimized_method3(const arma::vec& rVals, const arma::vec& zVals)
{
    Chrono local("optimized_method3");
    const int zSize(zVals.size());
    const int rSize(rVals.size());
    FactorisationHelper<struct quantum_numbers, int> nza_factor(select_nza, nuclear_filter);
    /* Rather than making things hard, let's just use the most naive method */
    for (int m = 0; m<basis.mMax; m++) {
        for (int n = 0; n<basis.nMax(m); n++) {
            for (int n_z = 0; n_z<basis.n_zMax(m, n); n_z++) {
                for (int np = 0; np<basis.nMax(m); np++) {
                    for (int n_zp = 0; n_zp<basis.n_zMax(m, np); n_zp++) {
                        if (np<n && n_zp<n_z) {
                            nza_factor.add({m, n, n_z, m, np, n_zp, 2});
                        }
                        else if (np<=n || n_zp<=n_z) {
                            nza_factor.add({m, n, n_z, m, np, n_zp, 1});
                        }
                    }
                }
            }
        }
    }
    /* Thread safe accumulator to which we can send the results once we stop using openMP
     * Is has its own buffer when the acc is locked and a spinlock mechanism as the operations
     * are expected to be quite fast and not happen a lot */
    std::shared_ptr<ThreadSafeAccumulator<arma::mat>> builder(std::make_shared<ThreadSafeAccumulator<arma::mat>>(arma::zeros(rSize, zSize), operation_type::Add));
    const arma::colvec unit(rSize, arma::fill::ones);
    const std::vector<factored<quantum_numbers, int>> nza_entries(nza_factor.get_vfactored());
    const Basis basis_shared(br, bz, N, Q, rVals, zVals); /* Its easier if each thread has its own Basis class */
#pragma omp parallel for default(none) shared(nza_entries, builder, unit, basis_shared, zSize, rSize)
    /* nza_zpart is the loop constant */
    for (const auto& nza_entry : nza_entries) {
        Basis basis_local(basis_shared); /* We copy the basis in each thread, strangely the private openMP dont work */
        const arma::rowvec nza_zpart(basis_local.zPart_mem(nza_entry.factor).as_row());
        arma::mat tmp(arma::zeros(rSize, zSize));
        /* We factor out nzb of the sum left to compute */
        FactorisationHelper<quantum_numbers, int> nzb_Factor(nza_entry.factored_out, select_nzb, nuclear_filter);
        const std::list<factored<quantum_numbers, int>> nzb_entries(nzb_Factor.get_factored());
        /* nzb_zpart is the loop constant */
        for (const factored<quantum_numbers, int>& nzb_entry : nzb_entries) {
            const arma::rowvec nzb_zpart(basis_local.zPart_mem(nzb_entry.factor).as_row());
            arma::colvec all_rpart(arma::zeros(rSize));
            /* We factor out the pair ma na of the sum that is left to compute */
            FactorisationHelper<quantum_numbers, m_n_pair> ma_na_factor(nzb_entry.factored_out, select_ma_na, nuclear_filter);
            const std::list<factored<quantum_numbers, m_n_pair>> ma_na_entries(ma_na_factor.get_factored());
            /* mana_rpart is the loop constant */
            for (const factored<quantum_numbers, m_n_pair>& ma_na_entry : ma_na_entries) {
                const arma::colvec mana_rpart(basis_local.rPart_mem(ma_na_entry.factor.m_a, ma_na_entry.factor.n_a));
                arma::colvec mbnb_rpart(arma::zeros(rSize));
                /* We could factor out the pair mb and nb but its useless */
                for (const quantum_numbers& e : ma_na_entry.factored_out) {
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

/**
 *
 */
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
