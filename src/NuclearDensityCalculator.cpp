#include "NuclearDensityCalculator.h"
#include "Chrono.h"
#include "CallbackResultBuilder.hpp"
#include <list>
#include <omp.h>
#include <memory>

#include "FactorisationFinder.hpp"

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
    int m_a = 0;
    auto builder = std::make_shared<CallbackResultBuilder<arma::mat>>(arma::zeros(rVals.size(), zVals.size()), operation_type::Add);
#pragma omp parallel for default(none) private(m_a) shared(builder, zVals, rVals)
    for (m_a = 0; m_a<basis.mMax; m_a++) {
        std::list<opt2_pair> list;
        for (int n_a = 0; n_a<basis.nMax(m_a); n_a++) {
            for (int n_z_a = 0; n_z_a<basis.n_zMax(m_a, n_a); n_z_a++) {
                list.push_back({n_a, n_z_a});
            }
        }
        std::list<opt2_pair>::iterator a;
#pragma omp parallel default(none) private(a) shared(list, m_a, rVals, zVals, builder) //num_threads(5)
        for (a = list.begin(); a!=list.end(); a++) {
            arma::mat tmp = arma::zeros(rVals.size(), zVals.size());
            auto basisptr = std::make_shared<Basis>(br, bz, N, Q);
#pragma omp single nowait
            for (auto b : list) {
                tmp += basisptr->basisFunc_mem(m_a, b.n, b.nz, rVals, zVals)*rho(m_a, a->n, a->nz, m_a, b.n, b.nz);
            }
            builder->push(basisptr->basisFunc_mem(m_a, a->n, a->nz, rVals, zVals)%tmp);
        }
    }
    return builder->GetResult();
}

/**
 * If we check some hadamard product properties, we can conclude that
 * colA * rowA %  colB * rowB == (colA%colB) * (rowA%rowB) which allows
 * us to do even mode factorisations !
 * @param rVals
 * @param zVals
 * @return
 */
arma::mat NuclearDensityCalculator::optimized_method3(const arma::vec& rVals, const arma::vec& zVals)
{
    Chrono local("optimized_method3");
    FactorisationHelper<struct nuclear_sum_entry, int> nza_factor(nuclear_filter, select_nza);
    /* Rather than making things hard, let's just use the most naive method */
    arma::mat result = arma::zeros(rVals.size(), zVals.size());
    for (int m = 0; m<basis.mMax; m++) {
        for (int n = 0; n<basis.nMax(m); n++) {
            for (int n_z = 0; n_z<basis.n_zMax(m, n); n_z++) {
                for (int mp = 0; mp<basis.mMax; mp++) {
                    for (int np = 0; np<basis.nMax(mp); np++) {
                        for (int n_zp = 0; n_zp<basis.n_zMax(mp, np); n_zp++) {
                            nza_factor.add({m, n, n_z, mp, np, n_zp});
                        }
                    }
                }
            }
        }
    }

    auto nza_entries = nza_factor.get_factored();
    for (auto& nza_entry : nza_entries) {
        /* nza_zpart is the loop constant */
        arma::rowvec nza_zpart = basis.zPart_mem(zVals, nza_entry.common).as_row();
        arma::mat tmp = arma::zeros(rVals.size(), zVals.size());
        /* We factor out nzb of the sum to do */
        FactorisationHelper<struct nuclear_sum_entry, int> nzb_Factor(nza_entry.to_sum, nuclear_filter, select_nzb);
        auto nzb_entries = nzb_Factor.get_factored();

        for (auto& nzb_entry : nzb_entries) {
            /* nzb_zpart is the loop constant */
            arma::rowvec nzb_zpart = basis.zPart_mem(zVals, nzb_entry.common).as_row();
            arma::colvec all_rpart = arma::zeros(rVals.size());
            /* We factor out the pair ma_na of the sum that is left to compute */
            FactorisationHelper<struct nuclear_sum_entry, struct m_n_pair> ma_na_factor(nzb_entry.to_sum, nuclear_filter, select_ma_na);
            auto ma_na_entries = ma_na_factor.get_factored();

            for (auto& ma_na_entry : ma_na_entries) {
                /* mana_rpart is the loop constant */
                arma::colvec mana_rpart = basis.rPart_mem(rVals, ma_na_entry.common.m_a, ma_na_entry.common.n_a).as_col();
                arma::colvec mbnb_rpart = arma::zeros(rVals.size());
                /* We could factor out the pair mb and nb but its useless */

                for (auto& e : ma_na_entry.to_sum) {
                    mbnb_rpart += basis.rPart_mem(rVals, e.m_b, e.n_b).as_col()*rho(e.m_a, e.n_a, e.nz_a, e.m_b, e.n_b, e.nz_b);
                }
                all_rpart += mbnb_rpart%mana_rpart;
            }
            tmp += all_rpart*nzb_zpart;
        }
        result += tmp%(arma::colvec(rVals.size(), arma::fill::ones)*nza_zpart);
    }
    return result;
}

/**
 *
 */
NuclearDensityCalculator::NuclearDensityCalculator()
{
    imported_rho_values.load("src/rho.arma", arma::arma_ascii);

#ifdef DEBUG
    std::cout << "[src/rho.arma defs imported]" << std::endl;
#endif
    basis = Basis(br, bz, N, Q);

    ind = arma::Cube<arma::sword>(basis.n_zMax.max(), basis.nMax.max(), basis.mMax);
    int i = 0;
    for (int m = 0; m<basis.mMax; m++) {
        for (int n = 0; n<basis.nMax.at(m); n++) {
            for (int n_z = 0; n_z<basis.n_zMax.at(m, n); n_z++) {
                ind.at(n_z, n, m) = i;
                i++;
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
double NuclearDensityCalculator::rho(int m, int n, int n_z, int mp, int np, int n_zp)
{
    auto a = ind.at(n_z, n, m);
    auto b = ind.at(n_zp, np, mp);
    return imported_rho_values.at(a, b);
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
