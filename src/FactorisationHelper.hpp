/**
 * @file FactorisationHelper.hpp
 */

#ifndef PROJET_IPS1_FACTORISATIONFINDER_H
#define PROJET_IPS1_FACTORISATIONFINDER_H

#include <list>
#include <algorithm>

/**
 * Struct to hold the values of a factorisation
 * @tparam fa the type of what we factor out (can be a pair of values, of whatever)
 * @tparam T the type of the entries we added
 */
template<typename T, typename fa>
struct factored {
  fa factor;
  std::vector<T> factored_out;
  //std::list<T> factored_out;
};

/**
 * @class FactorisationHelper
 * @tparam fa the type of what we factor out (can be a pair of values, of whatever)
 * @tparam T the type of the entries we added
 */
template<typename T, typename f>
class FactorisationHelper {
public:
    typedef bool(* input_filter)(T& a);
    typedef f (* selector_function)(const T& a);

    /**
     * Constructor of the factorisation helper.
     * @param filter to filter some values out
     * @param selector to select the values we want to factor out;
     */
    explicit FactorisationHelper(selector_function selector, input_filter filter = [](T& a) { return true; });

    /**
     * Constructor that imports a list of entries
     * @param input list of type T entries
     * @param filter to filter some values out
     * @param selector to select the values we want to factor out;
     */
    FactorisationHelper(std::list<T> input, selector_function selector, input_filter filter = [](T& a) { return true; });
    FactorisationHelper(std::vector<T> input, selector_function selector, input_filter filter = [](T& a) { return true; });

    /**
     * Adds an entry to the factoriser
     * @param entry
     */
    inline void add(T entry);

    /**
     * Called at the end to get the factorisation result
     * @return a list of structs. Each struct contains a "factor" which is the key that was factored out
     * of the values that are in the list "factored_out"
     */
    std::list<struct factored<T, f>> get_factored();

    std::vector<struct factored<T, f>> get_vfactored();
private:
    const input_filter filter;
    const selector_function selector;
    std::list<struct factored<T, f>> out{};

    /**
     * As we add entries, they are placed in the right categories
     * to factor them out given a selector that "selects" a peculiar value in our entry
     * @tparam fa the type of what we factor out (can be a pair of values, of whatever)
     * @tparam T the type of the entries we added
     * @param entry
     */
    inline void dispatch_entry(T entry);
    void remove_duplicates();
    void remove_symetric_elts();
};

/**
 * As we add entries, they are placed in the right categories
 * to factor them out given a selector that "selects" a peculiar value in our entry
 * @tparam fa the type of what we factor out (can be a pair of values, of whatever)
 * @tparam T the type of the entries we added
 * @param entry
 */
template<typename T, typename f>
inline void FactorisationHelper<T, f>::dispatch_entry(T entry)
{
    f fac = selector(entry);
    auto it = std::find_if(out.begin(), out.end(), [&fac](const struct factored<T, f>& x) { return x.factor==fac; });
    if (it!=out.end()) {
        it->factored_out.push_back(entry);
    }
    else {
        out.push_back({fac, std::vector<T>{entry}});
    //    out.push_back({fac, std::list<T>{entry}});
    }
}

/**
 *
 * @tparam T
 * @tparam f
 * @param selec
 * @param filt
 */
template<typename T, typename f>
FactorisationHelper<T, f>::FactorisationHelper(FactorisationHelper::selector_function selec, FactorisationHelper::input_filter filt)
        :filter(filt), selector(selec) {}

/**
 *
 * @tparam T
 * @tparam f
 * @param input
 * @param select
 * @param filt
 */
template<typename T, typename f>
FactorisationHelper<T, f>::FactorisationHelper(std::list<T> input, FactorisationHelper::selector_function select, FactorisationHelper::input_filter filt)
        : FactorisationHelper<T, f>::FactorisationHelper(select, filt) {
    for (auto &in : input) {
        add(in);
    }
}

template<typename T, typename f>
FactorisationHelper<T, f>::FactorisationHelper(std::vector<T> input, FactorisationHelper::selector_function select, FactorisationHelper::input_filter filt)
        : FactorisationHelper<T, f>::FactorisationHelper(select, filt) {
    for (auto &in : input) {
        add(in);
    }
}


template<typename T, typename f>
inline void FactorisationHelper<T, f>::add(T entry) {
    if (filter(entry)) {
        dispatch_entry(entry);
    }
}

template<typename T, typename f>
std::list<struct factored<T, f>> FactorisationHelper<T, f>::get_factored()
{
    remove_duplicates();
    return out;
}

/**
 *
 * @tparam T
 * @tparam f
 */
template<typename T, typename f>
void FactorisationHelper<T, f>::remove_duplicates() {
    //TODO ? Same huge impact :/
}

template<typename T, typename f>
std::vector<factored<T, f>> FactorisationHelper<T, f>::get_vfactored() {
    std::list<factored<T, f>> list = get_factored();
    std::vector<factored<T, f>> result(list.size());
    std::copy(list.begin(), list.end(), result.begin());
    return result;
}


/**
 * Specitic types and functions defined to fit our problem
 */
typedef struct quantum_numbers {
    int m_a, n_a, nz_a, m_b, n_b, nz_b, count;
} quantum_numbers;

typedef struct m_n_pair {
    int m_a, n_a;
} m_n_pair;


static inline bool symmetry_filter(quantum_numbers& entry)
{
    if (entry.n_b<entry.n_a && entry.nz_b<entry.nz_a) {
        entry.count *= 2;
        return true;
    }
    else if (entry.n_b<=entry.n_a || entry.nz_b<=entry.nz_a) {
        return true;
    }
    else {
        return false;
    }
}

static inline int select_nza(const quantum_numbers& entry) { return entry.nz_a; }
// or we could use [](const quantum_numbers & q){return q.nz_a;}


static inline int select_nzb(const quantum_numbers& entry) { return entry.nz_b; }

static inline struct m_n_pair select_ma_na(const quantum_numbers& entry) { return {entry.m_a, entry.n_a}; }

static inline bool operator==(const struct m_n_pair l, const struct m_n_pair r) { return l.m_a==r.m_a && l.n_a==r.n_a; }

#endif //PROJET_IPS1_FACTORISATIONFINDER_H