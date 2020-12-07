/**
 * @file FactorisationHelper.hpp
 */

#ifndef PROJET_IPS1_FACTORISATIONFINDER_H
#define PROJET_IPS1_FACTORISATIONFINDER_H


#ifndef likely
#define likely(x) __builtin_expect((x), 1)
#endif


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
  inline factored(fa factor, std::vector<T>&& vec)
          :factor(factor), factored_out(std::move(vec))
  {
      factored_out.reserve(50);
  }
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
     * @param input vector of type T entries
     * @param filter to filter some values out
     * @param selector to select the values we want to factor out;
     */
    FactorisationHelper(std::vector<T>& input, selector_function selector, input_filter filter = [](T& a) { return true; });

    /**
     * Adds an entry to the factoriser
     * @param entry
     */
    inline void add(T&& entry);
    // inline void add(T entry);

    /**
     * Called at the end to get the factorisation result
     * @return a list of structs. Each struct contains a "factor" which is the key that was factored out
     * of the values that are in the list "factored_out"
     */
    //std::vector<struct factored<T, f>> get_factored();

    inline factored<T, f>& operator[](int index) { return out[index]; };

    inline typename std::vector<struct factored<T, f>>::iterator begin() { return out.begin(); };

    inline typename std::vector<struct factored<T, f>>::iterator end() { return out.end(); };

    inline size_t size() { return out.size(); };

private:
    const input_filter filter;
    const selector_function selector;
    std::vector<struct factored<T, f>> out{};

    /**
     * As we add entries, they are placed in the right categories
     * to factor them out given a selector that "selects" a peculiar value in our entry
     * @tparam fa the type of what we factor out (can be a pair of values, of whatever)
     * @tparam T the type of the entries we added
     * @param entry
     */
    inline void dispatch_entry(const T& entry);
};

/**
 * As we add entries, they are placed in the right categories
 * to factor them out given a selector that "selects" a peculiar value in our entry
 * @tparam fa the type of what we factor out (can be a pair of values, of whatever)
 * @tparam T the type of the entries we added
 * @param entry
 */
template<typename T, typename f>
inline void FactorisationHelper<T, f>::dispatch_entry(const T& entry)
{
    f fac = selector(entry);
    auto it = std::find_if(out.begin(), out.end(), [&fac](const struct factored<T, f>& x) { return x.factor==fac; });
    if (it!=out.end()) {
        it->factored_out.push_back(std::move(entry));
    }
    else {
        out.push_back({fac, std::vector<T>{std::move(entry)}});
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
FactorisationHelper<T, f>::FactorisationHelper(selector_function selec, input_filter filt)
        :filter(filt), selector(selec)
{
    out.reserve(500);
}

template<typename T, typename f>
FactorisationHelper<T, f>::FactorisationHelper(std::vector<T>& input, selector_function select, input_filter filt)
        : FactorisationHelper<T, f>::FactorisationHelper(select, filt)
{
    out.reserve(10+input.size()/5);
    for (auto& in : input) {
        if (likely(filter(in))) {
            dispatch_entry(std::move(in));
        }
    }
}

template<typename T, typename f>
inline void FactorisationHelper<T, f>::add(T&& entry)
{
    if (likely(filter(entry))) {
        dispatch_entry(entry);
    }
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

static inline struct m_n_pair select_mana(const quantum_numbers& entry) { return {entry.m_a, entry.n_a}; }

static inline bool operator==(const struct m_n_pair l, const struct m_n_pair r) { return l.m_a==r.m_a && l.n_a==r.n_a; }

#endif //PROJET_IPS1_FACTORISATIONFINDER_H