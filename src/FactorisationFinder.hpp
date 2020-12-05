#ifndef PROJET_IPS1_FACTORISATIONFINDER_H
#define PROJET_IPS1_FACTORISATIONFINDER_H

#include <list>
#include <algorithm>

struct nuclear_sum_entry {
  int m_a, n_a, nz_a, m_b, n_b, nz_b;
};

struct m_n_pair {
  int m_a, n_a;
};

static inline bool nuclear_symetry(nuclear_sum_entry a, nuclear_sum_entry b)
{
    return a.n_a==b.n_b && a.m_a==b.m_b && a.nz_a==b.nz_b && a.n_b==b.n_a && a.m_b==b.m_a && a.nz_b==b.nz_a;
}

static inline bool nuclear_filter(nuclear_sum_entry entry)
{
    return entry.m_a==entry.m_b;
}

static inline int select_nza(nuclear_sum_entry entry)
{
    return entry.nz_a;
}

static inline int select_nzb(nuclear_sum_entry entry)
{
    return entry.nz_b;
}

static inline struct m_n_pair select_ma_na(nuclear_sum_entry entry)
{
    return {entry.m_a, entry.n_a};
}


static inline bool operator==(const struct m_n_pair l, const struct m_n_pair r)
{
    return l.m_a==r.m_a && l.n_a==r.n_a;
}

template<typename T, typename fa>
struct factored {
  fa factor;
  std::list<T> factored_out;
};

template<typename T, typename f>
class FactorisationHelper {
public:
    typedef bool(* input_filter)(T a);
    typedef bool (* symetry_function)(T a, T b);
    typedef f (* selector_function)(T a);
    FactorisationHelper(input_filter filter, selector_function selector);
    FactorisationHelper(std::list<T> input, input_filter filter, selector_function selector);
    inline void add(T entry);
    std::list<struct factored<T, f>> get_factored();
    //  void apply_symetry(symetry_function fun);
private:
    input_filter filter;
    selector_function selector;
    std::list<struct factored<T, f>> out{};
    void dispatch_entry(T entry);
    void remove_duplicates();
    void remove_symetric_elts();
};

template<typename T, typename f>
void FactorisationHelper<T, f>::dispatch_entry(T entry)
{
    f fac = selector(entry);
    auto it = std::find_if(out.begin(), out.end(), [&fac](const struct factored<T, f>& x) { return x.factor==fac; });
    if (it!=out.end()) {
        it->factored_out.push_back(entry);
    }
    else {
        out.push_back({fac, std::list<T>{entry}});
    }
}
template<typename T, typename f>
FactorisationHelper<T, f>::FactorisationHelper(FactorisationHelper::input_filter filt, FactorisationHelper::selector_function selec)
        :filter(filt), selector(selec) { }

template<typename T, typename f>
void FactorisationHelper<T, f>::add(T entry)
{
    if (filter(entry)) {
        dispatch_entry(entry);
    }
}

template<typename T, typename f>
std::list<struct factored<T, f>> FactorisationHelper<T, f>::get_factored()
{
    remove_duplicates();
    remove_symetric_elts();
    return out;
}

template<typename T, typename f>
FactorisationHelper<T, f>::FactorisationHelper(std::list<T> input, input_filter filt, FactorisationHelper::selector_function select)
        : FactorisationHelper<T, f>::FactorisationHelper(filt, select)
{
    for (auto& in : input) {
        add(in);
    }
}

template<typename T, typename f>
void FactorisationHelper<T, f>::remove_symetric_elts()
{
    //TODO ? But can have a huge negative impact of perfs
}

template<typename T, typename f>
void FactorisationHelper<T, f>::remove_duplicates()
{
    //TODO ? Same huge impact :/
}

#endif //PROJET_IPS1_FACTORISATIONFINDER_H