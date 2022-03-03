//
// Created by rja on 07/06/2021.
//

#ifndef M7_ABELIANGROUP_H
#define M7_ABELIANGROUP_H

#include <functional>
#include <utility>
#include <src/core/linalg/Dense.h>
#include <vector>

struct AbelianGroup {
    typedef std::function<size_t(const size_t &, const size_t &)> direct_prod_fn_t;
    typedef const direct_prod_fn_t &cr_direct_prod_fn_t;
    const std::vector<std::string> m_labels;
    dense::Matrix<size_t> m_products;

    size_t nirrep() const {
        return m_labels.size();
    }

    const size_t &product(const size_t &iirrep, const size_t &jirrep) const {
        return m_products(iirrep, jirrep);
    }

    const std::string &product_label(const size_t &iirrep, const size_t &jirrep) const {
        return m_labels[product(iirrep, jirrep)];
    }

    template<size_t nirrep>
    bool is_conservative(const std::array<size_t, nirrep> &iirreps) const {
        static_assert(nirrep > 1, "function applies to direct products of irreps");
        size_t tmp = iirreps[0];
        for (size_t i = 1ul; i < nirrep - 1; ++i) tmp = product(tmp, iirreps[i]);
        return tmp == iirreps.back();
    }

    template<size_t nirrep>
    bool is_conservative(const std::array<size_t, nirrep> &iirreps1, const std::array<size_t, nirrep> &iirreps2) const {
        static_assert(nirrep > 0, "function applies to direct products of irreps");
        size_t tmp1 = iirreps1[0];
        size_t tmp2 = iirreps2[0];
        for (size_t i = 1ul; i < nirrep; ++i) {
            tmp1 = product(tmp1, iirreps1[i]);
            tmp2 = product(tmp2, iirreps2[i]);
        }
        return tmp1 == tmp2;
    }

    AbelianGroup(std::vector<std::string> labels, direct_prod_fn_t direct_prod_fn) :
            m_labels(std::move(labels)), m_products({nirrep(), nirrep()}) {
        for (size_t iirrep = 0ul; iirrep < nirrep(); ++iirrep)
            for (size_t jirrep = 0ul; jirrep < nirrep(); ++jirrep)
                m_products(iirrep, jirrep) = direct_prod_fn(iirrep, jirrep);
    }

    /**
     * default group has only one irrep label
     */
    AbelianGroup(): AbelianGroup({"A"}, [](const size_t&, const size_t&){return 0ul;}){}

    AbelianGroup(const AbelianGroup& other):
            m_labels(other.m_labels), m_products(other.m_products){
    }

    AbelianGroup(AbelianGroup&& other):
            m_labels(other.m_labels), m_products(other.m_products){}

    static std::vector<std::string> merge_labels(const AbelianGroup &g1, const AbelianGroup &g2) {
        std::vector<std::string> tmp;
        tmp.reserve(g1.nirrep() * g2.nirrep());
        for (const auto &l1: g1.m_labels)
            for (const auto &l2: g2.m_labels)
                tmp.emplace_back(l1 + " " + l2);
        return tmp;
    }

    AbelianGroup(const AbelianGroup &g1, const AbelianGroup &g2) :
            AbelianGroup(merge_labels(g1, g2),
             [&](const size_t &iirrep, const size_t &jirrep) {
                 // combine the product tables of the two groups
                 const size_t iirrep1 = iirrep / g2.nirrep();
                 const size_t iirrep2 = iirrep - iirrep1 * g2.nirrep();
                 const size_t jirrep1 = jirrep / g2.nirrep();
                 const size_t jirrep2 = jirrep - jirrep1 * g2.nirrep();
                 return g1.m_products(iirrep1, jirrep1) * g2.nirrep() + g2.m_products(iirrep2, jirrep2);
             }) {};

};

/**
 * FCIDUMPs don't specify exactly which symmetry labels are identified by each index. All supported point groups have
 * the same XOR multiplication table so knowing the exact point group is only necessary if we want to accurately report
 * the symmetries of orbitals.
 */
struct PointGroup : AbelianGroup {
    PointGroup() : AbelianGroup({"0", "1", "2", "3", "4", "5", "6", "7"},
               [](const size_t& i, const size_t& j){return i^j;}){}
};
// D2h; "A1g", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"

struct AbelianGroupMap {
    const AbelianGroup m_grp;
    const defs::inds m_site_irreps;
    const size_t m_nsite;

    AbelianGroupMap(AbelianGroup grp, defs::inds site_irreps) :
            m_grp(grp), m_site_irreps(site_irreps), m_nsite(m_site_irreps.size()) {
        ASSERT(std::all_of(m_site_irreps.cbegin(), m_site_irreps.cend(),
                           [&](const size_t& i){return i<m_grp.nirrep();}));
    }

    AbelianGroupMap(size_t nsite): AbelianGroupMap({}, defs::inds(nsite, 0)){}
};

#endif //M7_ABELIANGROUP_H
