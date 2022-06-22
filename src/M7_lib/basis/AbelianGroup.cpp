//
// Created by Robert J. Anderson on 07/06/2021.
//

#include "AbelianGroup.h"

size_t AbelianGroup::nirrep() const {
    return m_labels.size();
}

const size_t &AbelianGroup::product(const size_t &iirrep, const size_t &jirrep) const {
    return m_products(iirrep, jirrep);
}

const std::string &AbelianGroup::product_label(const size_t &iirrep, const size_t &jirrep) const {
    return m_labels[product(iirrep, jirrep)];
}

AbelianGroup::AbelianGroup(std::vector<std::string> labels, AbelianGroup::direct_prod_fn_t direct_prod_fn) :
        m_labels(std::move(labels)), m_products({nirrep(), nirrep()}) {
    for (size_t iirrep = 0ul; iirrep < nirrep(); ++iirrep)
        for (size_t jirrep = 0ul; jirrep < nirrep(); ++jirrep)
            m_products(iirrep, jirrep) = direct_prod_fn(iirrep, jirrep);
}

AbelianGroup::AbelianGroup() : AbelianGroup({"A"}, [](const size_t&, const size_t&){return 0ul;}){}

AbelianGroup::AbelianGroup(const AbelianGroup &other) :
        m_labels(other.m_labels), m_products(other.m_products){
}

AbelianGroup::AbelianGroup(AbelianGroup &&other) :
        m_labels(other.m_labels), m_products(other.m_products){}

std::vector<std::string> AbelianGroup::merge_labels(const AbelianGroup &g1, const AbelianGroup &g2) {
    std::vector<std::string> tmp;
    tmp.reserve(g1.nirrep() * g2.nirrep());
    for (const auto &l1: g1.m_labels)
        for (const auto &l2: g2.m_labels)
            tmp.emplace_back(l1 + " " + l2);
    return tmp;
}

AbelianGroup::AbelianGroup(const AbelianGroup &g1, const AbelianGroup &g2) :
        AbelianGroup(merge_labels(g1, g2),
                     [&](const size_t &iirrep, const size_t &jirrep) {
                         // combine the product tables of the two groups
                         const size_t iirrep1 = iirrep / g2.nirrep();
                         const size_t iirrep2 = iirrep - iirrep1 * g2.nirrep();
                         const size_t jirrep1 = jirrep / g2.nirrep();
                         const size_t jirrep2 = jirrep - jirrep1 * g2.nirrep();
                         return g1.m_products(iirrep1, jirrep1) * g2.nirrep() + g2.m_products(iirrep2, jirrep2);
                     }) {}

bool AbelianGroup::operator==(const AbelianGroup &other) const {
    return m_labels==other.m_labels && m_products==other.m_products;
}

PointGroup::PointGroup() : AbelianGroup({"0", "1", "2", "3", "4", "5", "6", "7"},
                                        [](const size_t& i, const size_t& j){return i^j;}){}

AbelianGroupMap::AbelianGroupMap(AbelianGroup grp, defs::inds_t site_irreps) :
        m_grp(grp), m_site_irreps(site_irreps), m_nsite(m_site_irreps.size()) {
    REQUIRE_TRUE(std::all_of(m_site_irreps.cbegin(), m_site_irreps.cend(),
                       [&](size_t i){return i<m_grp.nirrep();}), "irrep label OOB");
}

AbelianGroupMap::AbelianGroupMap(size_t nsite) : AbelianGroupMap({}, defs::inds_t(nsite, 0)){}

bool AbelianGroupMap::operator==(const AbelianGroupMap &other) const {
    return (m_grp == other.m_grp) && (m_site_irreps == other.m_site_irreps) && (m_nsite == other.m_nsite);
}

AbelianGroupMap::operator bool() const {
    return !m_site_irreps.empty();
}
