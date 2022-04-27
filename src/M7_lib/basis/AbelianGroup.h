//
// Created by Robert J. Anderson on 07/06/2021.
//

#ifndef M7_ABELIANGROUP_H
#define M7_ABELIANGROUP_H

#include <functional>
#include <utility>
#include <vector>

#include <M7_lib/linalg/Dense.h>

struct AbelianGroup {
    typedef std::function<size_t(const size_t &, const size_t &)> direct_prod_fn_t;
    typedef const direct_prod_fn_t &cr_direct_prod_fn_t;
    const std::vector<std::string> m_labels;
    dense::Matrix<size_t> m_products;

    size_t nirrep() const;

    const size_t &product(const size_t &iirrep, const size_t &jirrep) const;

    const std::string &product_label(const size_t &iirrep, const size_t &jirrep) const;

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

    AbelianGroup(std::vector<std::string> labels, direct_prod_fn_t direct_prod_fn);

    /**
     * default group has only one irrep label
     */
    AbelianGroup();

    AbelianGroup(const AbelianGroup& other);

    AbelianGroup(AbelianGroup&& other);

    static std::vector<std::string> merge_labels(const AbelianGroup &g1, const AbelianGroup &g2);

    AbelianGroup(const AbelianGroup &g1, const AbelianGroup &g2);;

    bool operator==(const AbelianGroup& other) const;
};

/**
 * FCIDUMPs don't specify exactly which symmetry labels are identified by each index. All supported point groups have
 * the same XOR multiplication table so knowing the exact point group is only necessary if we want to accurately report
 * the symmetries of orbitals.
 */
struct PointGroup : AbelianGroup {
    PointGroup();
};
// D2h; "A1g", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"

struct AbelianGroupMap {
    const AbelianGroup m_grp;
    const defs::inds m_site_irreps;
    const size_t m_nsite;

    AbelianGroupMap(AbelianGroup grp, defs::inds site_irreps);

    AbelianGroupMap(size_t nsite);

    bool operator==(const AbelianGroupMap& other) const;

    operator bool() const;
};

#endif //M7_ABELIANGROUP_H
