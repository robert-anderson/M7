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
    typedef std::function<uint_t(const uint_t &, const uint_t &)> direct_prod_fn_t;
    typedef const direct_prod_fn_t &cr_direct_prod_fn_t;
    const std::vector<std::string> m_labels;
    dense::Matrix<uint_t> m_products;

    uint_t nirrep() const;

    const uint_t &product(const uint_t &iirrep, const uint_t &jirrep) const;

    const std::string &product_label(const uint_t &iirrep, const uint_t &jirrep) const;

    template<uint_t nirrep>
    bool is_conservative(const uinta_t<nirrep> &iirreps) const {
        static_assert(nirrep > 1, "function applies to direct products of irreps");
        uint_t tmp = iirreps[0];
        for (uint_t i = 1ul; i < nirrep - 1; ++i) tmp = product(tmp, iirreps[i]);
        return tmp == iirreps.back();
    }

    template<uint_t nirrep>
    bool is_conservative(const uinta_t<nirrep> &iirreps1, const uinta_t<nirrep> &iirreps2) const {
        static_assert(nirrep > 0, "function applies to direct products of irreps");
        uint_t tmp1 = iirreps1[0];
        uint_t tmp2 = iirreps2[0];
        for (uint_t i = 1ul; i < nirrep; ++i) {
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
    const uintv_t m_site_irreps;
    const uint_t m_nsite;

    AbelianGroupMap(AbelianGroup grp, uintv_t site_irreps);

    AbelianGroupMap(uint_t nsite);

    bool operator==(const AbelianGroupMap& other) const;

    operator bool() const;
    
};

#endif //M7_ABELIANGROUP_H
