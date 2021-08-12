//
// Created by rja on 11/08/2021.
//

#ifndef M7_MAETABLE_H
#define M7_MAETABLE_H

#include "src/core/field/Fields.h"

struct MaeInds : fields::Numbers<defs::mev_ind_t, 1> {
    typedef fields::Numbers<defs::mev_ind_t, 1> base_t;
    const size_t m_exsig;
    const std::array<size_t, 4> m_nops;
    const std::array<size_t, 4> m_nop_offsets;

private:
    std::array<size_t, 4> make_nops() const;

    std::array<size_t, 4> make_nop_offsets() const;

public:
    MaeInds(Row *row, size_t exsig, std::string name="indices");

    MaeInds& operator=(const std::array<defs::inds, 4> &inds);

    const size_t& nop(const size_t& iop) const;

    const defs::mev_ind_t& get_ref(const size_t& iop, const size_t& iind) const;

    defs::mev_ind_t& get_ref(const size_t& iop, const size_t& iind);

    template<size_t iop>
    const defs::mev_ind_t& get_ref(const size_t& iind) const{
        static_assert(iop<4, "index OOB");
        DEBUG_ASSERT_LT(iind, m_nops[iop], "operator index OOB");
        return (*this)[m_nop_offsets[iop]+iind];
    }

    template<size_t iop>
    defs::mev_ind_t& get_ref(const size_t& iind) {
        static_assert(iop<4, "index OOB");
        DEBUG_ASSERT_LT(iind, m_nops[iop], "operator index OOB");
        return (*this)[m_nop_offsets[iop]+iind];
    }

private:
    bool is_ordered(const size_t& iop, bool strict) const;
public:
    bool is_ordered() const;

    void common_frm_inds(defs::inds &common) const;

    std::string get_exsig_string() const;
};

struct SpecMomInds : MultiField<MaeInds, MaeInds> {
    typedef MultiField<MaeInds, MaeInds> base_t;
    const size_t m_exsig;
    const std::string m_name;
    MaeInds &m_left, m_right;

    SpecMomInds(Row *row, size_t exsig, std::string name = "indices");

    SpecMomInds(const SpecMomInds &other);
};

#endif //M7_MAETABLE_H
