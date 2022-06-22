//
// Created by Robert J. Anderson on 12/08/2021.
//

#ifndef M7_MAEINDSFIELD_H
#define M7_MAEINDSFIELD_H

#include <M7_lib/connection/FrmBosOnvConnection.h>

#include "NumberField.h"

/**
 * accesses an operator-type specific partition of a field containing all indices of the MAE element
 */
class MaeIndsPartition {
    NdNumberField<defs::mev_ind_t, 1> &m_field;
    const size_t m_offset, m_size;
public:
    MaeIndsPartition(NdNumberField<defs::mev_ind_t, 1> &field, size_t offset, size_t size) :
        m_field(field), m_offset(offset), m_size(size) {}

    const defs::mev_ind_t &operator[](const size_t &iind) const {
        DEBUG_ASSERT_LT(iind, m_size, "operator index OOB in partition");
        return m_field[m_offset + iind];
    };

    defs::mev_ind_t &operator[](const size_t &iind) {
        DEBUG_ASSERT_LT(iind, m_size, "operator index OOB in partition");
        return m_field[m_offset + iind];
    };

    const size_t &size() const {
        return m_size;
    }

    void zero() {
        for (size_t i=0ul; i<m_size; ++i) (*this)[i]=0ul;
    }

    /**
     * can't simply copy inds_t since the mev_ind_t is not the same as size_t in general, so must do a narrowing loop
     * @param inds
     *  indices in unsigned system words
     * @return
     *  reference to this
     */
    MaeIndsPartition &operator=(const defs::inds_t &inds) {
        DEBUG_ASSERT_EQ(inds.size(), m_size, "inds_t vector size incompatible");
        size_t iind = 0ul;
        for (auto &ind: inds) (*this)[iind++] = ind;
        return *this;
    }
};

/**
 * a pair of partitions, one for creation operators, and one for annihilations
 */
struct MaeIndsPair {
    MaeIndsPartition m_cre;
    MaeIndsPartition m_ann;

    MaeIndsPair(NdNumberField<defs::mev_ind_t, 1> &field, size_t cre_offset,
                size_t cre_size, size_t ann_offset, size_t ann_size) :
            m_cre(field, cre_offset, cre_size),
            m_ann(field, ann_offset, ann_size) {}

    MaeIndsPair &operator=(const std::pair<defs::inds_t, defs::inds_t> &inds) {
        m_cre = inds.first;
        m_ann = inds.second;
        return *this;
    }

    void zero() {
        m_cre.zero();
        m_ann.zero();
    }
};

/**
 * Field for storage of general MAE indices,
 * Constructs two MaeIndsPair objects as views on the fermion and boson sectors of the indices respectively.
 * Unlike the ONV case, in this situation it is overhead-free to handle fermions and bosons at the same time, because
 * they do not demand fundamentally different datastructures. If the calculation is purely fermionic, the m_bos view
 * will simply not be accessed, since neither view is considered in hashing and comparison operations, the underlying
 * field is.
 */
struct MaeIndsField : NdNumberField<defs::mev_ind_t, 1> {
    typedef NdNumberField<defs::mev_ind_t, 1> base_t;
    using base_t::operator=;
    const size_t m_exsig;
    const std::array<size_t, 4> m_nops;
    const std::array<size_t, 4> m_nop_offsets;
    MaeIndsPair m_frm;
    MaeIndsPair m_bos;

private:
    std::array<size_t, 4> make_nops() const;

    std::array<size_t, 4> make_nop_offsets() const;

public:
    MaeIndsField(Row *row, size_t exsig, std::string name = "indices");

    MaeIndsField(const MaeIndsField& other);

    MaeIndsField &operator=(const MaeIndsField &other);
    MaeIndsField &operator=(const FrmOnvConnection& conn);
    MaeIndsField &operator=(const BosOnvConnection& conn);
    MaeIndsField &operator=(const FrmBosOnvConnection& conn);

private:
    bool is_ordered(const size_t &iop, bool strict) const;

public:
    bool is_ordered() const;

    void common_frm_inds(defs::inds_t &common) const;

    std::string get_exsig_string() const;
};


#endif //M7_MAEINDSFIELD_H
