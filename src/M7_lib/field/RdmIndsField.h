//
// Created by Robert J. Anderson on 12/08/2021.
//

#ifndef M7_RDMINDSFIELD_H
#define M7_RDMINDSFIELD_H

#include <M7_lib/connection/FrmBosOnvConnection.h>

#include "NumberField.h"

/**
 * accesses an operator-type specific partition of a field containing all indices of the RDM element
 */
class RdmIndsPartition {
    NdNumberField<rdm_ind_t, 1> &m_field;
    const uint_t m_offset, m_size;
public:
    RdmIndsPartition(NdNumberField<rdm_ind_t, 1> &field, uint_t offset, uint_t size) :
        m_field(field), m_offset(offset), m_size(size) {}

    const rdm_ind_t &operator[](const uint_t &iind) const {
        DEBUG_ASSERT_LT(iind, m_size, "operator index OOB in partition");
        return m_field[m_offset + iind];
    };

    rdm_ind_t &operator[](const uint_t &iind) {
        DEBUG_ASSERT_LT(iind, m_size, "operator index OOB in partition");
        return m_field[m_offset + iind];
    };

    bool operator==(const RdmIndsPartition& other) const {
        DEBUG_ASSERT_EQ(m_size, other.m_size, "comparing incompatible partitions");
        return std::memcmp(m_field.ctbegin()+m_offset, other.m_field.ctbegin()+m_offset, m_size)==0;
    }

    const uint_t &size() const {
        return m_size;
    }

    void zero() {
        for (uint_t i=0ul; i<m_size; ++i) (*this)[i]=0ul;
    }

    /**
     * can't simply copy uintv_t since the rdm_ind_t is not the same as uint_t in general, so must do a narrowing loop
     * @param inds
     *  indices in unsigned system words
     * @return
     *  reference to this
     */
    RdmIndsPartition &operator=(const uintv_t &inds) {
        DEBUG_ASSERT_EQ(inds.size(), m_size, "uintv_t vector size incompatible");
        uint_t iind = 0ul;
        for (auto &ind: inds) (*this)[iind++] = ind;
        return *this;
    }
};

/**
 * a pair of partitions, one for creation operators, and one for annihilations
 */
struct RdmIndsPair {
    RdmIndsPartition m_cre;
    RdmIndsPartition m_ann;

    RdmIndsPair(NdNumberField<rdm_ind_t, 1> &field, uint_t cre_offset,
                uint_t cre_size, uint_t ann_offset, uint_t ann_size) :
            m_cre(field, cre_offset, cre_size), m_ann(field, ann_offset, ann_size) {}

    RdmIndsPair &operator=(const std::pair<uintv_t, uintv_t> &inds) {
        m_cre = inds.first;
        m_ann = inds.second;
        return *this;
    }

    void zero() {
        m_cre.zero();
        m_ann.zero();
    }

    void conjugate() {
        DEBUG_ASSERT_EQ(m_ann.size(), m_cre.size(), "can't conjugate unless numbers of cre and ann ops are the same");
        for (uint_t i=0ul; i<m_cre.size(); ++i) std::swap(m_cre[i], m_ann[i]);
    }
};

/**
 * Field for storage of general MAE indices,
 * Constructs two RdmIndsPair objects as views on the fermion and boson sectors of the indices respectively.
 * Unlike the MBF case, in this situation it is overhead-free to handle fermions and bosons at the same time, because
 * they do not demand fundamentally different datastructures. If the calculation is purely fermionic, the m_bos view
 * will simply not be accessed, since neither view is considered in hashing and comparison operations, the underlying
 * field is.
 */
struct RdmIndsField : NdNumberField<rdm_ind_t, 1> {
    typedef NdNumberField<rdm_ind_t, 1> base_t;
    using base_t::operator=;
    const OpSig m_exsig;
    const uinta_t<4> m_nops;
    const uinta_t<4> m_nop_offsets;
    RdmIndsPair m_frm;
    RdmIndsPair m_bos;

private:
    uinta_t<4> make_nops() const;

    uinta_t<4> make_nop_offsets() const;

public:
    RdmIndsField(Row *row, OpSig exsig, str_t name = "indices");

    RdmIndsField(const RdmIndsField& other);

    virtual ~RdmIndsField() = default;

    RdmIndsField &operator=(const RdmIndsField &other);
    RdmIndsField &operator=(const FrmOnvConnection& conn);
    RdmIndsField &operator=(const BosOnvConnection& conn);
    RdmIndsField &operator=(const FrmBosOnvConnection& conn);

private:
    bool is_ordered(const uint_t &iop, bool strict) const;

public:
    bool is_ordered() const;

    void common_frm_inds(uintv_t &common) const;

    str_t get_exsig_string() const;
};


#endif //M7_RDMINDSFIELD_H
