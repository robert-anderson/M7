//
// Created by rja on 12/08/2021.
//

#ifndef M7_MAEINDSFIELD_H
#define M7_MAEINDSFIELD_H

#include <src/core/connection/FrmBosOnvConnection.h>
#include "NumberField.h"


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

    /**
     * can't simply copy inds since the mev_ind_t is not the same as size_t in general, so must do a narrowing loop
     * @param inds
     *  indices in unsigned system words
     * @return
     *  reference to this
     */
    MaeIndsPartition &operator=(const defs::inds &inds) {
        DEBUG_ASSERT_EQ(inds.size(), m_size, "inds vector size incompatible");
        size_t iind = 0ul;
        for (auto &ind: inds) (*this)[iind++] = ind;
        return *this;
    }
};

struct MaeIndsPartitionPair {
    MaeIndsPartition m_cre;
    MaeIndsPartition m_ann;

    MaeIndsPartitionPair(NdNumberField<defs::mev_ind_t, 1> &field, size_t cre_offset,
                         size_t cre_size, size_t ann_offset, size_t ann_size) :
            m_cre(field, cre_offset, cre_size),
            m_ann(field, ann_offset, ann_size) {}

    MaeIndsPartitionPair &operator=(const std::pair<defs::inds, defs::inds> &inds) {
        m_cre = inds.first;
        m_ann = inds.second;
        return *this;
    }
};

struct MaeIndsField : NdNumberField<defs::mev_ind_t, 1> {
    typedef NdNumberField<defs::mev_ind_t, 1> base_t;
    using base_t::operator=;
    const size_t m_exsig;
    const std::array<size_t, 4> m_nops;
    const std::array<size_t, 4> m_nop_offsets;
    MaeIndsPartitionPair m_frm;
    MaeIndsPartitionPair m_bos;

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

    void common_frm_inds(defs::inds &common) const;

    std::string get_exsig_string() const;
};


#endif //M7_MAEINDSFIELD_H
