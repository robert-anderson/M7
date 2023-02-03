//
// Created by Robert J. Anderson on 09/02/2021.
//

#ifndef M7_TABLE_H
#define M7_TABLE_H

#include <stack>

#include <M7_lib/sort/ExtremalIndices.h>
#include <M7_lib/field/Row.h>

#include "TableBase.h"

template<typename row_t>
struct Table : TableBase {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
    row_t m_row;

    Table(const row_t &row) :
            TableBase(static_cast<const Row &>(row).m_size), m_row(row) {
        static_cast<Row &>(m_row).m_table = this;
        static_cast<Row &>(m_row).select_null();
    }

    Table& operator=(const Table& other) {
        TableBase::operator=(other);
        return *this;
    }

    Table(const Table<row_t> &other) : Table(other.m_row) {
        ASSERT(static_cast<Row &>(m_row).m_table == this);
        // only set the layout, not the buffer, since the result is not yet associated with a buffer
    }

    virtual ~Table(){}

    str_t to_string(const uintv_t *ordering = nullptr) const override {
        str_t tmp = m_row.field_names_string();
        if (m_bw.cend() == m_bw.cbegin()) return tmp;
        const auto n = ordering ? std::min(ordering->size(), nrow_in_use()) : nrow_in_use();
        for (uint_t iirow = 0ul; iirow < n; ++iirow) {
            auto irow = ordering ? (*ordering)[iirow] : iirow;
            m_row.jump(irow);
            tmp += std::to_string(irow) + ". " + m_row.to_string() + "\n";
        }
        DEBUG_ASSERT_TRUE(bool(m_row), "row is not dereferencable");
        return tmp;
    }

    /**
     * set the table pointer of the given row to this table
     */
    void associate(row_t& row) {
        static_cast<Row&>(row).m_table = this;
    }

    bool associated(const row_t& row) const {
        return static_cast<const Row&>(row).m_table == this;
    }

protected:
    /**
     * @param parent
     *  HDF5 node within which to save this Table's fields as datasets
     * @param name
     *  name to be given to the group in the parent node
     * @param field_names
     *  pairs of field names, first in pair is FieldBase::m_name of field selected for saving, second is intended
     *  dataset name
     * @param this_rank
     *  true if this MPI rank participates in writing
     * @param max_nitem_per_op
     *  maximum number of items (in-use rows) saved in a single write operation
     */
    virtual void save_fn(const hdf5::NodeWriter& parent, const str_t& name,
                         strm_t field_names, bool this_rank, uint_t max_nitem_per_op) const {
        hdf5::GroupWriter gw(parent, name);
        for (auto field : m_row.m_fields) {
            auto it = std::find_if(field_names.cbegin(), field_names.cend(),
                       [&field](const strp_t& pair){return field->m_name==pair.first;});
            if (it != field_names.cend()) field->save(gw, it->second, this_rank, max_nitem_per_op);
        }
    }

    /**
     * @param parent
     *  HDF5 node from which to load this Table's fields as datasets
     * @param name
     *  name to be given to the group in the parent node
     * @param max_nitem_per_op
     *  maximum number of items (in-use rows) loaded in a single read operation
     * @param field_names
     *  pairs of field names, first in pair is FieldBase::m_name of field selected for saving, second is intended
     *  dataset name
     * @param part
     *  true each MPI rank only loads part of the dataset
     * @param this_rank
     *  true if this MPI rank participates in writing
     */
    virtual void load_fn(const hdf5::NodeReader& parent, const str_t& name, uint_t max_nitem_per_op,
                         strm_t field_names, bool part, bool this_rank) {
        clear();
        hdf5::GroupReader gr(parent, name);
        for (auto field : m_row.m_fields) {
            auto it = std::find_if(field_names.cbegin(), field_names.cend(),
                                   [&field](const strp_t& pair){return field->m_name==pair.first;});
            if (it != field_names.cend()) field->load(gr, it->second, max_nitem_per_op, part, this_rank);
        }
    }


public:
    void save(const hdf5::NodeWriter& parent, const str_t& name,
              strm_t field_names, bool this_rank, uint_t max_nitem_per_op) const {
        save_fn(parent, name, field_names, this_rank, max_nitem_per_op);
    }

    /**
     * overload in the case that the FieldBase::m_name is identical to the intended HDF5 dataset names
     */
    void save(const hdf5::NodeWriter& parent, const str_t& name,
              const strv_t &field_names, bool this_rank, uint_t max_nitem_per_op) const {
        strm_t pairs;
        for (const auto& field_name : field_names) pairs.insert({field_name, field_name});
        save(parent, name, pairs, this_rank, max_nitem_per_op);
    }

    /**
     * overload in the case that all fields are to be saved with their FieldBase::m_name as dataset name
     */
    void save(const hdf5::NodeWriter& parent, const str_t& name, bool this_rank, uint_t max_nitem_per_op) const {
        save(parent, name, m_row.all_field_names(), this_rank, max_nitem_per_op);
    }

    void load(const hdf5::NodeReader& parent, const str_t& name, uint_t max_nitem_per_op,
              strm_t field_names, bool part, bool this_rank) {
        load_fn(parent, name, max_nitem_per_op, field_names, part, this_rank);
    }

    /**
     * overload in the case that the FieldBase::m_name is identical to the intended HDF5 dataset names
     */
    void load(const hdf5::NodeReader& parent, const str_t& name, uint_t max_nitem_per_op,
              const strv_t &field_names, bool part, bool this_rank) {
        strm_t pairs;
        for (const auto& field_name : field_names) pairs.insert({field_name, field_name});
        load(parent, name, max_nitem_per_op, pairs, part, this_rank);
    }

    /**
     * overload in the case that all fields are to be loaded with their FieldBase::m_name as dataset name
     */
    void load(const hdf5::NodeReader& parent, const str_t& name, uint_t max_nitem_per_op, bool part, bool this_rank) {
        load(parent, name, max_nitem_per_op, m_row.all_field_names(), part, this_rank);
    }

    /*
     * the following overloads simply assume the default max_nitem_per_op
     */

    void save(const hdf5::NodeWriter& parent, const str_t& name, strm_t field_names, bool this_rank) const {
        save(parent, name, field_names, this_rank, hdf5::c_default_max_nitem_per_op);
    }

    void save(const hdf5::NodeWriter& parent, const str_t& name, const strv_t &field_names, bool this_rank) const {
        save(parent, name, field_names, this_rank, hdf5::c_default_max_nitem_per_op);
    }

    void save(const hdf5::NodeWriter& parent, const str_t& name, bool this_rank) const {
        save(parent, name, this_rank, hdf5::c_default_max_nitem_per_op);
    }

    void load(const hdf5::NodeReader& parent, const str_t& name, strm_t field_names, bool part, bool this_rank) {
        load(parent, name, hdf5::c_default_max_nitem_per_op, field_names, part, this_rank);
    }

    void load(const hdf5::NodeReader& parent, const str_t& name, const strv_t &field_names, bool part, bool this_rank) {
        load(parent, name, hdf5::c_default_max_nitem_per_op, field_names, part, this_rank);
    }

    void load(const hdf5::NodeReader& parent, const str_t& name, bool part, bool this_rank) {
        load(parent, name, hdf5::c_default_max_nitem_per_op, part, this_rank);
    }
};


#endif //M7_TABLE_H
