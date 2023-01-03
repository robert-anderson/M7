//
// Created by Robert J. Anderson on 09/02/2021.
//

#ifndef M7_TABLE_H
#define M7_TABLE_H

#include <stack>

#include <M7_lib/field/Fields.h>
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
        if (m_bw.m_hwm == m_bw.m_begin) return tmp;
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

private:
    uint_t nrow_to_write() const {
        uint_t n = 0ul;
        for (m_row.restart(); m_row; ++m_row) {
            n += !static_cast<const Row&>(m_row).is_h5_write_exempt();
        }
        return n;
    }

public:
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
     */
    virtual void save(const hdf5::NodeWriter& parent, str_t name, strm_t field_names, bool this_rank) const {
        hdf5::GroupWriter gw(parent, name);
        for (auto field : m_row.m_fields) {
            auto it = std::find_if(field_names.cbegin(), field_names.cend(),
                       [&field](const strp_t& pair){return field->m_name==pair.first;});
            if (it != field_names.cend()) field->save(gw, it->second, this_rank);
        }
    }
    /**
     * overload in the case that the FieldBase::m_name is identical to the intended HDF5 dataset names
     */
    virtual void save(const hdf5::NodeWriter& parent, str_t name, strv_t field_names, bool this_rank) const {
        strm_t pairs;
        for (const auto& field_name : field_names) pairs.insert({field_name, field_name});
        save(parent, name, pairs, this_rank);
    }
    /**
     * overload in the case that all fields are to be saved with their FieldBase::m_name as dataset name
     */
    void save(const hdf5::NodeWriter& parent, str_t name, bool this_rank) const {
        save(parent, name, m_row.all_field_names(), this_rank);
    }

    /**
     * @param parent
     *  HDF5 node from which to load this Table's fields as datasets
     * @param name
     *  name to be given to the group in the parent node
     * @param field_names
     *  pairs of field names, first in pair is FieldBase::m_name of field selected for saving, second is intended
     *  dataset name
     * @param this_rank
     *  true if this MPI rank participates in writing
     */
    virtual void load(const hdf5::NodeReader& parent, str_t name, strm_t field_names, bool part, bool this_rank) {
        hdf5::GroupReader gr(parent, name);
        for (auto field : m_row.m_fields) {
            auto it = std::find_if(field_names.cbegin(), field_names.cend(),
                                   [&field](const strp_t& pair){return field->m_name==pair.first;});
            if (it != field_names.cend()) field->load(gr, it->second, part, this_rank);
        }
    }
    /**
     * overload in the case that the FieldBase::m_name is identical to the intended HDF5 dataset names
     */
    virtual void load(const hdf5::NodeReader& parent, str_t name, strv_t field_names, bool part, bool this_rank) {
        strm_t pairs;
        for (const auto& field_name : field_names) pairs.insert({field_name, field_name});
        load(parent, name, pairs, part, this_rank);
    }
    /**
     * overload in the case that all fields are to be loaded with their FieldBase::m_name as dataset name
     */
    void load(const hdf5::NodeReader& parent, str_t name, bool part, bool this_rank) {
        load(parent, name, m_row.all_field_names(), part, this_rank);
    }

};


#endif //M7_TABLE_H
