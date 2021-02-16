//
// Created by rja on 09/02/2021.
//

#ifndef M7_ROWZ_H
#define M7_ROWZ_H

#include "src/core/util/utils.h"
#include "src/core/parallel/MPIAssert.h"
#include "src/core/table/Buffer.h"

struct FieldBaseZ;

struct RowZ {
    Buffer::Window *m_table_bw = nullptr;
    size_t *m_table_hwm = nullptr;
    mutable defs::data_t *m_dbegin = nullptr;
    mutable size_t m_i = 0ul;
    std::vector<FieldBaseZ *> m_fields;
    size_t m_size;
    size_t m_dsize;
    size_t m_current_offset = 0ul;
    mutable RowZ* m_child = nullptr;

//    void set_table(TableBaseX *table) {
//        m_table = table;
//        select_first();
//    }

    bool in_range() const {
        return m_i < *m_table_hwm;
    }

    /*
     * the 3 "cursor" methods
     */
    void restart() const {
        MPI_ASSERT(m_table_bw, "Row must be assigned to a Table");
        MPI_ASSERT(m_table_hwm, "Row must be assigned to a Table");
        MPI_ASSERT(m_table_bw->m_dbegin, "Row is assigned to Table buffer window without a beginning");
        MPI_ASSERT(m_table_bw->m_dend, "Row is assigned to Table buffer window without an end");
        MPI_ASSERT(in_range(), "Row is out of table bounds");
        m_dbegin = m_table_bw->m_dbegin;
    }

    void step() const {
        MPI_ASSERT(m_table_bw, "Row must be assigned to a Table");
        MPI_ASSERT(m_table_hwm, "Row must be assigned to a Table");
        MPI_ASSERT(m_table_bw->m_dbegin, "Row is assigned to Table buffer window without a beginning");
        MPI_ASSERT(m_table_bw->m_dend, "Row is assigned to Table buffer window without an end");
        MPI_ASSERT(in_range(), "Row is out of table bounds");
        m_dbegin += m_dsize;
        m_i++;
    }

    void jump(const size_t& i) const {
        MPI_ASSERT(m_table_bw, "Row must be assigned to a Table");
        MPI_ASSERT(m_table_hwm, "Row must be assigned to a Table");
        MPI_ASSERT(m_table_bw->m_dbegin, "Row is assigned to Table buffer window without a beginning");
        MPI_ASSERT(m_table_bw->m_dend, "Row is assigned to Table buffer window without an end");
        MPI_ASSERT(in_range(), "Row is out of table bounds");
        m_dbegin = m_table_bw->m_dbegin+m_dsize*i;
        m_i = i;
    }


//    defs::inds field_format() const {
//        defs::inds tmp;
//        tmp.reserve(m_fields.size());
//        for (auto field: m_fields) tmp.push_back(field->m_size * (size_t) &field->m_type_info);
//        return tmp;
//    }

    RowZ() {}

//    RowX(TableBaseX *table) {
//        set_table(table);
//    }
//
    RowZ(const RowZ &other) {
        m_table_bw = other.m_table_bw;
        m_table_hwm = other.m_table_hwm;
        m_dbegin = other.m_dbegin;
        //m_size = other.m_size;
        //m_dsize = other.m_dsize;
        //m_current_offset = other.m_current_offset;
        other.m_child = this;
        ASSERT(m_fields.empty())
    }

    std::string to_string() const;

    size_t add_field(FieldBaseZ *field);

//    void copy_keys(const RowX &other) {
//        ASSERT(other.m_key_field_inds.size() == m_key_field_inds.size());
//        for (const auto &ifield : m_key_field_inds) {
//            m_fields[ifield]->copy(*other.m_fields[ifield]);
//        }
//    }
//
//    defs::hash_t hash_keys() const {
//        defs::hash_t tmp = 0ul;
//        for (const auto &ifield: m_key_field_inds) {
//            const auto field = m_fields[ifield];
//            tmp ^= hashing::fnv_hash(field->begin(), field->m_size);
//        }
//        return tmp;
//    }
//
//    bool equal_keys(const RowX &other) const {
//        ASSERT(other.m_key_field_inds.size() == m_key_field_inds.size());
//        for (const auto &ifield : m_key_field_inds) {
//            if (!m_fields[ifield]->equals(*other.m_fields[ifield])) return false;
//        }
//        return true;
//    }

};


#endif //M7_ROWZ_H
