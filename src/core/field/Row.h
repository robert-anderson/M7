//
// Created by rja on 09/02/2021.
//

#ifndef M7_ROW_H
#define M7_ROW_H

#include "src/core/util/utils.h"
#include "src/core/io/HDF5Wrapper.h"
#include "src/core/parallel/MPIAssert.h"
#include "src/core/table/Buffer.h"

struct FieldBase;

struct Row {
    Buffer::Window *m_table_bw = nullptr;
    size_t *m_table_hwm = nullptr;
    mutable defs::data_t *m_dbegin = nullptr;
    mutable size_t m_i = 0ul;
    std::vector<FieldBase *> m_fields;
    size_t m_size = 0ul;
    size_t m_dsize = 0ul;
    size_t m_current_offset = 0ul;
    mutable Row* m_child = nullptr;

    bool in_range() const {
        return m_i < *m_table_hwm;
    }

    bool ptr_in_range() const {
        return (m_dbegin >= m_table_bw->m_dbegin) && (m_dbegin < m_table_bw->m_dbegin+*m_table_hwm*m_dsize);
    }

    defs::data_t* dbegin() {
        MPI_ASSERT(m_dbegin, "the row pointer is not set")
        MPI_ASSERT(ptr_in_range(), "the row is not pointing memory in the permitted range");
        return m_dbegin;
    }

    const defs::data_t* dbegin() const {
        return m_dbegin;
    }

    /*
     * the 3 "cursor" methods
     */
    void restart() const {
        MPI_ASSERT(m_table_bw, "Row must be assigned to a Table");
        MPI_ASSERT(m_table_hwm, "Row must be assigned to a Table");
        MPI_ASSERT(m_table_bw->m_dbegin, "Row is assigned to Table buffer window without a beginning");
        MPI_ASSERT(m_table_bw->m_dend, "Row is assigned to Table buffer window without an end");
        m_dbegin = m_table_bw->m_dbegin;
        m_i = 0ul;
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
        m_dbegin = m_table_bw->m_dbegin+m_dsize*i;
        m_i = i;
        MPI_ASSERT(in_range(), "Row is out of table bounds");
    }

    Row() {}

    Row(const Row &other) {
        m_table_bw = other.m_table_bw;
        m_table_hwm = other.m_table_hwm;
        m_dbegin = other.m_dbegin;
        other.m_child = this;
        ASSERT(m_fields.empty())
    }

    std::string to_string() const;

    size_t add_field(FieldBase *field);

    void clear();

    bool is_cleared() const;

};

template<typename row_t>
struct KeyField {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
    typedef typename std::remove_reference<typename std::result_of<decltype(&row_t::key_field)(row_t)>::type>::type type;
    static type& get(row_t& row){return row.key_field();}
    static const type& get(const row_t& row) {return const_cast<row_t&>(row).key_field();}
};


#endif //M7_ROW_H
