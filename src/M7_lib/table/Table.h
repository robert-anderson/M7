//
// Created by Robert J. Anderson on 09/02/2021.
//

#ifndef M7_TABLE_H
#define M7_TABLE_H

#include <stack>

#include <M7_lib/field/RowHdf5.h>
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

    Table(const Table<row_t> &other) : Table(other.m_row) {
        ASSERT(static_cast<Row &>(m_row).m_table == this);
    }

    virtual ~Table(){}

    std::string to_string(const defs::ivec_t *ordering = nullptr) const override {
        std::string tmp = m_row.field_names_string();
        if (!m_hwm) return tmp;
        const auto n = ordering ? std::min(ordering->size(), m_hwm) : m_hwm;
        auto row = m_row;
        for (size_t iirow = 0ul; iirow < n; ++iirow) {
            auto irow = ordering ? (*ordering)[iirow] : iirow;
            row.jump(irow);
            tmp += std::to_string(irow) + ". " + row.to_string() + "\n";
        }
        DEBUG_ASSERT_TRUE(row.in_range(), "row is out of range");
        return tmp;
    }

private:
    size_t nrow_to_write() const {
        auto row = m_row;
        size_t n = 0ul;
        for (row.restart(); row.in_range(); row.step()) {
            n += !static_cast<const Row&>(row).is_h5_write_exempt();
        }
        return n;
    }


public:

    virtual void save(hdf5::GroupWriter &parent, std::string name, std::vector<std::string> field_names) const {
        RowHdf5Writer<row_t>(m_row, parent, name, nrow_to_write(), field_names).write();
    }

    virtual void save(hdf5::GroupWriter &parent, std::string name) const {
        RowHdf5Writer<row_t>(m_row, parent, name, nrow_to_write()).write();
    }

    virtual void load(hdf5::GroupReader &parent, std::string name) {
        RowHdf5Reader<row_t>(m_row, parent, name).read();
    }

private:
    Row &base_row() {
        return static_cast<Row &>(m_row);
    }
};


#endif //M7_TABLE_H
