//
// Created by rja on 09/02/2021.
//

#ifndef M7_TABLE_H
#define M7_TABLE_H

#include <stack>
#include <src/core/field/RowHdf5.h>
#include <src/core/field/Fields.h>
#include <src/core/sort/ExtremalIndices.h>
#include "src/core/field/Row.h"
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

    std::string to_string(const defs::inds *ordering = nullptr) const override {
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

    virtual void write_rows(RowHdf5Writer<row_t>& row_writer) const {
        size_t iitem = 0ul;
        log::debug_("beginning HDF5 write loop over rows");
        for (row_writer.restart(); row_writer.in_range(); row_writer.step()) {
            if (!row_writer.is_h5_write_exempt()) row_writer.write(iitem++);
        }
        while (iitem<row_writer.m_nitem_max) row_writer.write(iitem++);
        mpi::barrier();
        log::debug_("ending HDF5 write loop over rows");
    }

    virtual void read_rows(RowHdf5Reader<row_t>& row_reader) {
        size_t iitem = 0ul;
        TableBase::clear();
        log::debug_("beginning HDF5 read loop over rows");
        push_back(row_reader.m_nitem);
        for (row_reader.restart(); row_reader.in_range(); row_reader.step()){
            row_reader.read(iitem++);
        }
        while (iitem<row_reader.m_nitem_max) row_reader.read(iitem++);
        mpi::barrier();
        log::debug_("ending HDF5 read loop over rows");
    }

public:

    virtual void save(hdf5::GroupWriter &parent, std::string name, std::vector<std::string> field_names) const {
        RowHdf5Writer<row_t> row_writer(m_row, parent, name, nrow_to_write(), field_names);
        write_rows(row_writer);
    }

    virtual void save(hdf5::GroupWriter &parent, std::string name) const {
        RowHdf5Writer<row_t> row_writer(m_row, parent, name, nrow_to_write());
        write_rows(row_writer);
    }

    virtual void load(hdf5::GroupReader &parent, std::string name) {
        RowHdf5Reader<row_t> row_reader(m_row, parent, name);
        read_rows(row_reader);
    }

private:
    Row &base_row() {
        return static_cast<Row &>(m_row);
    }
};


#endif //M7_TABLE_H
