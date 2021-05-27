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
            TableBase(static_cast<const Row &>(row).m_dsize), m_row(row) {
        static_cast<Row &>(m_row).m_table = this;
        static_cast<Row &>(m_row).m_dbegin = nullptr;
    }

    Table(const Table<row_t> &other) : Table(other.m_row) {
        ASSERT(static_cast<Row &>(m_row).m_table == this);
    }

    virtual ~Table(){}

    std::string to_string(const defs::inds *ordering = nullptr) const override {
        if (!m_hwm) return "";
        std::string tmp;
        const auto n = ordering ? std::min(ordering->size(), m_hwm) : m_hwm;
        auto row = m_row;
        for (size_t iirow = 0ul; iirow < n; ++iirow) {
            auto irow = ordering ? (*ordering)[iirow] : iirow;
            row.jump(irow);
            tmp += std::to_string(irow) + ". " + row.to_string() + "\n";
        }
        ASSERT(m_row.in_range());
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
        while (iitem<row_writer.m_nitem_max)
            row_writer.write(iitem++);
        mpi::barrier();
        log::debug_("ending HDF5 write loop over rows");
    }

public:

    void write(hdf5::GroupWriter &parent, std::string name, std::vector<std::string> field_names) const {
        RowHdf5Writer<row_t> row_writer(m_row, parent, name, nrow_to_write(), field_names);
        write_rows(row_writer);
    }

    void write(hdf5::GroupWriter &parent, std::string name) const {
        RowHdf5Writer<row_t> row_writer(m_row, parent, name, nrow_to_write());
        write_rows(row_writer);
    }

    virtual void read(hdf5::GroupReader &parent, std::string name) {
        RowHdf5Reader<row_t> row_reader(m_row, parent, name);
        size_t iitem = 0ul;
        clear();
        push_back(row_reader.m_nitem);
        for (row_reader.restart(); row_reader.in_range(); row_reader.step()){
            row_reader.read(iitem++);
        }
    }

private:
    Row &base_row() {
        return static_cast<Row &>(m_row);
    }
};


#endif //M7_TABLE_H
