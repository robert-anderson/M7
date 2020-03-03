//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_TABLEARRAY_H
#define M7_TABLEARRAY_H

#include "src/defs.h"
#include "Table.h"
#include <iostream>

template<typename table_T>
class TableArray {
    static_assert(std::is_base_of<Table, table_T>::value, "Template type parameter must derive from Table");
    std::vector<defs::data_t> m_data{};
    std::vector<table_T> m_tables{};
    defs::inds m_offsets;
public:
    TableArray(size_t ntable, const Specification &spec, size_t nrow_per_table) :
        m_offsets(ntable) {
        for (size_t itable = 0ul; itable < ntable; ++itable)
            m_tables.emplace_back(spec, 0, nullptr);
        grow(nrow_per_table);
    }

    void grow(const size_t &nrow_per_table) {
        m_data.resize(m_tables.size() * nrow_per_table * row_length(), 0);
        /*
         * move tables in reverse order to avoid corruption due to overlap
         */
        for (size_t itable = m_tables.size() - 1; itable != ~0ul; --itable) {
            m_offsets[itable] = itable * nrow_per_table * row_length();
            m_tables[itable].grow(data() + m_offsets[itable], nrow_per_table);
        }
    }

    void print() {
        std::cout << std::endl << "Table Array (" << m_tables.size() << ") tables" << std::endl;
        for (size_t itable = 0ul; itable < m_tables.size(); ++itable) {
            std::cout << "Table " << itable << std::endl;
            m_tables[itable].print();
        }
    }

    table_T &operator[](const size_t itable) {
        assert(itable < m_tables.size());
        return m_tables[itable];
    }

    defs::inds high_water_marks() const {
        defs::inds out;
        for (auto &table:m_tables) out.push_back(table.high_water_mark());
        return out;
    }

    const size_t row_length() const {
        return m_tables[0].row_length();
    }

    defs::data_t *data() {
        return m_data.data();
    }

    defs::inds &offsets() {
        return m_offsets;
    }

    void zero(){
        for (auto &table:m_tables) table.zero();
    }
};


#endif //M7_TABLEARRAY_H
