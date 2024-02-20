//
// Created by Robert John Anderson on 14/02/2024.
//

#ifndef M7_OPENADDRESSEDTABLE_H
#define M7_OPENADDRESSEDTABLE_H

#include <M7_lib/field/Fields.h>
#include "Table.h"
#include "OpenAddresser.h"

template<typename row_t>
struct OpenAddressedTable : Table<row_t> {
    /**
     * won't compile unless the row defines a key_field_t;
     */
    typedef typename row_fields::Key<row_t>::type key_field_t;

    using Table<row_t>::to_string;
    using Table<row_t>::m_row;
    using TableBase::m_bw;

    row_t m_lookup_row;
    row_t m_insert_row;
    row_t m_erase_row;

    OpenAddresser m_oa;

    OpenAddressedTable(const row_t &row, double fmax) :
            Table<row_t>(row), m_lookup_row(m_row), m_insert_row(m_row), m_erase_row(m_row),
            m_oa(*this, (*this).owner(), row_fields::key(m_row).row_offset(), row_fields::key(m_row).m_size, fmax){}

    OpenAddressedTable(const row_t &row) : OpenAddressedTable(row, 0.5){}

    OpenAddressedTable& operator=(const OpenAddressedTable& other) {
        Table<row_t>::operator=(other);
        return *this;
    }

    OpenAddressedTable(const OpenAddressedTable<row_t> &other) : OpenAddressedTable(other.m_row, other.m_oa.m_fmax) {}

    bool operator==(const MappedTable<row_t> &other) const {
        return TableBase::operator==(other);
    }

    /**
     * attempt to find the element identified by key in the hash map.
     * @param key
     *  key to lookup
     * @return
     *  Row moved to the found position if it is found else it's moved to a null position
     */
    void lookup(const key_field_t &key, const row_t& row) const {
        const auto i_row = m_oa.lookup(key.cbegin());
        if (i_row == ~0ul) static_cast<const Row&>(row).select_null();
        else static_cast<const Row&>(row).jump(i_row);
    }

    const row_t& lookup(const key_field_t& key) const {
        lookup(key, m_lookup_row);
        return m_lookup_row;
    }

    row_t& lookup(const key_field_t& key) {
        lookup(key, m_lookup_row);
        return m_lookup_row;
    }

    void clear() override {
        TableBase::clear();
        m_oa.clear();
    }

    void clear_row(row_t& row) {
        DEBUG_ASSERT_TRUE(Table<row_t>::associated(row), "row object not associated with this table");
        auto lookup = this->lookup(row_fields::key(row), row);
        erase(lookup);
    }

    /**
     * remove the record pointed to by row by:
     *  1. erasing from the OpenAddresser
     *  1. physically clearing by zeroing its record in the buffer,
     *  2. adding its index to the empty records stack
     * @param lookup
     *  points to the bucket and element within it to be deleted
     */
    void erase(const key_field_t &key) {
        lookup(key, m_erase_row);
        if (!m_erase_row) return;
        m_oa.erase(m_erase_row.index());
        /*
         * do steps 2 and 3.
         */
        TableBase::free(m_erase_row.index());
    }

    /**
     * insert into the hash table a new mapped element defined by the key arg
     * @param key
     *  hash table key object
     * @param row
     *  row object which points to the newly inserted record upon return of the method
     */
    void insert(const key_field_t &key, row_t& row) {
        DEBUG_ASSERT_TRUE(Table<row_t>::associated(row), "row object not associated with this table");
        const auto cap = TableBase::capacity();
        auto irow = TableBase::get_free_row();
        row.jump(irow);
        row.key_field() = key;
        if (TableBase::capacity() > cap) m_oa.remap();
        m_oa.insert(irow);
    }

    row_t& insert(const key_field_t &key) {
        insert(key, m_insert_row);
        return m_insert_row;
    }

    void insert(const row_t& src, row_t& row) {
        insert(row_fields::key(src), row);
        row.copy_in(src);
    }

    row_t& insert(const row_t& src) {
        insert(src, m_insert_row);
        return m_insert_row;
    }


    row_t& lookup_or_insert(const key_field_t& key) {
        lookup(key, m_lookup_row);
        if (m_lookup_row) return m_lookup_row;
        return insert(key);
    }

    row_t& lookup_or_insert(const row_t& src) {
        lookup(src, m_lookup_row);
        if (m_lookup_row) return m_lookup_row;
        return insert(src);
    }

    void resize(uint_t nrec, double factor=-1.0) override {
        TableBase::resize(nrec, factor);
        m_oa.remap();
    }
};

namespace buffered {
    template<typename row_t>
    struct OpenAddressedTable : BufferedTable<row_t, ::OpenAddressedTable<row_t>> {
        OpenAddressedTable(str_t name, const row_t &row, double fmax, Owner owner = Owner::local()) :
            BufferedTable<row_t, ::OpenAddressedTable<row_t>>(name, ::OpenAddressedTable<row_t>(row, fmax), owner) {}
        OpenAddressedTable(const row_t &row, Owner owner = Owner::local()) :
            OpenAddressedTable("", row, owner) {}
        OpenAddressedTable(str_t name, const row_t &row, Owner owner = Owner::local()) :
            OpenAddressedTable(name, row, 0.5, owner) {}
    };
}


#endif //M7_OPENADDRESSEDTABLE_H
