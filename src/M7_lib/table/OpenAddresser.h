//
// Created by Robert John Anderson on 12/02/2024.
//

#ifndef M7_OPENADDRESSER_H
#define M7_OPENADDRESSER_H

#include "M7_lib/util/Hash.h"
#include "M7_lib/table/TableBase.h"

struct OpenAddresser {
    /**
     * Reference to the Table object on which this object provides access
     */
    const TableBase& m_table;
    /**
     * Offset of the key field from the begin of a given row
     */
    const uint_t m_key_offset;
    /**
     * Size of the key field in bytes
     */
    const uint_t m_key_size;
    /**
     * The maximum load factor.
     * The load factor of an open addressing hash map is:
     *  f = (number of slots for keys) / (number of inserted keys)
     * So the maximum such factor is:
     *  fmax = (number of slots for keys) / (number of rows allocated in the mapped table)
     */
    const double m_fmax;
    /**
     * Shared or local memory storage of the memory
     */
    Buffer m_addrs;
    Buffer::Window m_addrs_window;
    /**
     * number of inserted addresses
     */
    uint_t m_naddr_inserted;

    OpenAddresser(const TableBase& table, size_t key_offset, size_t key_size, double fmax = 0.3);

private:
    uint_t naddr() const;

    const buf_t* get_key(uint_t irow) const;

    uint_t get_addr(uint_t iaddr) const;

    void set_addr(uint_t iaddr, uint_t addr);

    bool eq(uint_t irow, const buf_t* key) const;

    bool eq(uint_t irow, uint_t j_row) const;

public:
    /**
     * @param irow
     *  row index of the table at which the key has either been written or is about to be written
     * @param key
     *  pointer to beginning of key to be inserted
     * @return
     *  true if the insertion was successful i.e. a new element was written to m_addrs, false if the key already exists
     */
    bool insert(uint_t irow, const buf_t* key);
    /**
     * Overload for the case that the key is already written to m_table.
     * @param irow
     *  row index of the table at which the key has been written
     * @return
     *  true if the insertion was successful i.e. a new element was written to m_addrs, false if the key already exists
     */
    bool insert(uint_t irow);

    /**
     * Find the key in the map and return its position within m_table.
     * @param key
     *  key to lookup in m_table
     * @return
     *  Row index corresponding to key if found, else ~0ul
     */
    size_t lookup(const buf_t* key);
    /**
     * Set the size of the map correctly given the capacity of the table and the max load factor
     */
    void resize();
    /**
     * resize the map and reenter all in-use rows of m_table
     */
    void remap();
};



#endif //M7_OPENADDRESSER_H
