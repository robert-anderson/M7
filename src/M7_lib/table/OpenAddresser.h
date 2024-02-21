//
// Created by Robert John Anderson on 12/02/2024.
//

#ifndef M7_OPENADDRESSER_H
#define M7_OPENADDRESSER_H

#include "M7_lib/util/Hash.h"
#include "M7_lib/table/Table.h"

struct OpenAddresser : Table<SingleFieldRow<field::Number<uint_t>>> {
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
     *  fmax = (number of slots for keys) / (number of rows allocated in the table)
     */
    const double m_fmax;

    OpenAddresser(const TableBase& table, Owner owner, size_t key_offset, size_t key_size, double fmax = 0.3);

private:
    const buf_t* get_key(uint_t irow) const;

    uint_t get_addr(uint_t iaddr) const;

    void set_addr(uint_t iaddr, uint_t addr);

    bool eq(uint_t irow, const buf_t* key) const;

    bool eq(uint_t irow, uint_t j_row) const;
    /**
     * Find the index of the key in the m_addrs table
     * @param key
     *  key to lookup in m_addrs
     * @return
     *  Row index in m_addrs corresponding to key if found, else ~0ul
     */
    size_t lookup_iaddr(const buf_t* key) const;

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
    size_t lookup(const buf_t* key) const;
    /**
     * Lazily delete the given key by assigning a tombstone value
     * @param key
     *  pointer to beginning of key to be inserted
     * @return
     *  true if the deletion was successful, false if the key does not exist
     */
    bool erase(const buf_t* key);
    /**
     * @param irow
     *  Row index to erase
     * @return
     *  true if the deletion was successful, false if the key does not exist
     */
    bool erase(uint_t irow);
    /**
     * Set the size of the map correctly given the capacity of the table and the max load factor
     */
    void resize();
    /**
     * resize the map and reenter all in-use rows of m_table
     */
    void remap();
    /**
     * set all addrs to the unused value ~0ul
     */
    void clear() override;
};



#endif //M7_OPENADDRESSER_H
