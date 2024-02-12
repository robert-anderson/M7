//
// Created by Robert John Anderson on 12/02/2024.
//

#include "OpenAddresser.h"

OpenAddresser::OpenAddresser(const TableBase &table, size_t key_offset, size_t key_size, double fmax) :
        m_table(table), m_key_offset(key_offset), m_key_size(key_size), m_fmax(fmax),
        m_addrs(m_table.name()+" open addresser", 1, table.owner()),
        m_addrs_window(&m_addrs, sizeof(uint_t)){
    remap();
}

uint_t OpenAddresser::naddr() const {
    return m_addrs_window.m_nrow;
}

const buf_t *OpenAddresser::get_key(uint_t irow) const {
    return m_table.cbegin(irow) + m_key_offset;
}

uint_t OpenAddresser::get_addr(uint_t iaddr) const {
    return reinterpret_cast<const uint_t*>(m_addrs_window.cbegin())[iaddr];
}

void OpenAddresser::set_addr(uint_t iaddr, uint_t addr) {
    if (m_addrs_window.i_can_modify())
        reinterpret_cast<uint_t*>(m_addrs_window.begin())[iaddr] = addr;
}

bool OpenAddresser::eq(uint_t irow, const buf_t *key) const {
    return !std::memcmp(get_key(irow), key, m_key_size);
}

bool OpenAddresser::eq(uint_t irow, uint_t j_row) const {
    return eq(irow, get_key(j_row));
}

bool OpenAddresser::insert(uint_t irow, const buf_t *key) {
    DEBUG_ASSERT_TRUE(naddr(), "Must have a non-zero number of addresses allocated");
    const auto hash = hash::fnv(key, m_key_size) % naddr();

    auto insert_fn = [&](size_t iaddr) -> int {
        const auto addr = get_addr(iaddr);
        if (addr == ~0ul) {
            set_addr(iaddr, irow);
            return 1;
        }
        else if (eq(irow, addr)) return 0; // key exists
        // key neither found nor inserted
        return -1;
    };
    // probe from the hash position to the end of the map
    for (auto iaddr = hash; iaddr < naddr(); ++iaddr) {
        auto tmp = insert_fn(iaddr);
        if (tmp >= 0) return tmp;
    }
    // wrap around and start from the top
    for (auto iaddr = 0; iaddr < hash; ++iaddr) {
        auto tmp = insert_fn(iaddr);
        if (tmp >= 0) return tmp;
    }
    ABORT("Should always be able to insert into an open addressing hashmap");
    return false;
}

bool OpenAddresser::insert(uint_t irow) {
    return insert(irow, get_key(irow));
}

size_t OpenAddresser::lookup(const buf_t *key) {
    const auto hash = hash::fnv(key, m_key_size) % naddr();
    auto lookup_fn = [&](size_t iaddr) -> int {
        const auto addr = get_addr(iaddr);
        if (addr == ~0ul) return 0; // key not in map
        else if (eq(addr, key)) return 1; // key exists
        // key doesn't match addr
        return -1;
    };
    // probe from the hash position to the end of the map
    for (auto iaddr = hash; iaddr < naddr(); ++iaddr) {
        const auto tmp = lookup_fn(iaddr);
        if (tmp < 0) continue;
        else return tmp ? get_addr(iaddr) : ~0ul;
    }
    // wrap around and start from the top
    for (auto iaddr = 0; iaddr < hash; ++iaddr) {
        const auto tmp = lookup_fn(iaddr);
        if (tmp < 0) continue;
        else return tmp ? get_addr(iaddr) : ~0ul;
    }
    return ~0ul;
}

void OpenAddresser::resize() {
    const auto size = m_addrs_window.m_size;
    m_addrs_window.resize(sizeof(uint_t) * uint_t(double(m_table.capacity()) / m_fmax));
    if (m_addrs_window.i_can_modify())
        std::fill(m_addrs_window.begin() + size, m_addrs_window.begin() + m_addrs_window.m_size, 0xff);
}

void OpenAddresser::remap() {
    resize();
    m_addrs = 0xff;
    for (size_t irow = 0; irow < m_table.nrow_in_use(); ++irow) {
        if (m_table.is_freed(irow)) continue;
        insert(irow);
    }
    DEBUG_ASSERT_EQ(m_naddr_inserted, m_table.nrow_in_use(),
                    "All rows in use in the Table should have an associated entry in the hashmap");
}
