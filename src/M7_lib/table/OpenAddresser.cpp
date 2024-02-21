//
// Created by Robert John Anderson on 12/02/2024.
//

#include "OpenAddresser.h"

str_t name(const str_t& name) {
    return name.empty() ? "" : name + " open addresser";
}

OpenAddresser::OpenAddresser(const TableBase &table, size_t key_offset, size_t key_size, double fmax) :
        Table<SingleFieldRow<field::Number<uint_t>>>({}),
        m_table(table), m_key_offset(key_offset), m_key_size(key_size), m_fmax(fmax) {
    remap();
}

const buf_t *OpenAddresser::get_key(uint_t irow) const {
    return m_table.cbegin(irow) + m_key_offset;
}

uint_t OpenAddresser::get_addr(uint_t iaddr) const {
    m_row.jump(iaddr);
    return m_row.m_field;
}

void OpenAddresser::set_addr(uint_t iaddr, uint_t addr) {
    if (i_can_modify()) {
        m_row.jump(iaddr);
        m_row.m_field = addr;
    }
}

bool OpenAddresser::eq(uint_t irow, const buf_t *key) const {
    return !std::memcmp(get_key(irow), key, m_key_size);
}

bool OpenAddresser::eq(uint_t irow, uint_t j_row) const {
    return eq(irow, get_key(j_row));
}

size_t OpenAddresser::lookup_iaddr(const buf_t *key) const {
    if (empty()) return ~0ul;
    const auto hash = hash::fnv(key, m_key_size) % capacity();
    auto lookup_fn = [&](size_t iaddr) -> int {
        const auto addr = get_addr(iaddr);
        if (addr == ~0ul) return 0; // key not in map
        else if (addr == capacity() || !eq(addr, key)) return -1; // key doesn't match addr
        return 1; // key exists
    };
    // probe from the hash position to the end of the map
    for (auto iaddr = hash; iaddr < capacity(); ++iaddr) {
        const auto tmp = lookup_fn(iaddr);
        if (tmp < 0) continue;
        else return tmp ? iaddr : ~0ul;
    }
    // wrap around and start from the top
    for (auto iaddr = 0; iaddr < hash; ++iaddr) {
        const auto tmp = lookup_fn(iaddr);
        if (tmp < 0) continue;
        else return tmp ? iaddr : ~0ul;
    }
    return ~0ul;
}

bool OpenAddresser::insert(uint_t irow, const buf_t *key) {
    DEBUG_ASSERT_TRUE(capacity(), "Must have a non-zero number of addresses allocated");
    const auto hash = hash::fnv(key, m_key_size) % capacity();

    auto insert_fn = [&](size_t iaddr) -> int {
        const auto addr = get_addr(iaddr);
        if (addr >= capacity()) {
            set_addr(iaddr, irow);
            return 1;
        }
        else if (eq(irow, addr)) return 0; // key exists
        // key neither found nor inserted
        return -1;
    };
    // probe from the hash position to the end of the map
    for (auto iaddr = hash; iaddr < capacity(); ++iaddr) {
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

size_t OpenAddresser::lookup(const buf_t *key) const {
    const auto hash = hash::fnv(key, m_key_size) % capacity();
    auto iaddr = lookup_iaddr(key);
    return iaddr == ~0ul ? ~0ul : get_addr(iaddr);
}

bool OpenAddresser::erase(const buf_t *key) {
    auto iattr = lookup_iaddr(key);
    if (iattr >= capacity()) {
        // not found
        return false;
    }
    set_addr(iattr, capacity());
    return true;
}

bool OpenAddresser::erase(uint_t irow) {
    return erase(get_key(irow));
}

void OpenAddresser::resize() {
    if (!m_table.capacity()) return;
    auto size = push_back(uint_t(double(m_table.capacity()) / m_fmax));
    if (i_can_modify()) {
        // set all the newly created addrs to ~0ul (unused)
        for (m_row.restart(size); m_row; ++m_row) m_row.m_field = ~0ul;
    }
}

void OpenAddresser::remap() {
    resize();
    for (m_row.restart(); m_row; ++m_row) m_row.m_field = ~0ul;
    for (size_t irow = 0; irow < m_table.nrow_in_use(); ++irow) {
        if (m_table.is_freed(irow)) continue;
        insert(irow);
    }
}

void OpenAddresser::clear() {
    if (i_can_modify()) std::fill(begin(), begin() + bw_size(), 0xff);
}
