//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_PERFORABLEMAPPEDTABLE_H
#define M7_PERFORABLEMAPPEDTABLE_H

#include <src/data/MappedTable.h>
#include <stack>

template<typename T>
class PerforableMappedTable : public MappedTable<T> {
    using Table::m_segment_mutex;
    using MappedTable<T>::m_ihashentry;
    using MappedTable<T>::m_map;
    std::vector<std::stack<size_t>> m_stack;

public:
    PerforableMappedTable(Specification spec, size_t nrow_initial, size_t m_nsegment = 1,
                          float nrow_growth_factor = 2.0, size_t nrow_mutex_blocks = 0, size_t ihashentry = 0) :
            MappedTable<T>(spec, nrow_initial, m_nsegment, nrow_growth_factor, nrow_mutex_blocks, ihashentry),
            m_stack{m_nsegment} {}

    size_t remove(const size_t &isegment, const size_t &irow) {
        if (irow != ~0ul) {
            m_stack[isegment].push(irow);
            auto key_view = Table::view<T>(isegment, irow, m_ihashentry);
            m_map[isegment].erase(key_view);
            key_view.zero();
        }
        return irow;
    }

    size_t remove(const size_t &isegment, const T &key) {
        return remove(isegment, MappedTable<T>::lookup(key, isegment));
    }

    size_t safe_remove(const size_t &isegment, const size_t &irow) {
        m_segment_mutex.acquire_lock(isegment);
        auto tmp = remove(isegment, irow);
        m_segment_mutex.release_lock(isegment);
        return tmp;
    }

    size_t safe_remove(const size_t &isegment, const T &key) {
        m_segment_mutex.acquire_lock(isegment);
        auto tmp = remove(isegment, key);
        m_segment_mutex.release_lock(isegment);
        return tmp;
    }

    size_t push(const size_t &isegment, const T &key) {
        size_t irow;
        if (m_stack[isegment].empty()) {
            irow = MappedTable<T>::push(isegment, key);
        } else {
            irow = m_stack[isegment].top();
            m_stack[isegment].pop();
        }
        return irow;
    }

    size_t safe_push(const size_t &isegment, const T &key) {
        size_t irow;
        if (m_stack[isegment].empty()) {
            irow = MappedTable<T>::safe_push(isegment, key);
        } else {
            m_segment_mutex.acquire_lock(isegment);
            irow = m_stack[isegment].top();
            m_stack[isegment].pop();
            m_segment_mutex.release_lock(isegment);
        }
        return irow;
    }

};


#endif //M7_PERFORABLEMAPPEDTABLE_H
