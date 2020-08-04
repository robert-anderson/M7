//
// Created by Robert John Anderson on 2020-04-02.
//

#ifndef M7_PERFORABLEMAPPEDLIST_H
#define M7_PERFORABLEMAPPEDLIST_H

#include <stack>
#include <src/core/hash/ConcurrentStack.h>
#include "MappedList.h"

/*
 * The PerforableMappedList extends the MappedList by supporting a
 * removal mechanism which does not leak rows.
 * In the base class, List, new elements are always inserted at the
 * high-water-mark, and so if removal were allowed on such objects,
 * deleted entries would become abandoned rows.
 *
 * The obvious solution here is to employ a last-in-first-out stack
 * to keep track of the available rows below the high-water-mark,
 * however, ensuring threadsafety of such an object is tricky.
 *
 * It is noted, however that the present algorithm does not demand
 * a complete stack implementation, we only require that the container
 * does not leak rows unsustainably.
 *
 * Therefore, in this class the free row indices FROM THE PREVIOUS
 * ITERATION are stored in a vector, whose next is referenced by an
 * atomically updated counter. If this counter reaches or exceeds the
 * size of the free_rows, the push method will continue to increment
 * the high-water-mark.
 *
 * At the same time, another vector is filled up in the same manner, this
 * will constitute the free_rows of the next iteration.
 *
 */

template<typename T>
class PerforableMappedList : public MappedList<T> {
    typedef typename T::Field_T Field_T;
    ConcurrentStack<size_t> m_free_rows;

    using MappedList<T>::m_map;

public:

    PerforableMappedList(std::string name, Field_T &key_field, size_t nbucket) :
            MappedList<T>(name, key_field, nbucket) {}

    size_t push(const T &key) override {
        size_t irow = ~0ul;
        if (m_free_rows.pop(irow)) {
#ifdef VERBOSE_DEBUGGING
            std::cout << consts::verb << consts::chevs << "POPPED AN EMPTY ROW FROM THE FREE STACK" << std::endl;
            std::cout << consts::verb << "free rows remaining:    " << m_free_rows.size() << std::endl;
#endif
        } else {
            irow = List::push();
#ifdef VERBOSE_DEBUGGING
            std::cout << consts::verb << consts::chevs << "NO ROWS IN FREE STACK - PUSHING BACK THE HIGHWATER MARK"
                      << std::endl;
#endif
        }
        MappedList<T>::m_key_field(irow) = key;
        MappedList<T>::m_map.insert(key, irow);
        return irow;
    }

    void remove(const size_t &irow) {
        ASSERT(!MappedList<T>::m_key_field(irow, 0).is_zero());
        m_free_rows.append(irow);
        m_map.mark_for_delete(irow);
    }

#if 0
    void clear_tombstones() {
        m_map.clear_tombstones();
    }
#endif

    size_t nzero_rows(size_t isegment = 0) const {
        // debugging only
        size_t result = 0ul;
        for (size_t irow = 0ul; irow < MappedList<T>::high_water_mark(isegment); ++irow) {
            result += MappedList<T>::m_key_field(irow, isegment).is_zero();
        }
        return result;
    }

    size_t nrow_in_free_stack(size_t isegment = 0) const {
        // debugging only
        return m_free_rows.size();
    }

    void print() override {
        Table::print();
        std::cout << "Free Rows:" << utils::to_string(m_free_rows.to_vector()) << std::endl;
        std::cout << "Mapped Rows:" << utils::to_string(m_map.to_vector()) << std::endl;
    }

};


#endif //M7_PERFORABLEMAPPEDLIST_H
