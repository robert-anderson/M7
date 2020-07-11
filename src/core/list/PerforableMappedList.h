//
// Created by Robert John Anderson on 2020-04-02.
//

#ifndef M7_PERFORABLEMAPPEDLIST_H
#define M7_PERFORABLEMAPPEDLIST_H

#include <stack>
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
    std::stack<size_t> m_free_rows;

public:

    PerforableMappedList(std::string name, Field_T& key_field, size_t nbucket):
        MappedList<T>(name, key_field, nbucket) {}

    size_t push(const T &key) override {
        size_t irow = ~0ul;
        // first see if there are any free rows left
        if (!m_free_rows.empty()){
            // free row was available, grab from top of stack
            irow = m_free_rows.top();
            m_free_rows.pop();
        }
        else {
            irow = List::push();
        }
        MappedList<T>::m_map.insert(key, irow);
        return irow;
    }

    size_t remove(const T &key, const size_t &irow) {
        auto irow_rm = MappedList<T>::m_map.remove(key, irow);
        ASSERT(irow != ~0ul);
        m_free_rows.push(irow_rm);
        MappedList<T>::m_key_field(irow_rm).zero();
        Table::zero_row(irow, 0);
        ASSERT(MappedList<T>::m_key_field(irow_rm).is_zero());
        ASSERT(irow_rm==irow)
        return irow_rm;
    }

    size_t nzero_rows(size_t isegment=0) const {
        // debugging only
        size_t result=0ul;
        for (size_t irow = 0ul; irow<MappedList<T>::high_water_mark(isegment); ++irow) {
            result += MappedList<T>::m_key_field(irow, isegment).is_zero();
        }
        return result;
    }

};


#endif //M7_PERFORABLEMAPPEDLIST_H
