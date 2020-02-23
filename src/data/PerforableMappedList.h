//
// Created by Robert John Anderson on 2020-02-23.
//

#ifndef M7_PERFORABLEMAPPEDLIST_H
#define M7_PERFORABLEMAPPEDLIST_H


#include "MappedList.h"

/*
 * The PerforableMappedList extends the MappedList by supporting a
 * deletion mechanism which does not leak rows.
 * In the base class, List, new elements are always inserted at the
 * high-water-mark, and so if deletion were allowed on such objects,
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
    defs::inds m_free_rows;
    size_t m_nfree_row;
    size_t m_free_rows_used;
    defs::inds m_deleted_rows;
    size_t m_ndeleted_row;
public:

    PerforableMappedList(Specification spec, size_t nrow, size_t key_entry) :
    MappedList<T>(spec, nrow, key_entry),
    m_free_rows(nrow, 0ul), m_nfree_row(0ul), m_free_rows_used(0ul),
    m_deleted_rows(nrow, 0ul), m_ndeleted_row(0){}

    PerforableMappedList(Specification spec, size_t nrow, size_t key_entry, defs::data_t *data_external) :
    MappedList<T>(spec, nrow, key_entry, data_external){}

    void next_cycle(){

    }

    size_t push(Mutex &mutex, const T& key){
        size_t irow;
        irow = MappedList<T>::m_map.lookup(mutex, key);
        if (irow!=~0ul) {
            return irow;
        }
        // key not found in table, so see if there are any free rows left
#pragma omp atomic capture
        irow = m_free_rows_used++;
        if (irow<m_nfree_row) {
            // free_row was available
            irow=m_free_rows[irow];
        }
        else {
            // free_rows is empty, need to push back high-water-mark
            irow = List::push();
        }
        MappedList<T>::m_map.insert(mutex, key, irow);
        return irow;
    }

    size_t push(const T& key) {
        auto mutex = MappedList<T>::m_map.find_mutex(key);
        return push(mutex, key);
    }


};


#endif //M7_PERFORABLEMAPPEDLIST_H
