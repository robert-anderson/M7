//
// Created by Robert John Anderson on 2020-02-23.
//

#ifndef M7_PERFORABLEMAPPEDLIST_H
#define M7_PERFORABLEMAPPEDLIST_H


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
    using spec_t = typename MappedList<T>::spec_t;
    defs::inds m_free;
    defs::inds m_removed;

    size_t m_nfree = 0ul;
    size_t m_nfree_used = 0ul;
    size_t m_nremoved = 0ul;
public:

    PerforableMappedList(spec_t spec, size_t nrow, size_t key_entry) :
        MappedList<T>(spec, nrow, key_entry),
        m_free(nrow, 0ul), m_removed(nrow, 0ul) {}

    PerforableMappedList(spec_t spec, size_t nrow, size_t key_entry, defs::data_t *data_external) :
        MappedList<T>(spec, nrow, key_entry, data_external) {}

    void synchronize() {
        /*
         * if the number of attempted accesses of the free stack exceeded the number
         * of elements in the stack, then the stack was emptied.
         * Now move the indices of the newly removed rows into the free stack
         */
        m_nfree = m_nfree_used>m_nfree? 0 : m_nfree-m_nfree_used;
        std::move(m_removed.begin(), m_removed.end(), m_free.begin()+m_nfree);
        m_nfree+=m_nremoved;
        m_nremoved = 0ul;
        m_nfree_used = 0ul;
    }

    size_t push(Mutex &mutex, const T &key) {
        size_t irow = MappedList<T>::m_map.lookup(mutex, key);
        if (irow != ~0ul) return irow;
        /*
         * key not found in table, so see if there are any free rows left
         */
        if (m_nfree_used<m_nfree) {
#pragma omp atomic capture
            irow = m_nfree_used++;
            if (irow < m_nfree) {
                // free_row was available
                irow = m_free[irow];
            }
            else irow = ~0ul;
        }
        if (irow==~0ul) {
            // free_rows is empty, need to push back high-water-mark
            irow = List::push();
        }
        MappedList<T>::m_map.insert(mutex, key, irow);
        return irow;
    }

    size_t push(const T &key) {
        auto mutex = MappedList<T>::m_map.key_mutex(key);
        return push(mutex, key);
    }

    size_t remove(Mutex mutex, const T &key) {
        auto irow = MappedList<T>::m_map.remove(mutex, key);
        Table::zero(irow);
        size_t iremoved;
#pragma omp atomic capture
        iremoved = m_nremoved++;
        assert(iremoved<m_removed.size());
        m_removed[iremoved] = irow;
        return irow;
    }

    size_t remove(const T &key) {
        auto mutex = MappedList<T>::key_mutex(key);
        return remove(mutex, key);
    }


};


#endif //M7_PERFORABLEMAPPEDLIST_H
