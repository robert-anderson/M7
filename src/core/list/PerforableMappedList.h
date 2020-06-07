//
// Created by Robert John Anderson on 2020-04-02.
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
    typedef typename T::Field_T Field_T;
    defs::inds m_free;
    defs::inds m_removed;

    size_t m_nfree = 0ul;
    size_t m_nfree_used = 0ul;
    size_t m_nremoved = 0ul;
public:

    PerforableMappedList(Field_T& key_field, size_t nbucket):
        MappedList<T>(key_field, nbucket) {}

    void synchronize() {
        /*
         * if the number of attempted accesses of the free stack exceeded the number
         * of elements in the stack, then the stack was emptied.
         * Now move the indices of the newly removed rows into the free stack
         */
        //std::cout << "PerforableMappedList high water mark: " << List::high_water_mark(0) <<std::endl;
        //std::cout << "PerforableMappedList free rows available: " << m_nfree <<std::endl;
        //std::cout << "PerforableMappedList free rows used: " << m_nfree_used <<std::endl;
        m_nfree = m_nfree_used>m_nfree? 0 : m_nfree-m_nfree_used;
        //std::cout << "PerforableMappedList free rows left over: " << m_nfree <<std::endl;
        //std::cout << "PerforableMappedList rows removed: " << m_nremoved <<std::endl;
        //std::cout << "PerforableMappedList zero rows: " << nzero_rows(0) <<std::endl;
        std::move(m_removed.begin(), m_removed.begin()+m_nremoved, m_free.begin()+m_nfree);
        m_nfree+=m_nremoved;
        //ASSERT(m_nfree == nzero_rows(0));
        m_nremoved = 0ul;
        m_nfree_used = 0ul;
    }

    size_t push(Mutex &mutex, const T &key) override {
        size_t irow = ~0ul;
        // first see if there are any free rows left
        if (m_nfree_used<m_nfree) {
#pragma omp atomic capture
            irow = ++m_nfree_used;
            if (irow <= m_nfree) {
                // free_row was available, grab from back of stack
                irow = m_free[m_nfree-irow];
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

    size_t push(const T &key) override {
        auto mutex = MappedList<T>::m_map.key_mutex(key);
        return push(mutex, key);
    }

    size_t remove(Mutex &mutex, const size_t &key_index) {
        auto irow = MappedList<T>::m_map.remove(mutex, key_index);
        ASSERT(irow != ~0ul);
        size_t iremoved;
#pragma omp atomic capture
        iremoved = m_nremoved++;
        ASSERT(iremoved < m_removed.size());
        m_removed[iremoved] = irow;
        MappedList<T>::m_key_field(irow).zero();
        Table::zero_row(irow, 0);
        ASSERT(MappedList<T>::m_key_field(irow).is_zero());
        return irow;
    }

    size_t remove(const T &key, const size_t &key_index) {
        auto mutex = MappedList<T>::key_mutex(key);
        return remove(mutex, key_index);
    }

    size_t remove(const T &key) {
        auto mutex = MappedList<T>::key_mutex(key);
        auto key_index = MappedList<T>::lookup(mutex, key);
        if (key_index==~0ul) return key_index;
        return remove(mutex, key_index);
    }

    size_t nfilled() const{
        return MappedList<T>::high_water_mark(0)-(m_nfree+m_nremoved);
    }

    size_t nzero_rows(size_t isegment=0) const {
        // debugging only
        size_t result=0ul;
        for (size_t irow = 0ul; irow<MappedList<T>::high_water_mark(isegment); ++irow) {
            result += MappedList<T>::m_key_field(irow, isegment).is_zero();
        }
        return result;
    }

    void expand(size_t delta_rows) override {
        Table::expand(delta_rows);
        m_free.resize(m_free.size()+delta_rows);
        m_removed.resize(m_removed.size()+delta_rows);
    }

};


#endif //M7_PERFORABLEMAPPEDLIST_H
