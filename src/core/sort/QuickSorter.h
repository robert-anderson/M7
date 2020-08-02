#include <utility>

//
// Created by Robert John Anderson on 2020-08-02.
//

#ifndef M7_QUICKSORTER_H
#define M7_QUICKSORTER_H


#include <cstddef>
#include <functional>
#include <src/core/thread/ThreadPrivate.h>

template<typename T>
class QuickSorter {
protected:
    static const size_t c_cutoff = 1000;
    const size_t m_row_dsize;
    T *m_array;
    ThreadPrivate<std::vector<T>> m_tmp;
    ThreadPrivate<std::vector<T>> m_pivot;
    /*
     * comparator returns:
     * -1 if key(row1) < key(row2)
     *  1 if key(row1) > key(row2)
     *  0 if key(row1) ==key(row2)
     */
    virtual int cmp(const T *row1, const T *row2) = 0;

    T *row_ptr(const size_t &irow) {
        return m_array + irow * m_row_dsize;
    }

    void swap_rows(const size_t &irow, const size_t &jrow) {
        memcpy(m_tmp.get().data(), row_ptr(irow), m_row_dsize * sizeof(T));
        memcpy(row_ptr(irow), row_ptr(jrow), m_row_dsize * sizeof(T));
        memcpy(row_ptr(jrow), m_tmp.get().data(), m_row_dsize * sizeof(T));
    }

    void tasked_partition(const size_t &left, const size_t &right) {
        size_t i = left, j = right;
        memcpy(m_pivot.get().data(), row_ptr((left + right) / 2), m_row_dsize * sizeof(T));

        while (i <= j) {
            while (cmp(row_ptr(i), m_pivot.get().data()) > 0) i++;
            while (cmp(row_ptr(j), m_pivot.get().data()) < 0) j--;
            if (i <= j) {
                swap_rows(i, j);
                i++;
                j--;
            }
        }

        if (((right - left) < c_cutoff)) {
            if (left < j) { tasked_partition(left, j); }
            if (i < right) { tasked_partition(i, right); }

        } else {
#pragma omp task
            { tasked_partition(left, j); }
#pragma omp task
            { tasked_partition(i, right); }
        }
    }

public:
    QuickSorter(const size_t &row_dsize, T *array):
        m_row_dsize(row_dsize), m_array(array), m_tmp(row_dsize) {}

    void sort(const size_t &nrow) {
#pragma omp parallel num_threads(omp_get_max_threads())
        {
#pragma omp single nowait
            {
                tasked_partition(0, nrow - 1);
            }
        }

    }

};

#endif //M7_QUICKSORTER_H
