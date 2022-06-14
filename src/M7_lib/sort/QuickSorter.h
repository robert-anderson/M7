//
// Created by Robert J. Anderson on 29/11/2020.
//

#ifndef M7_QUICKSORTER_H
#define M7_QUICKSORTER_H

#include <functional>

#include <M7_lib/defs.h>
#include <M7_lib/table/Table.h>

template<typename Fn>
struct QuickSorter {

    defs::inds m_inds;
    Fn m_cmp_fn;

    QuickSorter(Fn cmp_fn): m_cmp_fn(cmp_fn){}

    const size_t& operator[](const size_t& i) const {
        ASSERT(i < m_inds.size());
        return m_inds[i];
    }

    void preserve_sort(const size_t &hwm){
        if (m_inds.size() < hwm) m_inds.reserve(hwm);
        m_inds.clear();
        // reset ordering
        for (size_t i = 0; i < hwm; ++i) m_inds.push_back(i);
        qs(0, hwm - 1);
        ASSERT(is_preserve_sorted(hwm));
    }

    void preserve_sort(const TableBase &table){
        preserve_sort(table.m_hwm);
    }

    bool is_preserve_sorted(const size_t &hwm) {
        for (size_t irow = 1ul; irow < hwm; ++irow) {
            if (m_cmp_fn(m_inds[irow], m_inds[irow-1]) &&!m_cmp_fn(m_inds[irow - 1], m_inds[irow])) return false;
        }
        return true;
    }

    bool is_preserve_sorted(const TableBase &table){
        return is_preserve_sorted(table.m_hwm);
    }

    void reorder_sort(TableBase &table){
        qs(0, table.m_hwm - 1, table);
        ASSERT(is_reorder_sorted(table));
    }

    /**
     * run through rows and return false if an example is found where two consecutive rows do not conform to the
     * ordering prescribed by m_cmp_fn.
     * e.g. if the sorting criteria is "ascending", and the consecutive values are (1, 2), then this is test passes
     * since the possible criteria which express this ordering (v_i > v_i-1 and v_i >= v_i-1) are both met
     *
     * @param table
     *  table to check for correct ordering
     * @return
     *  true if the rows in the given table are sorted according to m_cmp_fn
     */
    bool is_reorder_sorted(const TableBase &table){
        for (size_t irow = 1ul; irow < table.m_hwm; ++irow) {
            if (m_cmp_fn(irow, irow-1) && !m_cmp_fn(irow - 1, irow)) return false;
        }
        return true;
    }

private:
    void swap(size_t ii1, size_t ii2){
        auto i2 = m_inds[ii2];
        m_inds[ii2] = m_inds[ii1];
        m_inds[ii1] = i2;
    }

    size_t partition(size_t iilo, size_t iihi) {
        auto ip = m_inds[iihi];
        auto ii = iilo - 1;

        for (size_t ij = iilo; ij <= iihi - 1; ij++) {
            if (m_cmp_fn(m_inds[ij], ip)) {
                ii++;
                swap(ii, ij);
            }
        }
        swap(ii + 1, iihi);
        return ii + 1;
    }

    void qs(size_t iilo, size_t iihi) {
        if (iihi != ~0ul && iilo < iihi) {
            auto iip = partition(iilo, iihi);
            qs(iilo, iip - 1);
            qs(iip + 1, iihi);
        }
    }

    void swap(size_t ii1, size_t ii2, TableBase &table) {
        table.swap_rows(ii1, ii2);
    }

    size_t partition(size_t iilo, size_t iihi, TableBase &table) {
        auto ip = iihi;
        auto ii = iilo - 1;

        for (size_t ij = iilo; ij <= iihi - 1; ij++) {
            if (m_cmp_fn(ij, ip)) {
                ii++;
                swap(ii, ij, table);
            }
        }
        swap(ii + 1, iihi, table);
        return ii + 1;
    }

    void qs(size_t iilo, size_t iihi, TableBase &table){
        if (iihi != ~0ul && iilo < iihi) {
            auto iip = partition(iilo, iihi, table);
            qs(iilo, iip - 1, table);
            qs(iip + 1, iihi, table);
        }
    }
};

struct LambdaQuickSortFn {
    comparators::index_cmp_fn_t m_fn;
    LambdaQuickSortFn(comparators::index_cmp_fn_t fn): m_fn(fn){}
    bool operator()(const size_t& i, const size_t& j){
        return m_fn(i, j);
    }
};

struct LambdaQuickSorter2 : QuickSorter<LambdaQuickSortFn> {
    LambdaQuickSorter2(comparators::index_cmp_fn_t fn): QuickSorter<LambdaQuickSortFn>(fn){}
};


#endif //M7_QUICKSORTER_H
