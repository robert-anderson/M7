//
// Created by Robert J. Anderson on 29/11/2020.
//

#ifndef M7_QUICKSORT_H
#define M7_QUICKSORT_H

#include <functional>
#include <utility>

#include <M7_lib/defs.h>
#include <M7_lib/util/Functor.h>
#include <M7_lib/table/Table.h>

namespace quicksort {
    /**
     * Quicksort implementation with element comparison given by the statically-dispatched bool(uint_t, uint_t) function
     * provided
     * @tparam fn_t
     *  comparator functor type
     */
    template<typename fn_t>
    struct SorterFn {

        uintv_t m_inds;
        fn_t m_cmp_fn;

        SorterFn(fn_t cmp_fn) : m_cmp_fn(cmp_fn) {
            functor::assert_prototype<bool(uint_t, uint_t)>(m_cmp_fn);
        }

        const uint_t& operator[](const uint_t& i) const {
            ASSERT(i < m_inds.size());
            return m_inds[i];
        }

        void preserve_sort(uint_t nrow) {
            if (m_inds.size() < nrow) m_inds.reserve(nrow);
            m_inds.clear();
            // reset ordering
            for (uint_t i = 0; i < nrow; ++i) m_inds.push_back(i);
            qs(0, nrow - 1);
            ASSERT(is_preserve_sorted(nrow));
        }

        void preserve_sort(const TableBase& table) {
            preserve_sort(table.nrow_in_use());
        }

        bool is_preserve_sorted(uint_t nrow) {
            for (uint_t irow = 1ul; irow < nrow; ++irow) {
                if (m_cmp_fn(m_inds[irow], m_inds[irow - 1]) && !m_cmp_fn(m_inds[irow - 1], m_inds[irow])) return false;
            }
            return true;
        }

        bool is_preserve_sorted(const TableBase& table) {
            return is_preserve_sorted(table.nrow_in_use());
        }

        void reorder_sort(TableBase& table) {
            qs(0, table.nrow_in_use() - 1, table);
            REQUIRE_TRUE(is_reorder_sorted(table), "reorder sort unsuccessful");
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
        bool is_reorder_sorted(const TableBase& table) {
            for (uint_t irow = 1ul; irow < table.nrow_in_use(); ++irow) {
                if (m_cmp_fn(irow, irow - 1) && !m_cmp_fn(irow - 1, irow)) return false;
            }
            return true;
        }

    private:
        void swap(uint_t ii1, uint_t ii2) {
            auto i2 = m_inds[ii2];
            m_inds[ii2] = m_inds[ii1];
            m_inds[ii1] = i2;
        }

        uint_t partition(uint_t iilo, uint_t iihi) {
            auto ip = m_inds[iihi];
            auto ii = iilo - 1;

            for (uint_t ij = iilo; ij <= iihi - 1; ij++) {
                if (m_cmp_fn(m_inds[ij], ip)) {
                    ii++;
                    swap(ii, ij);
                }
            }
            swap(ii + 1, iihi);
            return ii + 1;
        }

        void qs(uint_t iilo, uint_t iihi) {
            if (iihi != ~0ul && iilo < iihi) {
                auto iip = partition(iilo, iihi);
                qs(iilo, iip - 1);
                qs(iip + 1, iihi);
            }
        }

        void swap(uint_t ii1, uint_t ii2, TableBase& table) {
            table.swap_records(ii1, ii2);
        }

        uint_t partition(uint_t iilo, uint_t iihi, TableBase& table) {
            auto ip = iihi;
            auto ii = iilo - 1;

            for (uint_t ij = iilo; ij <= iihi - 1; ij++) {
                if (m_cmp_fn(ij, ip)) {
                    ii++;
                    swap(ii, ij, table);
                }
            }
            swap(ii + 1, iihi, table);
            return ii + 1;
        }

        void qs(uint_t iilo, uint_t iihi, TableBase& table) {
            if (iihi != ~0ul && iilo < iihi) {
                auto iip = partition(iilo, iihi, table);
                qs(iilo, iip - 1, table);
                qs(iip + 1, iihi, table);
            }
        }
    };

    template<typename fn_t>
    void sort(TableBase& table, fn_t fn){
        SorterFn<fn_t> sorter(fn);
        sorter.reorder_sort(table);
    }
    /**
     * allows the above class to call a std::function, achieving type-erasure but at cost of slower execution
     */
    struct Adaptor {
        comparators::index_cmp_fn_t m_fn;
        Adaptor(comparators::index_cmp_fn_t fn): m_fn(std::move(fn)){}

        bool operator()(uint_t i, uint_t j) {
            return m_fn(i, j);
        }
    };

    struct Sorter : SorterFn<Adaptor> {
        Sorter(comparators::index_cmp_fn_t fn): SorterFn<Adaptor>(fn){}
    };
}


#endif //M7_QUICKSORT_H
