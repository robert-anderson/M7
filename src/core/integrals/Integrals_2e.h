//
// Created by Robert John Anderson on 2020-01-10.
//

#ifndef M7_INTEGRALS_2E_H
#define M7_INTEGRALS_2E_H

#include <cstddef>
#include <vector>
#include <unordered_map>
#include <array>
#include "src/core/util/consts.h"
#include "src/defs.h"
#include "src/core/io/FcidumpFileReader.h"
#include "src/core/parallel/SharedArray.h"
#include "Integrals.h"

/*
 * given the case, is a complex conjugation necessary?
 */
constexpr std::array<bool, 11> case_to_tconj{
        false, false, false, false, false, false, true, true, false, true, true};
//      1      2             4/8    isym


template<typename T, size_t isym>
class Integrals_2e : public Integrals {
    //            none       i          ih         ihr
    static_assert(isym == 1 || isym == 2 || isym == 4 || isym == 8, "Invalid symmetry parameter specified.");
    const size_t m_nintind2, m_nintind3, m_nelem_8fold, m_nelem;
    SharedArray<T> m_data;
public:
    Integrals_2e(const size_t &nsite, bool spin_res) :
            Integrals(nsite, spin_res), m_nintind2(m_nintind * m_nintind), m_nintind3(m_nintind2 * m_nintind),
            m_nelem_8fold(trig(0, trig(0, m_nintind))), m_nelem(nelem(m_nintind)), m_data(m_nelem){
    }

    /*
    * given the indices and their identified case, return the flat indices
    */
    inline size_t flat_index(const size_t &icase, const size_t &i, const size_t &j,
                             const size_t &k, const size_t &l) const {
        switch (icase) {
            case 0:
                return i + j * m_nintind + k * m_nintind2 + l * m_nintind3;
            case 1:
                return trig(i + j * m_nintind, k + l * m_nintind);
            case 2:
                return trig(k + l * m_nintind, i + j * m_nintind);
            case 3:
                return trig(trig(i, j), trig(k, l));
            case 4:
                return trig(trig(k, l), trig(i, j));
            case 5:
                return trig(trig(i, j), trig(l, k)) + (isym == 4 ? m_nelem_8fold : 0);
            case 6:
                return trig(trig(l, k), trig(i, j)) + (isym == 4 ? m_nelem_8fold : 0);
            case 7:
                return trig(trig(j, i), trig(k, l)) + (isym == 4 ? m_nelem_8fold : 0);
            case 8:
                return trig(trig(k, l), trig(j, i)) + (isym == 4 ? m_nelem_8fold : 0);
            case 9:
                return trig(trig(j, i), trig(l, k));
            case 10:
                return trig(trig(l, k), trig(j, i));
        }
        return ~0ul;
    }

    void set(const size_t &i, const size_t &j, const size_t &k, const size_t &l, const T &value) {
        auto icase = get_case(i, j, k, l);
        auto iflat = flat_index(icase, i, j, k, l);
        ASSERT(iflat<m_nelem)
        auto conjd_value = case_to_tconj[icase] ? consts::conj(value) : value;
        if (mpi::on_node_i_am_root()) {
            if (consts::nearly_zero(m_data[iflat])) m_data.set(iflat, conjd_value);
            else {
                DEBUG_ASSERT_NEARLY_EQ(m_data[iflat], conjd_value, defs::integral_tol,
                                       "conjugated value conflicts with current non-zero value");
            }
        }
    }

    void set(
            const size_t &ispat, const size_t &ispin,
            const size_t &jspat, const size_t &jspin,
            const size_t &kspat, const size_t &kspin,
            const size_t &lspat, const size_t &lspin, const T &value) {
        set(spinorb(ispat, ispin), spinorb(jspat, jspin),
            spinorb(kspat, kspin), spinorb(lspat, lspin), value);
    }

    void set(const defs::inds &inds, const T &value) {
        DEBUG_ASSERT_EQ(inds.size(), 4ul, "incorrect number of indices");
        set(inds[0], inds[1], inds[2], inds[3], value);
    }

    T get(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
        DEBUG_ASSERT_TRUE(i < m_nintind && j < m_nintind && k < m_nintind && l < m_nintind,
                          "all indices must be in range");
        auto icase = get_case(i, j, k, l);
        auto iflat = flat_index(icase, i, j, k, l);
        return case_to_tconj[icase] ? consts::conj(m_data[iflat]) : m_data[iflat];
    }

    T element(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
        /*
         * return spin-orbital indexed integral from indices in chemists' ordering
         */
        if (m_spin_res) return get(i, j, k, l);
        else {
            // enforce spin symmetry
            if (((i < m_nintind) != (j < m_nintind)) || ((k < m_nintind) != (l < m_nintind))) return 0.0;
            return get(i % m_nintind, j % m_nintind, k % m_nintind, l % m_nintind);
        }
    }

    T element(
            const size_t &ispat, const size_t &ispin,
            const size_t &jspat, const size_t &jspin,
            const size_t &kspat, const size_t &kspin,
            const size_t &lspat, const size_t &lspin) const {
        return element(spinorb(ispat, ispin), spinorb(jspat, jspin),
                   spinorb(kspat, kspin), spinorb(lspat, lspin));
    }

    T phys_element(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
        /*
         * return spin-orbital indexed integral from indices in physicists' ordering
         */
        return element(i, k, j, l);
    }

    T phys_element(
            const size_t &ispat, const size_t &ispin,
            const size_t &jspat, const size_t &jspin,
            const size_t &kspat, const size_t &kspin,
            const size_t &lspat, const size_t &lspin) {
        /*
         * phys_element(i, j, k, l) == < i j | \hat{g} | k l >
         */
        return phys_element(spinorb(ispat, ispin), spinorb(jspat, jspin),
                        spinorb(kspat, kspin), spinorb(lspat, lspin));
    }

    T phys_antisym_element(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
        return phys_element(i, j, k, l) - phys_element(i, j, l, k);
    }

    T phys_antisym_element(
            const size_t &ispat, const size_t &ispin,
            const size_t &jspat, const size_t &jspin,
            const size_t &kspat, const size_t &kspin,
            const size_t &lspat, const size_t &lspin) const {
        return phys_antisym_element(spinorb(ispat, ispin), spinorb(jspat, jspin),
                                spinorb(kspat, kspin), spinorb(lspat, lspin));
    }

    static bool valid_inds(defs::inds inds) {
        return inds[3] < ~0ul;
    }

private:
    static inline size_t nelem(const size_t &nintind) {
        static_assert(isym == 1 || isym == 2 || isym == 4 || isym == 8, "Invalid symmetry parameter specified.");
        switch (isym) {
            case 1 :
                return nintind * nintind * nintind * nintind;
            case 2 :
                return trig(0, nintind * nintind);
            case 4 :
                return 2 * trig(0, trig(0, nintind));
            case 8 :
                return trig(0, trig(0, nintind));
        }
        return ~0ul;
    }

    static inline size_t get_case(const size_t &i, const size_t &j, const size_t &k, const size_t &l) {
        if (isym == 1) {
            return 0;
        } else if (isym == 2) {
            /*
             * Only assuming I symmetry
             * rect(i, k) <= rect(j, l)
             */
            return (i < j || (i == j && k < l)) ? 1 : 2;
        } else if (isym > 2) {
            /*
             * in IH and IHR symmetry, there are 8 identifiable cases
             */
            if (i <= j) {
                if (k <= l) return (trig(i, j) <= trig(k, l)) ? 3 : 4;
                else return (trig(i, j) <= trig(l, k)) ? 5 : 6;
            } else {
                if (k <= l) return (trig(j, i) <= trig(k, l)) ? 7 : 8;
                else return (trig(j, i) <= trig(l, k)) ? 9 : 10;
            }
        }
        return ~0ul;
    }
};


#endif //M7_INTEGRALS_2E_H
