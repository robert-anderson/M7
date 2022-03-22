//
// Created by Robert John Anderson on 2020-01-17.
//

#ifndef M7_INTEGRALS_1E_H
#define M7_INTEGRALS_1E_H

#include <cstddef>
#include <vector>

#include <M7_lib/defs.h>
#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/io/FcidumpFileReader.h>
#include <M7_lib/util/consts.h>

#include "Integrals.h"

/*
 * notational conventions used are as follows:
 * i, j: these refer to spin-orbitals if the integrals are spin-resolved,
 *       otherwise these are spatial orbitals
 * ispin, jspin: these refer to the spin channel [0, 1] of the spin orbital,
 *       which is discarded in the non-spin-resolved case
 * isat, jspat: these index the spatial orbital or kramers pair [0, nspatorb)
 * flat_index: the position of the symmetrically unique integral value within m_buffer.
 */

template<typename T, size_t isym>
class Integrals_1e : public Integrals {
    //            none       h
    static_assert(isym == 1 || isym == 2, "Invalid symmetry parameter specified.");
    const size_t m_nelem;
    std::vector<T> m_data;

public:
    Integrals_1e(const size_t &nsite, bool spin_res) :
            Integrals(nsite, spin_res), m_nelem(nelem(m_nintind)) {
        m_data.resize(m_nelem, 0.0);
    }

    /*
    Integrals_1e(std::string fname) :
            Integrals_1e(FcidumpFileReader<T>(fname).m_nintind, FcidumpFileReader<T>(fname).m_spin_res) {
        FcidumpFileReader<T> file_iterator(fname);
        defs::inds inds(4);
        T value;
        while (file_iterator.next(inds, value)) {
            if (valid_inds(inds)) {
                set(inds[0], inds[1], value);
            }
        }
    }
     */

    inline size_t flat_index(const size_t &i, const size_t &j) const {
        if (isym == 1) return i + j * m_nintind;
        else return i <= j ? trig(i, j) : trig(j, i);
    }

    void set(const size_t &i, const size_t &j, const T &value) {
        auto iflat = flat_index(i, j);
        auto conjd_value = (isym == 2 && i > j) ? consts::conj(value) : value;
        if (consts::nearly_zero(m_data[iflat])) m_data[iflat] = conjd_value;
        else {
            DEBUG_ASSERT_NEARLY_EQ(m_data[iflat], conjd_value, defs::integral_tol,
                                   "conjugated value conflicts with current non-zero value");
        }
    }

    void set(const size_t &ispat, const size_t &ispin,
             const size_t &jspat, const size_t &jspin,
             const T &value) {
        set(spinorb(ispat, ispin), spinorb(jspat, jspin), value);
    }

    void set(const defs::inds &inds, const T &value) {
        DEBUG_ASSERT_EQ(inds.size(), 4ul, "incorrect number of indices");
        DEBUG_ASSERT_TRUE(inds[2] == ~0ul && inds[3] == ~0ul, "only the first two indices should be valid");
        set(inds[0], inds[1], value);
    }

private:
    T get(const size_t &i, const size_t &j) const {
        //auto iflat = m_spin_res ? flat_index(i, j) : flat_index(i / 2, j / 2);
        auto iflat = flat_index(i, j);
        return (isym == 2 && i > j) ? consts::conj(m_data[iflat]) : m_data[iflat];
    }

    T get(const size_t &ispat, const size_t &ispin,
          const size_t &jspat, const size_t &jspin) const {
        return get(spinorb(ispat, ispin), spinorb(jspat, jspin));
    }

public:
    T operator()(const size_t &i, const size_t &j) const {
        /*
         * return the one-body integral between the two SPINORS indexed by i and j
         */
        if (!m_spin_res && ((i < m_nintind) != (j < m_nintind))) return 0.0; // spin conservation
        auto iflat = m_spin_res ? flat_index(i, j) : flat_index(i % m_nintind, j % m_nintind);
        return (isym == 2 && i > j) ? consts::conj(m_data[iflat]) : m_data[iflat];
    }

    static bool valid_inds(const defs::inds &inds) {
        return inds[1] < ~0ul && inds[2] == ~0ul;
    }

    bool spin_conserving() const {
        if (!m_spin_res) return true;
        /*
         * iterate through integrals looking for an example of a non-zero
         * spin non-conserving one-body integral
         */
        for (size_t i = 0ul; i < m_nintind; ++i) {
            for (size_t j = 0ul; i < m_nintind; ++i) {
                if (!consts::nearly_zero(get(i, 0, j, 1))) return false;
            }
        }
        return true;
    }

private:
    static inline size_t nelem(const size_t &nintind) {
        static_assert(isym == 1 || isym == 2, "Invalid symmetry parameter specified.");
        if (isym == 1) {
            return nintind * nintind;
        } else if (isym == 2) {
            return trig(0, nintind);
        }
        return ~0ul;
    }
};

#endif //M7_INTEGRALS_1E_H