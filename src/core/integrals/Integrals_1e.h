//
// Created by Robert John Anderson on 2020-01-17.
//

#ifndef M7_INTEGRALS_1E_H
#define M7_INTEGRALS_1E_H

#include <cstddef>
#include <vector>
#include <src/core/io/FcidumpFileReader.h>
#include "Integrals.h"
#include "src/core/util/consts.h"

/*
 * notational conventions used are as follows:
 * i, j: these refer to spin-orbitals if the integrals are spin-resolved,
 *       otherwise these are spatial orbitals
 * ispin, jspin: these refer to the spin channel [0, 1] of the spin orbital,
 *       which is discarded in the non-spin-resolved case
 * isat, jspat: these index the spatial orbital or kramers pair [0, nspatorb)
 * flat_index: the position of the symmetrically unique integral value within m_data.
 */

template<typename T, size_t isym>
class Integrals_1e : public Integrals {
    //            none       h
    static_assert(isym == 1 || isym == 2, "Invalid symmetry parameter specified.");
    const size_t m_nelem;
    std::vector<T> m_data;

public:
    Integrals_1e(const size_t &norb, bool spin_resolved) :
            Integrals(norb, spin_resolved), m_nelem(nelem(norb)) {
        m_data.resize(m_nelem, 0.0);
    }

    /*
    Integrals_1e(std::string fname) :
            Integrals_1e(FcidumpFileReader<T>(fname).m_norb, FcidumpFileReader<T>(fname).m_spin_resolved) {
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
        if (isym == 1) return i + j * m_norb;
        else return i <= j ? trig(i, j) : trig(j, i);
    }

    void set(const size_t &i, const size_t &j, const T &value) {
        auto iflat = flat_index(i, j);
        auto conjd_value = (isym == 2 && i > j) ? consts::conj(value) : value;
        if (consts::float_is_zero(m_data[iflat])) m_data[iflat] = conjd_value;
        else {
            ASSERT(consts::floats_nearly_equal(m_data[iflat], conjd_value));
        }
    }

    void set(const size_t &ispat, const size_t &ispin,
             const size_t &jspat, const size_t &jspin,
             const T &value) {
        set(spinorb(ispat, ispin), spinorb(jspat, jspin), value);
    }

    void set(const defs::inds &inds, const T &value) {
        ASSERT(inds.size() == 4);
        ASSERT(inds[2] == ~0ul && inds[3] == ~0ul);
        set(inds[0], inds[1], value);
    }

    T get(const size_t &i, const size_t &j) const {
        //auto iflat = m_spin_resolved ? flat_index(i, j) : flat_index(i / 2, j / 2);
        auto iflat = flat_index(i, j);
        return (isym == 2 && i > j) ? consts::conj(m_data[iflat]) : m_data[iflat];
    }

    T get(const size_t &ispat, const size_t &ispin,
          const size_t &jspat, const size_t &jspin) const {
        return get(spinorb(ispat, ispin), spinorb(jspat, jspin));
    }


    T operator()(const size_t &i, const size_t &j) const {
        /*
         * return the one-body integral between the two SPINORS indexed by i and j
         */
        if (!m_spin_resolved && ((i<m_norb)!=(j<m_norb))) return 0.0;
        auto iflat = m_spin_resolved ? flat_index(i, j) : flat_index(i % m_norb, j % m_norb);
        return (isym == 2 && i > j) ? consts::conj(m_data[iflat]) : m_data[iflat];
    }

    static bool valid_inds(const defs::inds &inds) {
        return inds[1] < (size_t) -1 && inds[2] == (size_t) -1;
    }

    bool spin_conserving() const {
        if (!m_spin_resolved) return true;
        /*
         * iterate through integrals looking for an example of a non-zero
         * spin non-conserving one-body integral
         */
        for (size_t i = 0ul; i < m_norb; ++i) {
            for (size_t j = 0ul; i < m_norb; ++i) {
                if (!consts::float_is_zero(get(i, 0, j, 1))) return false;
            }
        }
        return true;
    }

private:
    static const size_t nelem(const size_t &norb) {
        static_assert(isym == 1 || isym == 2, "Invalid symmetry parameter specified.");
        if (isym == 1) {
            return norb * norb;
        } else if (isym == 2) {
            return trig(0, norb);
        }
        return ~0ul;
    }
};

#endif //M7_INTEGRALS_1E_H
