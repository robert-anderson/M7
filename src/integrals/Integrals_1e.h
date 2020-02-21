//
// Created by Robert John Anderson on 2020-01-17.
//

#ifndef M7_INTEGRALS_1E_H
#define M7_INTEGRALS_1E_H

#include <cstddef>
#include <vector>
#include "Integrals.h"
#include "../io/FcidumpFileIterator.h"
#include "../consts.h"

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
        m_data.assign(m_nelem, 0.0);
    }

    Integrals_1e(std::string fname) :
            Integrals_1e(FcidumpFileIterator<T>(fname).m_norb, FcidumpFileIterator<T>(fname).m_spin_resolved) {
        FcidumpFileIterator<T> file_iterator(fname);
        defs::inds inds(4);
        T value;
        while (file_iterator.next(inds, value)) {
            if (valid_inds(inds)) {
                set(inds[0], inds[1], value);
            }
        }
    }

    inline size_t flat_index(const size_t &i, const size_t &j) const {
        if (isym == 1) return i + j * m_norb;
        else return i <= j ? trig(i, j) : trig(j, i);
    }

    void set(const size_t &i, const size_t &j, const T &value) {
        auto iflat = flat_index(i, j);
        auto conjd_value = (isym == 2 && i > j) ? consts::conj(value) : value;
        if (consts::float_is_zero(m_data[iflat])) m_data[iflat] = conjd_value;
        else {
            assert(consts::floats_nearly_equal(m_data[iflat], conjd_value));
        }
    }

    void set(const size_t &ispat, const size_t &ispin,
             const size_t &jspat, const size_t &jspin,
             const T &value) {
        set(spinorb(ispat, ispin), spinorb(jspat, jspin), value);
    }

    void set(const defs::inds &inds, const T &value) {
        assert(inds.size() == 4);
        assert(inds[2] == ((size_t) -1) && inds[3] == ((size_t) -1));
        set(inds[0], inds[1], value);
    }

    void set_from_fcidump(const defs::inds &inds, const T &value, bool spin_major = false) {
        /*
         * spin-resolved FCIDUMPs index in spinorbs, which may not may not be spin-major,
         * depending on the program they were generated for. E.g. NECI uses spatial-major
         * ordering throughout, so if the FCIDUMP supplied was intended for use with NECI,
         * spin_major should be passed in as false.
         */
        if (!m_spin_resolved) set(inds, value);
        else set(inds[0] / 2, inds[0] % 2, inds[1] / 2, inds[1] % 2, value);
    }

    T get(const size_t &i, const size_t &j) const {
        auto iflat = m_spin_resolved ? flat_index(i, j) : flat_index(i / 2, j / 2);
        return (isym == 2 && i > j) ? consts::conj(m_data[iflat]) : m_data[iflat];
    }

    T get(const size_t &ispat, const size_t &ispin,
          const size_t &jspat, const size_t &jspin) const {
        return get(spinorb(ispat, ispin), spinorb(jspat, jspin));
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
        for (auto i{0ul}; i < m_norb; ++i) {
            for (auto j{0ul}; i < m_norb; ++i) {
                if (!consts::float_is_zero(get(i, 0, j, 1))) return false;
            }
        }
        return true;
    }

private:
    constexpr size_t nelem(const size_t &norb) {
        static_assert(isym == 1 || isym == 2, "Invalid symmetry parameter specified.");
        if (isym == 1) {
            return norb * norb;
        } else if (isym == 2) {
            return trig(0, norb);
        }
    }
};

#endif //M7_INTEGRALS_1E_H
