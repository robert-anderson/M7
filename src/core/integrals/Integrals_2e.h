//
// Created by Robert John Anderson on 2020-01-10.
//

#ifndef M7_INTEGRALS_2E_H
#define M7_INTEGRALS_2E_H

#include <cstddef>
#include <vector>
#include <unordered_map>
#include <array>
#include <assert.h>
#include "src/core/multidim/Indexer.h"
#include "src/consts.h"
#include "src/core/io/FcidumpFileIterator.h"
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
    const size_t m_norb2, m_norb3, m_nelem_8fold, m_nelem;
    std::vector<T> m_data;
public:
    Integrals_2e(const size_t &norb, bool spin_resolved) :
            Integrals(norb, spin_resolved), m_norb2(m_norb * norb), m_norb3(m_norb2 * norb),
            m_nelem_8fold(trig(0, trig(0, norb))), m_nelem(nelem(norb)) {
        m_data.assign(m_nelem, 0.0);
    }

    Integrals_2e(std::string fname, bool spin_major = false) :
            Integrals_2e(FcidumpFileIterator<T>(fname).m_norb, FcidumpFileIterator<T>(fname).m_spin_resolved) {
        FcidumpFileIterator<T> file_iterator(fname);
        defs::inds inds(4);
        T value;
        while (file_iterator.next(inds, value)) {
            if (valid_inds(inds)) {
                set_from_fcidump(inds, value, spin_major);
            }
        }
    }

    /*
    * given the indices and their identified case, return the flat indices
    */
    inline size_t flat_index(const size_t &icase, const size_t &i, const size_t &j,
                             const size_t &k, const size_t &l) const {
        switch (icase) {
            case 0:
                return i + j * m_norb + k * m_norb2 + l * m_norb3;
            case 1:
                return trig(i + j * m_norb, k + l * m_norb);
            case 2:
                return trig(k + l * m_norb, i + j * m_norb);
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
    }

    void set(const size_t &i, const size_t &j, const size_t &k, const size_t &l, const T &value) {
        auto icase = get_case(i, j, k, l);
        auto iflat = flat_index(icase, i, j, k, l);
        auto conjd_value = case_to_tconj[icase] ? consts::conj(value) : value;
        if (consts::float_is_zero(m_data[iflat])) m_data[iflat] = conjd_value;
        else {
            assert(consts::floats_nearly_equal(m_data[iflat], conjd_value));
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
        assert(inds.size() == 4);
        set(inds[0], inds[1], inds[2], inds[3], value);
    }

    void set_from_fcidump(const defs::inds &inds, const T &value, bool spin_major = false) {
        /*
         * spin-resolved FCIDUMPs index in spinorbs, which may not may not be spin-major,
         * depending on the program they were generated for. E.g. NECI uses spatial-major
         * ordering throughout, so if the FCIDUMP supplied was intended for use with NECI,
         * spin_major should be passed in as false.
         */
        if (!m_spin_resolved) set(inds, value);
        else if (spin_major)
            set(inds, value);
        else
            set(
                    inds[0] / 2, inds[0] % 2, inds[1] / 2, inds[1] % 2,
                    inds[2] / 2, inds[2] % 2, inds[3] / 2, inds[3] % 2, value
            );

    }

    T get(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
        assert(i<m_norb && j<m_norb && k<m_norb && l<m_norb);
        auto icase = get_case(i, j, k, l);
        auto iflat = flat_index(icase, i, j, k, l);
        return case_to_tconj[icase] ? consts::conj(m_data[iflat]) : m_data[iflat];
    }

    T element(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
        /*
         * return spin-orbital indexed integral from indices in chemists' ordering
         */
        if (m_spin_resolved) return get(i, j, k, l);
        else {
            // enforce spin symmetry
            if (((i<m_norb)!=(j<m_norb)) || ((k<m_norb)!=(l<m_norb))) return 0.0;
            return get(i%m_norb, j%m_norb, k%m_norb, l%m_norb);
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
    constexpr size_t nelem(const size_t &norb) {
        static_assert(isym == 1 || isym == 2 || isym == 4 || isym == 8, "Invalid symmetry parameter specified.");
        switch (isym) {
            case 1 :
                return norb * norb * norb * norb;
            case 2 :
                return trig(0, norb * norb);
            case 4 :
                return 2 * trig(0, trig(0, norb));
            case 8 :
                return trig(0, trig(0, norb));
        }
    }

    constexpr size_t get_case(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
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
    }
};


#endif //M7_INTEGRALS_2E_H
