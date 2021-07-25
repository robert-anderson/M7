//
// Created by rja on 23/07/2021.
//

#ifndef M7_FRMONVCONNECTION_H
#define M7_FRMONVCONNECTION_H

#include "src/core/field/Fields.h"

class FrmOpString {
    defs::inds m_inds;
public:
    FrmOpString(size_t nsite) {
        m_inds.reserve(2*nsite);
    }

    operator const defs::inds &() const {
        return m_inds;
    }

    bool operator==(const defs::inds &v) const {
        return m_inds==v;
    }

    const size_t &operator[](const size_t &i) const {
        DEBUG_ASSERT_LT(i, m_inds.size(), "operator index is OOB");
        return m_inds[i];
    }

    void clear() {
        m_inds.clear();
    }

    size_t size() const {
        return m_inds.size();
    }

    size_t capacity() const {
        return m_inds.capacity();
    }

    void add(const size_t &i) {
        DEBUG_ASSERT_TRUE(m_inds.empty() || i > m_inds.back(),
                          "spin orbital indices should be added in ascending order");
        DEBUG_ASSERT_LT(i, capacity(), "spin orbital index is OOB");
        DEBUG_ASSERT_LT(size(), capacity(),
                        "should never have more fermion operators than spin orbitals");
        m_inds.push_back(i);
    }

    defs::inds::const_iterator cbegin() const {
        return m_inds.cbegin();
    }

    defs::inds::const_iterator cend() const {
        return m_inds.cend();
    }

    bool all_occ(const fields::FrmOnv &src) const {
        return std::any_of(cbegin(), cend(), [&src](size_t i) { return !src.get(i); });
    }

    bool all_vac(const fields::FrmOnv &src) const {
        return std::any_of(cbegin(), cend(), [&src](size_t i) { return src.get(i); });
    }

    bool is_sorted() const {
        return std::is_sorted(cbegin(), cend());
    }
};


struct FrmOnvConnection {
    /**
     * Annihilation and creation second quantized opertor strings defining the difference in occupation between src and
     * dst fermionic ONVs
     */
    FrmOpString m_ann, m_cre;
private:
    /**
     * efficient computation of phases without enumeration of common occupation requires the number of set bits between
     * the beginning of the bit string and the beginning of each each dataword to be cached
     */
    const size_t m_ndataword;
    std::vector<bool> m_dataword_phases;
public:

    explicit FrmOnvConnection(size_t nsite);

    void connect(const fields::FrmOnv &src, const fields::FrmOnv &dst);

    bool connect(const fields::FrmOnv &src, const fields::FrmOnv &dst, FrmOpString &com);

    void apply(const fields::FrmOnv &src, fields::FrmOnv &dst);

    bool apply(const fields::FrmOnv &src, FrmOpString &com);

    bool apply(const fields::FrmOnv &src, fields::FrmOnv &dst, FrmOpString &com);

    void clear();

    void add(const size_t &ann, const size_t &cre);

    void add(const size_t &ann1, const size_t &ann2, const size_t &cre1, const size_t &cre2);

    const defs::inds& ann() const;

    const defs::inds& cre() const;

private:
    void update_dataword_phases(const fields::FrmOnv &onv);

    /**
     * The individual phase of a bit position within a multi-word bit representation is the number of set bits
     * before that position.
     * @param onv
     * @param ibit
     * @return
     */
    bool independent_phase(const fields::FrmOnv &onv, const size_t &ibit);

public:
    /**
     * the overall phase of an excitation with respect to the "in" determinant is the product of the independent
     * phases only if the operators are applied in such a way that they do not interfere with one another
     * e.g. if the occupied set is [0, 1, 4, 6, 7, 9]
     * and the excitation is 9 -> 5
     * the independent phase of 9 is true, since there is an odd number of set bits before position 9
     * the independent phase of 5 is true, since there is an odd number of set bits before position 5
     * so overall, the phase is false, which is correct: the electron at position 9 moves past an even number of
     * others to reach its final position
     *
     * on the other hand, if the excitation is 4 -> 8
     * the independent phase of 4 is false, since there is an even number of set bits before position 4
     * the independent phase of 8 is true, since there is an odd number of set bits before position 8
     * so overall, the phase is true, which is clearly incorrect: the electron at position 4 moves past an even
     * number of others to reach its final position
     *
     * in the second example, the number of set positions before 8 is not correct since we have deleted an electron
     * at 4 first, so we need to compute the phase associated with doing the creation operation first.
     *
     * the algorithm is then to work through the creation and annihilation lists, and each time the smallest of the
     * current indices is a creation operator, negate the overall phase if there is an odd number of annihilations
     * remaining
     *
     * @param onv
     * @return
     */
    bool phase(const fields::FrmOnv &onv);
};


#endif //M7_FRMONVCONNECTION_H
