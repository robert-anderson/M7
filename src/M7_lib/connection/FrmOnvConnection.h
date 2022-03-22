//
// Created by rja on 23/07/2021.
//

#ifndef M7_FRMONVCONNECTION_H
#define M7_FRMONVCONNECTION_H

#include <M7_lib/field/FrmOnvField.h>

/**
 * A generic string of ordered, distinct spin-orbital indices.
 * the string is only modifiable via the clear and add methods, but (debug bounds checked) const access is provided to
 * the underlying index vector
 */
class FrmOps {
    const size_t m_nsite;
    defs::inds m_inds;
public:
    FrmOps(size_t nsite): m_nsite(nsite) {
        m_inds.reserve(2*nsite);
    }

    FrmOps(const FrmOps& other): FrmOps(other.m_nsite){}

    FrmOps& operator=(const FrmOps& other){
        DEBUG_ASSERT_EQ(other.m_inds.capacity(), m_inds.capacity(), "incompatible FrmOps instances");
        m_inds.assign(other.m_inds.cbegin(), other.m_inds.cend());
        return *this;
    }

    FrmOps& operator=(const defs::inds& inds){
        DEBUG_ASSERT_LT(inds.size(), m_inds.capacity(), "incompatible FrmOps instances");
        m_inds.assign(inds.cbegin(), inds.cend());
        return *this;
    }

    const defs::inds & inds() const {
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

    void add(std::pair<size_t, size_t> pair) {
        add(FrmOnvField::ibit(pair.first, pair.second, m_nsite));
    }

    void set(const defs::inds& orbs, const defs::inds& inds){
        clear();
        for (const auto& ind: inds) add(orbs[ind]);
        DEBUG_ASSERT_EQ(inds.size(), size(), "not all selected inds were added");
    }

    void add(std::pair<size_t, size_t> pair1, std::pair<size_t, size_t> pair2){
        add(pair1);
        add(pair2);
    }

    void set(std::pair<size_t, size_t> pair1, std::pair<size_t, size_t> pair2){
        clear();
        add(pair1, pair2);
    }

    defs::inds::const_iterator cbegin() const {
        return m_inds.cbegin();
    }

    defs::inds::const_iterator cend() const {
        return m_inds.cend();
    }
    /**
     * @param src
     *  the fermion ONV to check the occupation of
     * @return
     *  true if all elements of the index array are set bits in the ONV
     */
    bool all_occ(const FrmOnvField &src) const {
        return std::all_of(cbegin(), cend(), [&src](size_t i) { return src.get(i); });
    }
    /**
     * @param src
     *  the fermion ONV to check the occupation of
     * @return
     *  true if all elements of the index array are clear bits in the ONV
     */
    bool all_vac(const FrmOnvField &src) const {
        return std::all_of(cbegin(), cend(), [&src](size_t i) { return !src.get(i); });
    }
    /**
     * @return
     *  true if the stored indices are in ascending order and none are the same
     */
    bool is_valid() const {
        if (!size()) return true;
        if (!std::is_sorted(cbegin(), cend())) return false;
        for (auto it = cbegin()+1; it!=cend(); ++it) if (*it==*(it-1)) return false;
        return true;
    }

    size_t nalpha() const {
        size_t n=0ul;
        for (auto i: m_inds) n+=i<m_nsite;
        return n;
    }
};

/**
 * Stores two FrmOpProduct representing the difference in occupation between two determinants
 */
struct FrmOnvConnection {
    /**
     * Annihilation and creation second quantized opertor strings defining the difference in occupation between src and
     * dst fermionic ONVs
     */
    FrmOps m_ann, m_cre;
private:
    /**
     * efficient computation of phases without enumeration of common occupation requires the number of set bits between
     * the beginning of the bit string and the beginning of each each dataword to be cached
     */
    const size_t m_ndataword;
    mutable std::vector<bool> m_dataword_phases;
public:

    explicit FrmOnvConnection(size_t nsite);

    explicit FrmOnvConnection(BasisData bd);

    explicit FrmOnvConnection(const FrmOnvField& mbf);
    /**
     * update the internal state of the connection with the difference in occupation between the two given determinants
     * @param src
     *  fermion ONV to excite from
     * @param dst
     *  fermion ONV to excite to
     */
    void connect(const FrmOnvField &src, const FrmOnvField &dst);
    /**
     * update the internal state of the connection with the difference in occupation between the two given determinants
     * and populate the given FrmOpProduct with the indices occupied in both the src and dst
     * @param src
     *  fermion ONV to excite from
     * @param dst
     *  fermion ONV to excite to
     * @param com
     *  string of "common" spin orbital indices
     * @return
     *  antisymmetric phase of the connection
     */
    bool connect(const FrmOnvField &src, const FrmOnvField &dst, FrmOps &com);
    /**
     * update the dst determinant given the src determinant and the internal state of the connection
     * @param src
     *  fermion ONV to excite from
     * @param dst
     *  fermion ONV to excite to
     */
    void apply(const FrmOnvField &src, FrmOnvField &dst) const;
    /**
     * update only the common indices, and compute the phase in the process
     * @param src
     *  fermion ONV to excite from
     * @param com
     *  string of "common" spin orbital indices
     * @return
     *  antisymmetric phase of the connection
     */
    bool apply(const FrmOnvField &src, FrmOps &com) const;
    /**
     * update both the dst determinant and compute the associated phase given the src determinant and the internal
     * state of the connection
     * @param src
     *  fermion ONV to excite from
     * @param dst
     *  fermion ONV to excite to
     * @param com
     *  string of "common" spin orbital indices
     * @return
     *  antisymmetric phase of the connection
     */
    bool apply(const FrmOnvField &src, FrmOnvField &dst, FrmOps &com) const;
    /**
     * reset the internal state that to of a null excitation
     */
    void clear();
    /**
     * convenient wrapper for a fermion number-conserving single excitation
     * @param ann
     *  spin orbital index to annihilate from a src determinant
     * @param cre
     *  spin orbital index to create in a dst determinant
     */
    void add(const size_t &ann, const size_t &cre);
    /**
     * convenient wrapper for a fermion number-conserving double excitation
     * @param ann1
     *  first spin orbital index to annihilate from a src determinant
     * @param ann2
     *  second spin orbital index to annihilate from a src determinant
     * @param cre1
     *  first spin orbital index to create in a dst determinant
     * @param cre2
     *  second spin orbital index to create in a dst determinant
     */
    void add(const size_t &ann1, const size_t &ann2, const size_t &cre1, const size_t &cre2);

    void set(const size_t &ann, const size_t &cre) {
        clear();
        add(ann, cre);
    }
    void set(const size_t &ann1, const size_t &ann2, const size_t &cre1, const size_t &cre2){
        clear();
        add(ann1, ann2, cre1, cre2);
    }

    /**
     * convenient wrapper for a fermion number-conserving single excitation
     *
     * The pairs are denoting {spin index, site index} to directly index spin orbitals
     * @param ann
     *  first spin orbital index pair to annihilate from a src determinant
     * @param cre
     *  first spin orbital index pair to create in a dst determinant
     */
    void add(
        std::pair<size_t, size_t> ann,
        std::pair<size_t, size_t> cre
    ) {
        m_ann.add(ann);
        m_cre.add(cre);
    };


    void set(
        std::pair<size_t, size_t> ann,
        std::pair<size_t, size_t> cre
    ) {
        clear();
        add(ann, cre);
    };

    /**
     * convenient wrapper for a fermion number-conserving double excitation
     *
     * The pairs are denoting {spin index, site index} to directly index spin orbitals
     * @param ann1
     *  first spin orbital index pair to annihilate from a src determinant
     * @param ann2
     *  second spin orbital index pair to annihilate from a src determinant
     * @param cre1
     *  first spin orbital index pair to create in a dst determinant
     * @param cre2
     *  second spin orbital index pair to create in a dst determinant
     */
    void add(
        std::pair<size_t, size_t> ann1,
        std::pair<size_t, size_t> ann2,
        std::pair<size_t, size_t> cre1,
        std::pair<size_t, size_t> cre2
    ) {
        m_ann.add(ann1, ann2);
        m_cre.add(cre1, cre2);
    };


    void set(
        std::pair<size_t, size_t> ann1,
        std::pair<size_t, size_t> ann2,
        std::pair<size_t, size_t> cre1,
        std::pair<size_t, size_t> cre2
    ) {
        clear();
        add(ann1, ann2, cre1, cre2);
    };

    /**
     * @return
     *  the annihilation string cast to a vector
     */
    const defs::inds& ann() const;
    /**
     * @return
     *  the creation string cast to a vector
     */
    const defs::inds& cre() const;

    /**
     * @return
     *  total number of operators in connection product
     */
    size_t size() const {
        return m_cre.size()+m_ann.size();
    }

    bool kramers_conserve() const {
        if (m_cre.size()!=m_ann.size()) return false;
        return m_cre.nalpha() == m_ann.nalpha();
    }

private:
    /**
     * update the stored cache with the oddness or evenness of the number of set bits before each dataword in the
     * integer representation of the fermion ONV
     * @param src
     *  fermion ONV to excite from
     */
    void update_dataword_phases(const FrmOnvField &src) const;
    /**
     * compute the antisymmetric phase resulting if a second quantized operator with spin orbital index ibit were
     * operated on the fermion ONV src
     * @param src
     *  fermion ONV to excite from
     * @param ibit
     *  spin orbital index
     * @return
     *  antisymmetric phase of the projection of a single creation or annihilation operator with index ibit
     */
    bool independent_phase(const FrmOnvField &src, const size_t &ibit) const;

public:
    /**
     * the phase computed is the sign arising from the projection of the normal-ordered SQ operator product formed by
     * the product of creation operators a^+(cre_i) annihilation operators a(ann_i) acting on the ket determinant src
     *
     * both strings are applied in their stored order (ascending), but the creation string is hermitian-conjugated so
     * it algebraically should be applied in descending order. This method takes this rearrangement into account
     *
     * phase products are realised via the XOR operator on bools, since like-phases combine for positive (false) phase,
     * and differing phases combine for a negative (true) phase.
     *
     * the algorithm initially sets ann and cre iterators to the cbegin of their respective FrmOpProducts and loops
     * until the iterators have been incremented beyond the last added (highest) element of each string. If the next
     * lowest index is a creation operator, the phase associated with moving the operator past the remaining
     * annihilation operators is taken into account.
     *
     * at this point, the determinant and connection index-specific antisymmety has been dealt with. All that remains is
     * the include the excitation rank-associated details. The number of fermion exchanges required to invert a product
     * of distinct operators of length n is the number of pairs that can be drawn from n objects. In this case, there
     * is a phase associated with applying the creation operators in reverse order, and another associated with the fact
     * that each applied operator had to be moved through each previously handled one as it was being (virtually)
     * inserted into or removed from the determinant.
     *
     * Let m be the number of creation ops, and n the number of annihilation ops, then the total number of exchanges
     * entailed in the above process is the sum (or, equivalently the difference) of the two numbers of pairs:
     *
     * nexchange = npair(m+n) - npair(n)
     * i.e. 2*nexchange = (m+n)*(m+n-1) - n^2 + n = m^2 - m + 2mn
     * i.e. nexchange = npair(m) + mn
     *
     * The numbers of pairs alternate between odd and even with a period of 2:
     * 0->0   1->0   2->1   3->3   4->6   5->10   6->15
     *   E      E      O      O      E      E       O    ....
     *
     * so it is quick to compute whether these exchanges implies an inversion in the phase computed from the connection
     * indices, and this is done at the end of the implementation of this method
     *
     * @param src
     *  fermion ONV to excite from
     * @return
     *  true if the phase of the connection with respect to the src is negative
     */
    bool phase(const FrmOnvField &src) const;

    size_t exsig() const;

    /**
     * @param nop_insert
     *  number of pairs of same-indexed creation and annihilation operators to insert into the normal-ordered product
     * @return
     *  excitation signature of the resulting term
     */
    size_t exsig(const size_t& nop_insert) const;
};


#endif //M7_FRMONVCONNECTION_H