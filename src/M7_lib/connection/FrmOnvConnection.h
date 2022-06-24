//
// Created by Robert J. Anderson on 23/07/2021.
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
    const sys::frm::Size m_sites;
    defs::uintv_t m_inds;
public:
    FrmOps(size_t nsite): m_sites(nsite) {
        m_inds.reserve(2*m_sites);
    }

    FrmOps(const FrmOps& other): FrmOps(other.m_sites){}

    FrmOps& operator=(const FrmOps& other){
        m_inds = other.m_inds;
        return *this;
    }

    FrmOps& operator=(const defs::uintv_t& inds){
        m_inds = inds;
        return *this;
    }

    const defs::uintv_t & inds() const {
        return m_inds;
    }

    bool operator==(const defs::uintv_t &v) const {
        return m_inds==v;
    }

    const size_t &operator[](size_t i) const {
        DEBUG_ASSERT_LT(i, m_inds.size(), "operator index is OOB");
        return m_inds[i];
    }

    void clear() {
        m_inds.clear();
    }

    size_t size() const {
        return m_inds.size();
    }

    void add(size_t ispinorb) {
        DEBUG_ASSERT_TRUE(m_inds.empty() || ispinorb > m_inds.back(),
                          "spin orbital indices should be added in ascending order");
        DEBUG_ASSERT_LT(ispinorb, m_sites.m_nspinorb, "spin orbital index OOB");
        m_inds.push_back(ispinorb);
    }

    void add(std::pair<size_t, size_t> pair) {
        add(m_sites.ispinorb(pair));
    }

    /**
     * set a selection of spin orbital indices
     * @param orbs
     *  spin orbital indices
     * @param inds
     *  positions in the orbs vector of the selected spin orbitals
     */
    void set_selection(const defs::uintv_t& orbs, const defs::uintv_t& inds) {
        clear();
        for (const auto& ind: inds) {
            DEBUG_ASSERT_LT(ind, orbs.size(), "index OOB");
            add(orbs[ind]);
        }
        DEBUG_ASSERT_EQ(inds.size(), size(), "not all selected uintv_t were added");
    }

    void set(size_t ispinorb){
        clear();
        add(ispinorb);
    }

    void set(size_t ispinorb1, size_t ispinorb2){
        clear();
        add(ispinorb1);
        add(ispinorb2);
    }

    void set_in_order(size_t ispinorb1, size_t ispinorb2){
        clear();
        DEBUG_ASSERT_NE(ispinorb1, ispinorb2, "fermion indices in connection product may not coincide");
        if (ispinorb1 < ispinorb2) {
            add(ispinorb1);
            add(ispinorb2);
        }
        else {
            add(ispinorb2);
            add(ispinorb1);
        }
    }

    void set(std::pair<size_t, size_t> pair){
        clear();
        add(pair);
    }

    void set(std::pair<size_t, size_t> pair1, std::pair<size_t, size_t> pair2){
        clear();
        add(pair1);
        add(pair2);
    }

    void set_in_order(std::pair<size_t, size_t> pair1, std::pair<size_t, size_t> pair2){
        set_in_order(m_sites.ispinorb(pair1), m_sites.ispinorb(pair2));
    }

    defs::uintv_t::const_iterator cbegin() const {
        return m_inds.cbegin();
    }

    defs::uintv_t::const_iterator cend() const {
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
        for (auto i: m_inds) n+=!m_sites.ispin(i);
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

    explicit FrmOnvConnection(const sys::frm::Size& size);

    /*
     * universal interface (works with any compile-time value "MBF_TYPE")
     */
    explicit FrmOnvConnection(sys::Size size);

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
     * @return
     *  the annihilation string cast to a vector
     */
    const defs::uintv_t& ann() const;
    /**
     * @return
     *  the creation string cast to a vector
     */
    const defs::uintv_t& cre() const;

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
     * nexchange = npair(m+n) - nmode(n)
     * i.e. 2*nexchange = (m+n)*(m+n-1) - n^2 + n = m^2 - m + 2mn
     * i.e. nexchange = nmode(m) + mn
     *
     * The numbers of pairs alternate between odd and even with a period of 2:
     * 0->0   1->0   2->1   3->3   4->6   5->10   6->15
     *   E      E      O      O      E      E       O    ....
     *
     * so it is quick to compute whether these exchanges imply an inversion in the phase computed from the connection
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

static std::ostream &operator<<(std::ostream &os, const FrmOnvConnection &conn) {
    os << conn.m_ann.inds() << "->" << conn.m_cre.inds();
    return os;
}


#endif //M7_FRMONVCONNECTION_H
