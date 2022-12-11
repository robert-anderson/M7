//
// Created by Robert J. Anderson on 26/07/2021.
//

#ifndef M7_BOSONVCONNECTION_H
#define M7_BOSONVCONNECTION_H

#include <M7_lib/field/BosOnvField.h>
#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/connection/OpSig.h>

struct BosOpPair {
    const uint_t m_imode;
    uint_t m_nop;
    BosOpPair(uint_t imode, uint_t nop);
};

class BosOps {
    v_t<BosOpPair> m_pairs;
    /**
     * vector of length nmode which enables constant-time access to a pair
     */
    v_t<BosOpPair*> m_pair_ptrs;
    uint_t m_nop = 0ul;
public:
    BosOps(uint_t nmode);

    const v_t<BosOpPair>& pairs() const;

    /**
     * get vector of mode indices (multiple nop pairs will be duplicated)
     * @return
     *  vector of mode indices
     */
    uintv_t get() const;

    /**
     * set from vector of mode indices (adjacent like indices will be rolled together in the pair vector).
     * @param imodes
     *  vector of mode indices
     */
    void set(const uintv_t& imodes);

    /**
     * set from mode indices (adjacent like indices will be rolled together in the pair vector).
     * does the same as set but for specific number of modes without dynamic allocation of vector container,
     * which will be more efficient
     * @param i
     *  mode index
     */
    void set(uint_t i);

    /**
     * set from mode indices (adjacent like indices will be rolled together in the pair vector).
     * @param i
     *  mode index (lowest)
     * @param j
     *  mode index (highest)
     */
    void set(uint_t i, uint_t j);

    /**
     * set from mode indices (adjacent like indices will be rolled together in the pair vector)
     * @param i
     *  mode index (lowest)
     * @param j
     *  mode index
     * @param k
     *  mode index (highest)
     */
    void set(uint_t i, uint_t j, uint_t k);

    void add(uint_t imode, uint_t nop);

    /**
     * add the given number of operators
     * @param imode
     *  mode index which is taken to be already occupied in the connection: error if not
     * @param nop
     *  number of operators to add
     * @return
     *  new total nop of imode
     */
    uint_t add_to_nonempty(uint_t imode, uint_t nop=1ul) {
        auto ptr = m_pair_ptrs[imode];
        DEBUG_ASSERT_TRUE(ptr, "mode is not already part of the boson operator product");
        m_nop+=nop;
        return ptr->m_nop+=nop;
    }

    uint_t add_to(uint_t imode, uint_t nop=1ul) {
        auto ptr = m_pair_ptrs[imode];
        if (ptr) {
            m_nop+=nop;
            return ptr->m_nop+=nop;
        }
        add(imode, nop);
        return nop;
    }

    void clear();

    uint_t size() const;

    const BosOpPair& operator[](const uint_t& ipair) const;

    /**
     * @param iop
     *  operator index in the expanded product
     * @return
     *  mode index corresponding to the operator index
     */
    uint_t get_imode(uint_t iop) const;

    str_t to_string() const;
};

struct BosOnvConnection {
    BosOps m_ann, m_cre;

    explicit BosOnvConnection(uint_t nmode);

    explicit BosOnvConnection(sys::Size extents);

    explicit BosOnvConnection(const BosOnvField& mbf);

    void clear();

    uint_t size() const;

    uint_t nmode() const;

    void connect(const BosOnvField& src, const BosOnvField& dst);

    void connect(const BosOnvField& src, const BosOnvField& dst, BosOps& com);

    void apply(const BosOnvField &src, BosOnvField &dst) const;

    void apply(const BosOnvField &src, BosOps &com) const;

    void apply(const BosOnvField &src, BosOnvField &dst, BosOps &com) const;

    OpSig exsig() const;

    bool respects_occ_range(const BosOnvField &src, uint_t nboson_max) const;

};

static std::ostream &operator<<(std::ostream &os, const BosOps &ops) {
    os << ops.to_string();
    return os;
}

static std::ostream &operator<<(std::ostream &os, const BosOnvConnection &conn) {
    os << conn.m_ann << "->" << conn.m_cre;
    return os;
}

#endif //M7_BOSONVCONNECTION_H
