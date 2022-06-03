//
// Created by Robert J. Anderson on 26/07/2021.
//

#ifndef M7_BOSONVCONNECTION_H
#define M7_BOSONVCONNECTION_H

#include <M7_lib/field/BosOnvField.h>
#include <M7_lib/parallel/MPIAssert.h>

struct BosOpPair {
    const size_t m_imode;
    size_t m_nop;
    BosOpPair(size_t imode, size_t nop);
};

class BosOps {
    std::vector<BosOpPair> m_pairs;
    /**
     * vector of length nmode which enables constant-time access to a pair
     */
    std::vector<BosOpPair*> m_pair_ptrs;
    size_t m_nop = 0ul;
public:
    BosOps(size_t nmode);

    const std::vector<BosOpPair>& pairs() const;

    /**
     * get vector of mode indices (multiple nop pairs will be duplicated)
     * @return
     *  vector of mode indices
     */
    defs::inds get() const;

    /**
     * set from vector of mode indices (adjacent like indices will be rolled together in the pair vector).
     * @param imodes
     *  vector of mode indices
     */
    void set(const defs::inds& imodes);

    /**
     * set from mode indices (adjacent like indices will be rolled together in the pair vector).
     * does the same as set but for specific number of modes without dynamic allocation of vector container,
     * which will be more efficient
     * @param i
     *  mode index
     */
    void set(size_t i);

    /**
     * set from mode indices (adjacent like indices will be rolled together in the pair vector).
     * @param i
     *  mode index (lowest)
     * @param j
     *  mode index (highest)
     */
    void set(size_t i, size_t j);

    /**
     * set from mode indices (adjacent like indices will be rolled together in the pair vector)
     * @param i
     *  mode index (lowest)
     * @param j
     *  mode index
     * @param k
     *  mode index (highest)
     */
    void set(size_t i, size_t j, size_t k);

    void add(BosOpPair&& pair);

    void add(size_t imode, size_t nop);

    /**
     * add the given number of operators
     * @param imode
     *  mode index which is taken to be already occupied in the connection: error if not
     * @param nop
     *  number of operators to add
     * @return
     *  new total nop of imode
     */
    size_t add_to_nonempty(size_t imode, size_t nop=1ul) {
        auto ptr = m_pair_ptrs[imode];
        DEBUG_ASSERT_TRUE(ptr, "mode is not already part of the boson operator product");
        m_nop+=nop;
        return ptr->m_nop+=nop;
    }

    size_t add_to(size_t imode, size_t nop=1ul) {
        auto ptr = m_pair_ptrs[imode];
        if (ptr) {
            m_nop+=nop;
            return ptr->m_nop+=nop;
        }
        add(imode, nop);
        return nop;
    }

    void clear();

    size_t size() const;

    const BosOpPair& operator[](const size_t& ipair) const;

    /**
     * @param iop
     *  operator index in the expanded product
     * @return
     *  mode index corresponding to the operator index
     */
    size_t get_imode(size_t iop) const;
};

struct BosOnvConnection {
    BosOps m_ann, m_cre;

    explicit BosOnvConnection(size_t nmode);

    explicit BosOnvConnection(sys::Size extents);

    explicit BosOnvConnection(const BosOnvField& mbf);

    void clear();

    size_t size() const;

    size_t nmode() const;

    void connect(const BosOnvField& src, const BosOnvField& dst);

    void connect(const BosOnvField& src, const BosOnvField& dst, BosOps& com);

    void apply(const BosOnvField &src, BosOnvField &dst) const;

    void apply(const BosOnvField &src, BosOps &com) const;

    void apply(const BosOnvField &src, BosOnvField &dst, BosOps &com) const;

    size_t exsig() const;

    bool respects_occ_range(const BosOnvField &src, size_t nboson_max) const;

};


#endif //M7_BOSONVCONNECTION_H
