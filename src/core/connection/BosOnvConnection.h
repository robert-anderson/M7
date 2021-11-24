//
// Created by rja on 26/07/2021.
//

#ifndef M7_BOSONVCONNECTION_H
#define M7_BOSONVCONNECTION_H

#include "src/core/field/BosOnvField.h"
#include "src/core/parallel/MPIAssert.h"

struct BosOpPair {
    const size_t m_imode;
    const size_t m_nop;
    BosOpPair(size_t imode, size_t nop);
};

class BosOps {
    std::vector<BosOpPair> m_pairs;
    size_t m_nop = 0ul;
public:
    BosOps(size_t nmode);

    const std::vector<BosOpPair>& pairs() const;

    defs::inds to_vector() const;

    void clear();

    size_t size() const;

    void add(BosOpPair&& pair);

    void set(BosOpPair&& pair);

    void set(const size_t& imode);

    void set(const size_t& imode, const size_t& jmode);

    const BosOpPair& operator[](const size_t& ipair) const;
};

struct BosOnvConnection {
    BosOps m_ann, m_cre;

    BosOnvConnection(size_t nmode);

    BosOnvConnection(BasisDims bd);

    BosOnvConnection(const BosOnvField& mbf);

    void clear();

    size_t size() const;

    void connect(const BosOnvField& src, const BosOnvField& dst);

    void connect(const BosOnvField& src, const BosOnvField& dst, BosOps& com);

    void apply(const BosOnvField &src, BosOnvField &dst) const;

    void apply(const BosOnvField &src, BosOps &com) const;

    void apply(const BosOnvField &src, BosOnvField &dst, BosOps &com) const;

    size_t exsig() const;
};


#endif //M7_BOSONVCONNECTION_H
