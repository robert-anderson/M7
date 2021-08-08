//
// Created by rja on 26/07/2021.
//

#ifndef M7_BOSONVCONNECTION_H
#define M7_BOSONVCONNECTION_H

#include <src/core/field/Fields.h>
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

    void clear();

    size_t size() const;

    void add(BosOpPair&& pair);

    void set(BosOpPair&& pair);

    const BosOpPair& operator[](const size_t& ipair) const;
};

struct BosOnvConnection {
    BosOps m_ann, m_cre;
    BosOnvConnection(size_t nmode);

    void clear();

    size_t size() const;

    void connect(const fields::BosOnv& src, const fields::BosOnv& dst);

    void connect(const fields::BosOnv& src, const fields::BosOnv& dst, BosOps& com);

    void apply(const fields::BosOnv &src, fields::BosOnv &dst) const;

    void apply(const fields::BosOnv &src, BosOps &com) const;

    void apply(const fields::BosOnv &src, fields::BosOnv &dst, BosOps &com) const;
};


#endif //M7_BOSONVCONNECTION_H
