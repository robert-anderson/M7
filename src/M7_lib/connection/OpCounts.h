//
// Created by rja on 26/03/23.
//

#ifndef M7_OPCOUNTS_H
#define M7_OPCOUNTS_H

#include "OpSig.h"
#include "M7_lib/field/Fields.h"


/**
 * encodes a pair of creation and annihilation operator counts with unlimited range, unlike OpSig which has only a few
 * bits for each count
 */
struct OpCounts {
    uint_t m_nfrm_cre = 0ul;
    uint_t m_nfrm_ann = 0ul;
    uint_t m_nbos_cre = 0ul;
    uint_t m_nbos_ann = 0ul;

private:
    void set(const field::FrmOnv &src, const field::FrmOnv &dst);

    void set(const field::BosOnv &src, const field::BosOnv &dst);

    void set(const field::FrmBosOnv &src, const field::FrmBosOnv &dst);
public:

    template<typename mbf_t>
    OpCounts(const mbf_t& src, const mbf_t& dst) {
        set(src, dst);
    }

    OpSig opsig() const {
        return {{m_nfrm_cre, m_nfrm_ann}, {m_nbos_ann, m_nbos_cre}};
    }
};

#endif //M7_OPCOUNTS_H
