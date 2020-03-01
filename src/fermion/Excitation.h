//
// Created by Robert John Anderson on 2020-02-23.
//

#ifndef M7_EXCITATION_H
#define M7_EXCITATION_H


#include <src/defs.h>
#include "Determinant.h"

struct Excitation {
    const Determinant &m_det;
    const defs::inds m_removed;
    const defs::inds m_inserted;

    const defs::prob_t m_prob;
    const defs::ham_t m_helement;

    Excitation(const Determinant &det, const defs::inds &removed, const defs::inds &inserted,
               const defs::prob_t &prob, const defs::ham_t &helement);

    // null excitation
    Excitation(const Determinant &det);

    // single
    Excitation(const Determinant &det, const size_t &removed, const size_t &inserted,
               const defs::prob_t &prob, const defs::ham_t &helement);

    //double
    Excitation(const Determinant &det, const size_t &removed_1, const size_t &removed_2,
               const size_t &inserted_1, const size_t &inserted_2,
               const defs::prob_t &prob, const defs::ham_t &helement);

    bool is_null() const;

    Determinant get_connection() const {
        return m_det.get_excited_det(m_removed, m_inserted);
    }
};


#endif //M7_EXCITATION_H
