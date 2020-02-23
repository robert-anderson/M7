//
// Created by Robert John Anderson on 2020-02-23.
//

#ifndef M7_EXCITATION_H
#define M7_EXCITATION_H


#include <src/defs.h>
#include "Determinant.h"

class Excitation {
    const Determinant &m_det;
    defs::inds m_removed;
    defs::inds m_inserted;
    /*
     * m_factor is the ratio of the H element to the probability of
     * proposing the excitation
     *
     * the update to the weight of the connected determinant due to this
     * excitation is -tau * m_factor * (weight of m_det) / (number of spawning attempts)
     */
    defs::ham_t m_factor;
public:
    Excitation(const Determinant &det, const defs::inds &removed,
        const defs::inds &inserted, const defs::ham_t &factor);

    // null excitation
    Excitation(const Determinant &det);

    // single
    Excitation(const Determinant &det, const size_t &removed,
        const size_t &inserted, const defs::ham_t &factor);

    //double
    Excitation(const Determinant &det,
               const size_t &removed_1, const size_t &removed_2,
               const size_t &inserted_1, const size_t &inserted_2, const defs::ham_t &factor);

    bool is_null() const;

    Determinant get_connection() const {
    }
};


#endif //M7_EXCITATION_H
