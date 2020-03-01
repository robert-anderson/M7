//
// Created by Robert John Anderson on 2020-02-23.
//

#include <assert.h>
#include "Excitation.h"

Excitation::Excitation(const Determinant &det, const defs::inds &removed, const defs::inds &inserted,
                       const defs::prob_t &prob,
                       const defs::ham_t &helement) :
        m_det(det), m_removed(removed), m_inserted(inserted), m_helement(helement), m_prob(prob){
    assert(prob>0.0);
}

Excitation::Excitation(const Determinant &det) : m_det(det), m_helement(0.0), m_prob(0.0) {}

Excitation::Excitation(const Determinant &det, const size_t &removed, const size_t &inserted, const defs::prob_t &prob,
                       const defs::ham_t &helement) :
    Excitation(det, defs::inds{removed}, defs::inds{inserted}, prob, helement) {}

Excitation::Excitation(const Determinant &det, const size_t &removed_1, const size_t &removed_2,
                       const size_t &inserted_1,
                       const size_t &inserted_2, const defs::prob_t &prob, const defs::ham_t &helement) :
    Excitation(det, defs::inds{removed_1, removed_2}, defs::inds{inserted_1, inserted_2}, prob, helement) {}

bool Excitation::is_null() const {
    return m_removed.empty() && m_inserted.empty();
}
