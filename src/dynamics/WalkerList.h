//
// Created by Robert John Anderson on 2020-02-24.
//

#ifndef M7_WALKERLIST_H
#define M7_WALKERLIST_H

#include "src/data/PerforableMappedList.h"
#include <memory>

class WalkerList {
public:
    std::unique_ptr<PerforableMappedList<Determinant>> m_list = nullptr;
    Specification m_spec{};

    size_t m_idet;
    size_t m_iweight;
    size_t m_ihdiag;
    size_t m_iflag_reference_connection;
    size_t m_iflag_initiator;

    WalkerList(const size_t &nspinorb, const size_t &nrow) {
        m_idet = m_spec.add<Determinant>(nspinorb);
        m_iweight = m_spec.add<defs::ham_t>(1);
        m_ihdiag = m_spec.add<defs::ham_comp_t>(1);
        m_iflag_reference_connection = m_spec.add<bool>(1);
        m_iflag_initiator = m_spec.add<bool>(1);
        m_list = std::make_unique<PerforableMappedList<Determinant>>(m_spec, nrow, m_idet);
    }

    size_t high_water_mark(){
        return m_list->high_water_mark();
    }

    Determinant det(const size_t &irow){
        return m_list->view<Determinant>(irow, m_idet);
    }

    bool is_free(const size_t &irow){
        return det().is_zero();
    }

    defs::ham_t weight(const size_t &irow){
        return *m_list->view<defs::ham_t>(irow, m_iweight);
    }

};


#endif //M7_WALKERLIST_H
