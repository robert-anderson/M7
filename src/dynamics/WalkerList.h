//
// Created by Robert John Anderson on 2020-02-24.
//

#ifndef M7_WALKERLIST_H
#define M7_WALKERLIST_H

#include "src/data/PerforableMappedList.h"
#include <memory>

struct WalkerListFields {
    Specification m_spec;
    size_t idet;
    size_t iweight;
    size_t ihdiag;
    size_t iflag_reference_connection;
    size_t iflag_initiator;
    size_t iflag_deterministic;

    WalkerListFields(size_t nsite) {
        idet = m_spec.add<Determinant>(nsite);
        iweight = m_spec.add<defs::ham_t>(1);
        ihdiag = m_spec.add<defs::ham_comp_t>(1);
        iflag_reference_connection = m_spec.add<bool>(1);
        iflag_initiator = m_spec.add<bool>(1);
        iflag_deterministic = m_spec.add<bool>(1);
    }
};

class WalkerList : public PerforableMappedList<Determinant> {
    WalkerListFields m_fields;
public:
    WalkerList(size_t nsite, size_t nrow);

    Determinant get_determinant(const size_t &irow);

    NumericView<defs::ham_t> get_weight(const size_t &irow);

    NumericView<defs::ham_comp_t> get_hdiag(const size_t &irow);

    NumericView<bool> get_flag_reference_connection(const size_t &irow);

    NumericView<bool> get_flag_initiator(const size_t &irow);

    NumericView<bool> get_flag_deterministic(const size_t &irow);
};


#endif //M7_WALKERLIST_H
