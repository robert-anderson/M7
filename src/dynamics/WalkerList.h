//
// Created by Robert John Anderson on 2020-02-24.
//

#ifndef M7_WALKERLIST_H
#define M7_WALKERLIST_H

#include "src/data/PerforableMappedList.h"
#include <memory>

struct WalkerListSpecification : public Specification {
    size_t idet;
    size_t iweight;
    size_t ihdiag;
    size_t iflag_reference_connection;
    size_t iflag_initiator;
    size_t iflag_deterministic;
    WalkerListSpecification(size_t nsite);
};

class WalkerList : public PerforableMappedList<Determinant> {
public:
    using spec_T = WalkerListSpecification;
protected:
    const spec_T m_spec;
public:
    WalkerList(const spec_T &spec, size_t nrow);

    Determinant get_determinant(const size_t &irow);

    NumericView<defs::ham_t> get_weight(const size_t &irow);

    NumericView<defs::ham_comp_t> get_hdiag(const size_t &irow);

    NumericView<bool> get_flag_reference_connection(const size_t &irow);

    NumericView<bool> get_flag_initiator(const size_t &irow);

    NumericView<bool> get_flag_deterministic(const size_t &irow);
};


#endif //M7_WALKERLIST_H
