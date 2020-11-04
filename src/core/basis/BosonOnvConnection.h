//
// Created by RJA on 11/09/2020.
//

#ifndef M7_BOSONONVCONNECTION_H
#define M7_BOSONONVCONNECTION_H

#include "src/core/field/Fields.h"
#include "src/core/field/Views.h"

struct PermanentDiff {
    defs::inds m_changed_modes;
    size_t m_nchanged_mode = 0ul;
    std::vector<int> m_changes;
    PermanentDiff(size_t nmode);

    void zero();
};

class BosonOnvConnection {
    const size_t m_nmode;
    std::vector<int> m_com;
    PermanentDiff m_diff;

public:
    explicit BosonOnvConnection(const BosonOnvSpecifier& spec);

    const size_t & nchanged_mode() const;
    const size_t & changed_mode(const size_t& ichange) const;

    const int & changes(const size_t& ichange) const;

    const int & com(const size_t& icom) const;

    BosonOnvConnection(const views::BosonOnv &ket, const views::BosonOnv &bra);

    explicit BosonOnvConnection(const views::BosonOnv &ket);

    void connect(const views::BosonOnv &ket, const views::BosonOnv &bra);

    void apply(const views::BosonOnv &ket, views::BosonOnv& bra);
};


#endif //M7_BOSONONVCONNECTION_H
