//
// Created by RJA on 11/09/2020.
//

#ifndef M7_BOSONONVCONNECTION_H
#define M7_BOSONONVCONNECTION_H

#include "src/core/field/Fields.h"
#include "src/core/field/Views.h"


class BosonOnvConnection {
    struct Diff {
        defs::inds m_changed_modes;
        size_t m_nchanged_mode = 0ul;
        std::vector<int> m_changes;
        Diff(size_t nmode);
        void zero();
    };
    const size_t m_nmode;
    std::vector<int> m_com;
    Diff m_diff;

public:
    explicit BosonOnvConnection(const BosonOnvSpecifier& spec);

    const size_t & nchanged_mode() const;
    const size_t & changed_mode(const size_t& ichange) const;

    const int & changes(const size_t& ichange) const;

    const int & com(const size_t& icom) const;

    BosonOnvConnection(const views::BosonOnv &in, const views::BosonOnv &out);

    explicit BosonOnvConnection(const views::BosonOnv &in);

    void connect(const views::BosonOnv &in, const views::BosonOnv &out);

    void apply(const views::BosonOnv &in, views::BosonOnv& out);

    void add(const size_t imode, const int change);

    void zero(){
        m_diff.zero();
    }
};


#endif //M7_BOSONONVCONNECTION_H
