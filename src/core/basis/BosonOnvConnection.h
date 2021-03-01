//
// Created by RJA on 11/09/2020.
//

#ifndef M7_BOSONONVCONNECTION_H
#define M7_BOSONONVCONNECTION_H

#include "src/core/field/Fields.h"

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

    explicit BosonOnvConnection(const size_t& nmode);
    BosonOnvConnection(const fields::BosonOnv &in, const fields::BosonOnv &out);
    explicit BosonOnvConnection(const fields::BosonOnv &in);

    operator bool() const {return nchanged_mode();}

    const size_t & nchanged_mode() const;
    const size_t & changed_mode(const size_t& ichange) const;

    const int & changes(const size_t& ichange) const;

    const int & com(const size_t& icom) const;

    void connect(const fields::BosonOnv &in, const fields::BosonOnv &out);

    void apply(const fields::BosonOnv &in, fields::BosonOnv& out);

    void add(const size_t& imode, const int& change);

    size_t nexcit() const {
        size_t res = 0ul;
        for (size_t ichange=0ul; ichange<nchanged_mode(); ++ichange){
            res+=std::abs(changes(ichange));
        }
        return res;
    }

    void zero(){
        m_diff.zero();
    }
};


#endif //M7_BOSONONVCONNECTION_H
