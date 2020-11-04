//
// Created by RJA on 11/09/2020.
//

#ifndef M7_BOSONICCONNECTION_H
#define M7_BOSONICCONNECTION_H

#include "src/core/field/Fields.h"
#include "src/core/field/Views.h"

struct PermanentDiff {
    defs::inds m_changed_modes;
    size_t m_nchanged_mode = 0ul;
    std::vector<int> m_changes;
    PermanentDiff(size_t nmode) : m_changed_modes(nmode, 0ul), m_changes(nmode, 0){}

    void zero(){
        m_nchanged_mode = 0ul;
    }
};

class BosonicConnection {
    const size_t m_nmode;
    std::vector<int> m_com;
    PermanentDiff m_diff;

public:
    explicit BosonicConnection(const BosonOnvSpecifier& spec):
    m_nmode(spec.nmode()), m_com(m_nmode), m_diff(m_nmode) {}

    const size_t & nchanged_mode() const {
        return m_diff.m_nchanged_mode;
    }
    const size_t & changed_mode(const size_t& ichange) const {
        ASSERT(ichange < nchanged_mode());
        return m_diff.m_changed_modes[ichange];
    }

    const int & changes(const size_t& ichange) const {
        ASSERT(ichange < nchanged_mode());
        return m_diff.m_changes[ichange];
    }

    const int & com(const size_t& icom) const {
        ASSERT(icom < m_nmode - nchanged_mode());
        return m_com[icom];
    }

    BosonicConnection(const views::BosonOnv &ket, const views::BosonOnv &bra):
            BosonicConnection(ket.spec()){
        connect(ket, bra);
    }

    explicit BosonicConnection(const views::BosonOnv &ket):
            BosonicConnection(ket.spec()){
        connect(ket, ket);
    }

    void connect(const views::BosonOnv &ket, const views::BosonOnv &bra){
        m_diff.zero();
        ASSERT(!m_com.empty())
        for(size_t imode=0ul; imode<m_nmode; ++imode){
            int nket = ket(imode);
            int nbra = bra(imode);
            m_com[imode] = std::min(nket, nbra);
            if (nket!=nbra){
                m_diff.m_changed_modes[m_diff.m_nchanged_mode] = imode;
                m_diff.m_changes[m_diff.m_nchanged_mode] = nket-nbra;
                ++m_diff.m_nchanged_mode;
            }
        }
    }
};


#endif //M7_BOSONICCONNECTION_H
