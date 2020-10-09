//
// Created by RJA on 11/09/2020.
//

#ifndef M7_PERMANENTCONNECTION_H
#define M7_PERMANENTCONNECTION_H


#include "PermanentField.h"

struct PermanentDiff {
    defs::inds m_changed_modes;
    size_t m_nchanged_mode = 0ul;
    std::vector<int> m_changes;
    PermanentDiff(size_t nmode) : m_changed_modes(nmode, 0ul), m_changes(nmode, 0){}

    void zero(){
        m_nchanged_mode = 0ul;
    }
};

class PermanentConnection {
    const size_t m_nmode;
    std::vector<int> m_com;
    PermanentDiff m_diff;

public:
    explicit PermanentConnection(const PermanentField* field):
    m_nmode(field->nelement()), m_com(m_nmode), m_diff(m_nmode) {
    }

    const size_t & nchanged_mode() const {
        return m_diff.m_nchanged_mode;
    }
    const size_t & changed_modes(const size_t& ichange) const {
        ASSERT(ichange < nchanged_mode());
        return m_diff.m_changed_modes[ichange];
    }
    // todo: uncomment when current usage has been changed where appropriate.
//    const std::vector<size_t>& changed_modes() const {
//        return m_diff.m_changed_modes;
//    }
    const int & changes(const size_t& ichange) const {
        ASSERT(ichange < nchanged_mode());
        return m_diff.m_changes[ichange];
    }
//    const std::vector<int>& changes() const {
//        return m_diff.m_changes;
//    }
    const int & com(const size_t& icom) const {
        ASSERT(icom < m_nmode - nchanged_mode());
        return m_com[icom];
    }
//    const std::vector<int>& com() const {
//        ASSERT(icom < m_nmode - nchanged_mode());
//        return m_com;
//    }

    PermanentConnection(const PermanentElement &ket, const PermanentElement &bra):
            PermanentConnection(ket.field()){
        connect(ket, bra);
    }

    explicit PermanentConnection(const PermanentElement &ket):
            PermanentConnection(ket.field()){
        connect(ket, ket);
    }

    void connect(const PermanentElement &ket, const PermanentElement &bra){
        m_diff.zero();
        ASSERT(!m_com.empty())
        for(size_t imode=0ul; imode<m_nmode; ++imode){
            int nket = *ket(imode);
            int nbra = *bra(imode);
            m_com[imode] = std::min(nket, nbra);
            if (nket!=nbra){
                m_diff.m_changed_modes[m_diff.m_nchanged_mode] = imode;
                m_diff.m_changes[m_diff.m_nchanged_mode] = nket-nbra;
                ++m_diff.m_nchanged_mode;
            }
        }
    }
};


#endif //M7_PERMANENTCONNECTION_H
