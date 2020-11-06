//
// Created by RJA on 11/09/2020.
//

#include "BosonOnvConnection.h"

PermanentDiff::PermanentDiff(size_t nmode) : m_changed_modes(nmode, 0ul), m_changes(nmode, 0){}

void PermanentDiff::zero() {
    m_nchanged_mode = 0ul;
}

BosonOnvConnection::BosonOnvConnection(const BosonOnvSpecifier &spec) :
        m_nmode(spec.nmode()), m_com(m_nmode), m_diff(m_nmode) {}

const size_t &BosonOnvConnection::nchanged_mode() const {
    return m_diff.m_nchanged_mode;
}

const size_t &BosonOnvConnection::changed_mode(const size_t &ichange) const {
    ASSERT(ichange < nchanged_mode());
    return m_diff.m_changed_modes[ichange];
}

const int &BosonOnvConnection::changes(const size_t &ichange) const {
    ASSERT(ichange < nchanged_mode());
    return m_diff.m_changes[ichange];
}

const int &BosonOnvConnection::com(const size_t &icom) const {
    ASSERT(icom < m_nmode);
    return m_com[icom];
}

BosonOnvConnection::BosonOnvConnection(const views::BosonOnv &in, const views::BosonOnv &out) :
        BosonOnvConnection(in.spec()){
    connect(in, out);
}

BosonOnvConnection::BosonOnvConnection(const views::BosonOnv &in) :
        BosonOnvConnection(in.spec()){
    connect(in, in);
}

void BosonOnvConnection::connect(const views::BosonOnv &in, const views::BosonOnv &out) {
    m_diff.zero();
    ASSERT(!m_com.empty())
    for(size_t imode=0ul; imode<m_nmode; ++imode){
        int nin = in(imode);
        int nout = out(imode);
        m_com[imode] = std::min(nin, nout);
        if (nin!=nout){
            m_diff.m_changed_modes[m_diff.m_nchanged_mode] = imode;
            m_diff.m_changes[m_diff.m_nchanged_mode] = nin-nout;
            ++m_diff.m_nchanged_mode;
        }
    }
}

void BosonOnvConnection::apply(const views::BosonOnv &in, views::BosonOnv &out) {
    out = in;
    for (size_t ichange = 0ul; ichange<nchanged_mode(); ++ichange){
        out(changed_mode(ichange)) += changes(ichange);
    }
}
