//
// Created by RJA on 11/09/2020.
//

#include "BosonOnvConnection.h"

BosonOnvConnection::Diff::Diff(size_t nmode) : m_changed_modes(nmode, 0ul), m_changes(nmode, 0) {}

void BosonOnvConnection::Diff::zero() {
    m_nchanged_mode = 0ul;
}

BosonOnvConnection::BosonOnvConnection(const size_t& nmode) :
        m_nmode(nmode), m_com(nmode), m_diff(nmode) {}

BosonOnvConnection::BosonOnvConnection(const fields::BosonOnv &in, const fields::BosonOnv &out) :
        BosonOnvConnection(in.nelement()) {
    connect(in, out);
}

BosonOnvConnection::BosonOnvConnection(const fields::BosonOnv &in) :
        BosonOnvConnection(in.nelement()) {
    connect(in, in);
}

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

void BosonOnvConnection::connect(const fields::BosonOnv &in, const fields::BosonOnv &out) {
    m_diff.zero();
    ASSERT(!m_com.empty())
    for (size_t imode = 0ul; imode < m_nmode; ++imode) {
        int nin = in[imode];
        int nout = out[imode];
        m_com[imode] = std::min(nin, nout);
        if (nin != nout) {
            m_diff.m_changed_modes[m_diff.m_nchanged_mode] = imode;
            m_diff.m_changes[m_diff.m_nchanged_mode] = nout - nin;
            ++m_diff.m_nchanged_mode;
        }
    }
}

void BosonOnvConnection::apply(const fields::BosonOnv &in, fields::BosonOnv &out) {
    out = in;
    for (size_t imode = 0ul; imode < m_nmode; ++imode) {
        m_com[imode] = in[imode];
    }
    for (size_t ichange = 0ul; ichange < nchanged_mode(); ++ichange) {
        const auto imode = changed_mode(ichange);
        const auto change = changes(ichange);
        out[imode] += changes(ichange);
        if (change<0) m_com[imode] += change;
    }
}

void BosonOnvConnection::add(const size_t& imode, const int& change) {
    m_diff.m_changes[m_diff.m_nchanged_mode] = change;
    m_diff.m_changed_modes[m_diff.m_nchanged_mode] = imode;
    ++m_diff.m_nchanged_mode;
}