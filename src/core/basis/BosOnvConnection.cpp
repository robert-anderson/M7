//
// Created by RJA on 11/09/2020.
//

#include "BosOnvConnection.h"

BosOnvConnection::Diff::Diff(size_t nmode) : m_changed_modes(nmode, 0ul), m_changes(nmode, 0) {}

void BosOnvConnection::Diff::zero() {
    m_nchanged_mode = 0ul;
}

BosOnvConnection::BosOnvConnection(const size_t& nmode) :
        m_nmode(nmode), m_com(nmode), m_diff(nmode) {}

BosOnvConnection::BosOnvConnection(const fields::BosOnv &in, const fields::BosOnv &out) :
        BosOnvConnection(in.nelement()) {
    connect(in, out);
}

BosOnvConnection::BosOnvConnection(const fields::BosOnv &in) :
        BosOnvConnection(in.nelement()) {
    connect(in, in);
}

const size_t &BosOnvConnection::nchanged_mode() const {
    return m_diff.m_nchanged_mode;
}

const size_t &BosOnvConnection::changed_mode(const size_t &ichange) const {
    ASSERT(ichange < nchanged_mode());
    return m_diff.m_changed_modes[ichange];
}

const int &BosOnvConnection::changes(const size_t &ichange) const {
    ASSERT(ichange < nchanged_mode());
    return m_diff.m_changes[ichange];
}

const int &BosOnvConnection::com(const size_t &icom) const {
    ASSERT(icom < m_nmode);
    return m_com[icom];
}

void BosOnvConnection::connect(const fields::BosOnv &in, const fields::BosOnv &out) {
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

void BosOnvConnection::apply(const fields::BosOnv &in) {
    for (size_t imode = 0ul; imode < m_nmode; ++imode) {
        m_com[imode] = in[imode];
    }
    for (size_t ichange = 0ul; ichange < nchanged_mode(); ++ichange) {
        const auto imode = changed_mode(ichange);
        const auto change = changes(ichange);
        if (change<0) m_com[imode] += change;
    }
}

void BosOnvConnection::apply(const fields::BosOnv &in, fields::BosOnv &out) {
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

void BosOnvConnection::add(const size_t& imode, const int& change) {
    m_diff.m_changes[m_diff.m_nchanged_mode] = change;
    m_diff.m_changed_modes[m_diff.m_nchanged_mode] = imode;
    ++m_diff.m_nchanged_mode;
}