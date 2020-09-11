//
// Created by RJA on 11/09/2020.
//

#ifndef M7_BOSONCONNECTION_H
#define M7_BOSONCONNECTION_H


#include "PermanentField.h"

class BosonConnection {
    const size_t m_nmode;
    defs::inds m_changed_modes;
    size_t m_nchanged_mode = 0ul;
    std::vector<int> m_occ_diff;

    explicit BosonConnection(const PermanentField* field):
        m_nmode(field->nelement()), m_changed_modes(m_nmode, 0ul), m_occ_diff(m_nmode, 0) {
    }

    BosonConnection(const PermanentElement &ket, const PermanentElement &bra);

    explicit BosonConnection(const PermanentElement &ket);

    void connect(const PermanentElement &ket, const PermanentElement &bra) {
        m_nchanged_mode = 0ul;
        for (size_t imode=0ul; imode<m_nmode; ++imode) {
            //m_occ_diff[imode] = ((int) PermanentElement
        }
    }
};


#endif //M7_BOSONCONNECTION_H
