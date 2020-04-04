//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_CONNECTION_H
#define M7_CONNECTION_H

#include <src/core/fermion/Determinant.h>

/*
 * <bra| ...cre... ...des... |ket>
 * the m_des array contains the spin orbital indices in the ket but not the bra
 * the m_cre array contains the spin orbital indices in the bra but not the ket
 */

struct Connection {
    const size_t m_nbit;
    const size_t m_element_dsize;
    defs::inds m_des, m_cre;
    size_t m_ndes, m_ncre;

    Connection(const DeterminantElement &ket, const DeterminantElement &bra);

    virtual void update(const DeterminantElement &ket, const DeterminantElement &bra);

    const size_t &nexcit() const;

protected:
    void update_connection(const DeterminantElement &ket, const DeterminantElement &bra);
};

/*
 * a connection in which the common indices and antisymmetric phase is computed
 */
struct AntisymConnection : Connection {
    defs::inds m_com;
    size_t m_ncom;
    bool m_phase;

    AntisymConnection(const DeterminantElement &ket, const DeterminantElement &bra);

    void update(const DeterminantElement &ket, const DeterminantElement &bra) override;

protected:
    void update_antisym(const DeterminantElement &ket, const DeterminantElement &bra);
};

#endif //M7_CONNECTION_H
