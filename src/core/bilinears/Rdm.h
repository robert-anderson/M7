//
// Created by rja on 11/08/2021.
//

#ifndef M7_RDM_H
#define M7_RDM_H

#include "src/core/field/Fields.h"

struct Rdm {

    const size_t m_exsig;

    /**
     * some RDMs don't get contributions from the square of a CI coefficient, and so do not have "diagonal" elements
     * @return
     *  true if the operator definition conserves particle number
     */
    bool has_diagonals() const {
        return conn_utils::conserve_nfrm(m_exsig)&&conn_utils::conserve_nbos(m_exsig);
    }

    Rdm(size_t exsig): m_exsig(exsig){

    }
};

struct FrmRdm : Rdm {
    FrmRdm(size_t exsig): Rdm(exsig){
        REQUIRE_TRUE(conn_utils::pure_frm(exsig), "Fermion RDM must refer to purely fermionic excitations");
    }
};

struct BosRdm : Rdm {
    BosRdm(size_t exsig): Rdm(exsig){
        REQUIRE_TRUE(conn_utils::pure_bos(exsig), "Boson RDM must refer to purely bosonic excitations");
    }
};

struct FrmBosRdm : Rdm {
    FrmBosRdm(size_t exsig): Rdm(exsig){

    }
};

#endif //M7_RDM_H
