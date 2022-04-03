//
// Created by rja on 03/04/2022.
//

#ifndef M7_HAMOPTERM_H
#define M7_HAMOPTERM_H


#include <forward_list>
#include "M7_lib/foreach/ConnForeach.h"

/**
 * base class for the three currently implemented kinds of hamiltonian term based on the second quantised operators
 * which appear in their definitions
 *  - frm/FrmHam: fermion number conserving, boson-independent terms in H
 *  - bos/BosHam: boson number conserving or non-conserving, fermion-independent terms in H
 *  - frmbos/FrmBosHam: fermion number conserving, boson number conserving or non-conserving product terms in H
 */
struct HamOpTerm {

    virtual ~HamOpTerm(){}

    virtual void add_excitgens(){

    }

    typedef std::forward_list<std::unique_ptr<conn_foreach::Base>> conn_foreach_ptr_list_t;
    virtual conn_foreach_ptr_list_t conn_foreach(){
        return {};
    }
};


#endif //M7_HAMOPTERM_H
