//
// Created by rja on 03/04/2022.
//

#ifndef M7_HAMOPTERM_H
#define M7_HAMOPTERM_H


#include <forward_list>
#include "M7_lib/foreach/ConnForeach.h"
#include "M7_lib/excitgen2/ExcitGen2.h"

/**
 * base class for the three currently implemented kinds of hamiltonian term based on the second quantised operators
 * which appear in their definitions
 *  - frm/FrmHam: fermion number conserving, boson-independent terms in H
 *  - bos/BosHam: boson number conserving or non-conserving, fermion-independent terms in H
 *  - frmbos/FrmBosHam: fermion number conserving, boson number conserving or non-conserving product terms in H
 */
struct HamOpTerm {

    HamOpTerm(){}
    virtual ~HamOpTerm(){}

    typedef std::unique_ptr<ExcitGen2> excit_gen_ptr_t;
    typedef std::forward_list<excit_gen_ptr_t> excit_gen_list_t;
    virtual excit_gen_list_t make_excit_gens(PRNG& prng){
        return {};
    }

    typedef std::unique_ptr<conn_foreach::Base> conn_iter_ptr_t;
    typedef std::forward_list<conn_iter_ptr_t> conn_iter_ptr_list_t;
    virtual conn_iter_ptr_list_t make_conn_iters(){
        return {};
    }

    template<typename ham_t>
    static const ham_t* cast(const HamOpTerm* h) {
        static_assert(std::is_base_of<HamOpTerm, ham_t>::value, "template arg must be derived from HamOpTerm");
        return dynamic_cast<const ham_t*>(h);
    }
};


#endif //M7_HAMOPTERM_H
