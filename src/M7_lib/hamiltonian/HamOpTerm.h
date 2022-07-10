//
// Created by Robert J. Anderson on 03/04/2022.
//

#ifndef M7_HAMOPTERM_H
#define M7_HAMOPTERM_H


#include <forward_list>
#include "M7_lib/foreach/ConnForeach.h"
#include "M7_lib/excitgen/ExcitGen.h"
#include "M7_lib/conf/Conf.h"

/**
 * for multiple inheritance indicating that all elements of the HamOpTerm are zero
 */
struct NullOpTerm {};

/**
 * some Hamiltonian terms are only valid for expectation values computed between MBFs in a specific electron number
 * sector. This is indicated by multiple inheritance with SectoredTerm as a superclass
 */
struct ElecSpecTerm {
    sys::frm::Electrons m_elecs;
    ElecSpecTerm(sys::frm::Electrons elecs): m_elecs(std::move(elecs)){}
};

/**
 * base class for the three currently implemented kinds of hamiltonian term based on the second quantised operators
 * which appear in their definitions
 *  - frm/FrmHam: fermion number conserving, boson-independent terms in H
 *  - bos/BosHam: boson number conserving or non-conserving, fermion-independent terms in H
 *  - frmbos/FrmBosHam: fermion number conserving, boson number conserving or non-conserving product terms in H
 */
struct HamOpTerm {

    template<typename ham_opt_t>
    struct OptPair {
        static_assert(std::is_base_of<yaml::Section, ham_opt_t>::value,
                "template arg must be derived from yaml::Section");
        const ham_opt_t& m_ham;
        const conf::Basis& m_basis;
    };

    HamOpTerm(){}
    virtual ~HamOpTerm(){}

    template<typename T>
    const T* as() const {
        return dynamic_cast<const T*>(this);
    }

    template<typename T>
    T* as() {
        return dynamic_cast<T*>(this);
    }

    template<typename T>
    bool is() const {
        return as<T>();
    }

    operator bool() const {
        return !is<NullOpTerm>();
    }

    bool is_elec_spec() const {
        return is<ElecSpecTerm>();
    }

    using excit_gen_ptr_t = ExcitGen::excit_gen_ptr_t;
    using excit_gen_list_t = ExcitGen::excit_gen_list_t;
    using conn_foreach_ptr_t = conn_foreach::base_ptr_t;
    using conn_foreach_list_t = conn_foreach::base_list_t;
    /**
     * @param prng
     *  pseudorandom number generator to pass to the ctors of all associated excit gens
     * @param opts
     *  some derived classes of HamOpTerm may be compatible with multiple excitation generation methods, these can be
     *  are configured at input
     * @return
     *  forward linked list of excitation generators applicable to this HamOpTerm
     */
    virtual excit_gen_list_t make_excit_gens(PRNG& /*prng*/, const conf::Propagator& /*opts*/) const {
        return {};
    }
    /**
     * @return
     *  forward linked list of foreach iterators over the connections of a given MBF
     */
    virtual conn_foreach::base_list_t make_foreach_iters() const {
        return {};
    }
};


#endif //M7_HAMOPTERM_H
