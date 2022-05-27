//
// Created by Robert J. Anderson on 26/07/2021.
//

#ifndef M7_BOSHAM_H
#define M7_BOSHAM_H

#include "M7_lib/integrals/BosonCoeffs_1.h"
#include "M7_lib/integrals/BosonCoeffs_2.h"
#include "M7_lib/io/BosdumpFileReader.h"
#include "M7_lib/connection/Connections.h"
#include "M7_lib/parallel/SharedArray.h"
#include "M7_lib/field/Fields.h"

#include "M7_lib/hamiltonian/HamiltonianData.h"
#include "M7_lib/hamiltonian/HamOpTerm.h"

struct BosHam : HamOpTerm {
    /**
     * a convenient pair of references to the relevant Hamiltonian section and the Basis configuration section
     */
    typedef HamOpTerm::OptPair<conf::BosHam> opt_pair_t;
    /**
     * properties of the many-body basis
     */
    const sys::bos::Basis m_basis;

    ham_data::TermContribs m_contribs_0011;
    ham_data::TermContribs m_contribs_0022;

    BosHam(const sys::bos::Basis& basis):
            m_basis(basis), m_contribs_0011(exsig_utils::ex_0011), m_contribs_0022(exsig_utils::ex_0022) {}

private:
    /**
     * particle sector information can derive from the hamiltonian definition or from the configuration document, and
     * sometimes these can differ. the configuration document is treated as the authority, and any definitions it makes
     * override any that may have already been inferred from the hamiltonian definition
     * @param ham_elecs
     *  bosons inferred from the H definition
     * @param conf_elecs
     *  bosons provided in the configuration document (authoritative)
     * @return
     *  combined bosons object
     */
//    static sys::bos::Bosons make_bosons(const sys::bos::Bosons& from_ham, const sys::bos::Bosons& from_conf){
//        const size_t nboson = from_conf.has_value() ? from_conf : from_ham;
//        const bool conserve = from_ham.conserve();
//        const size_t occ_cutoff = from_conf.m_occ_cutoff;
//        return {nboson, conserve, occ_cutoff};
//    }

public:

//    BosHam(const sys::bos::Basis& sector, const sys::bos::Bosons& from_conf) :
//            BosHam(sys::bos::Particles(sector.m_basis, make_bosons(sector.m_bosons, from_conf))){}

	virtual ~BosHam(){}

    virtual defs::ham_t get_coeff_0011(size_t i, size_t j) const {return 0;}
    virtual defs::ham_t get_coeff_0022(size_t i, size_t j,
                                       size_t k, size_t l) const {return 0;}

    virtual defs::ham_t get_element_0000(const field::BosOnv &onv) const {return 0;}
    virtual defs::ham_t get_element_0011(const field::BosOnv &onv, const conn::BosOnv& conn) const {return 0;}
    virtual defs::ham_t get_element_0022(const field::BosOnv &onv, const conn::BosOnv& conn) const {return 0;}

    defs::ham_t get_element(const field::BosOnv &onv) const {
        return get_element_0000(onv);
    }

    defs::ham_comp_t get_energy(const field::BosOnv &onv) const {
        return consts::real(get_element(onv));
    }

    defs::ham_t get_element(const field::BosOnv &src, const conn::BosOnv& conn) const {
        switch (conn.size()) {
            case 0: return get_element_0000(src);
            case 2: return get_element_0011(src, conn);
            case 4: return get_element_0022(src, conn);
        }
        return 0;
    }

    virtual void log_data() const;

    virtual sys::bos::Bosons default_bosons() const {
        return {0ul, false};
    }

};

struct NullBosHam : BosHam {
    NullBosHam() : BosHam(0ul){}

    bool enabled() const override {
        return false;
    }
};


#endif //M7_BOSHAM_H
