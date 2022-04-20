//
// Created by rja on 26/07/2021.
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
     * properties of the single-particle basis
     */
    const BosBasisData m_bd;
    /**
     * properties of the many-body basis
     */
    const BosHilbertData m_hd;
    /**
     * total number of bosons in the system (only respected if H is number conserving)
     */
    const size_t m_nboson;
    /**
     * total number of bosons permitted to occupy any given mode
     */
    const size_t m_nboson_max;
    ham_data::TermContribs m_contribs_0011;
    ham_data::TermContribs m_contribs_0022;

    /**
     * @return
     *  true if H only consists of nonzero terms with rank signatures of the form 00xx
     */
    bool number_conserving() const {
        return true;
    }
    /**
     * @param bd
     *  parameters of the bosonic single-particle basis. Just as in FrmHam, nboson and nboson_max are not included in
     *  the definition of this object, since they are properties of the many-body basis
     * @param nboson
     *  number of bosons in the system (only respected if the Hamiltonian is boson-number conserving, else it is only
     *  referenced in setting initial or reference BosonONVs
     * @param nboson_max
     *  number of bosons permitted to occupy any given mode
     */
    BosHam(const BosBasisData& bd, size_t nboson, size_t nboson_max=defs::max_bos_occ);
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

    virtual bool enabled() const {
        return true;
    }

    bool disabled() const {
        return !enabled();
    }
};

struct NullBosHam : BosHam {
    NullBosHam() : BosHam({0ul, 0ul}){}

    bool enabled() const override {
        return false;
    }
};


#endif //M7_BOSHAM_H
