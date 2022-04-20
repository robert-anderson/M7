//
// Created by rja on 19/04/2022.
//

#ifndef M7_HILBERTDATA_H
#define M7_HILBERTDATA_H

#include <M7_lib/parallel/MPIAssert.h>

/**
 * information relating to the many-body attributes of the fermionic sector
 */
struct FrmHilbertData {
    /**
     * electron number sector in which the fermionic ONVs are restricted (fermion number non-conservation is not
     * implemented)
     */
    const size_t m_nelec;
    /**
     * nalpha-nbeta sector in which the fermionic ONVs are restricted (this is only respected if the Hamiltonian
     * commutes with Sz, as signalled by the kramers_attrs members.
     */
    const int m_ms2_restrict;
    /**
     * numbers of alpha and beta electrons implied by the above parameters
     */
    const size_t m_nalpha, m_nbeta;
private:
    size_t nalpha(){
        size_t spin_odd = std::abs(m_ms2_restrict) % 2;
        REQUIRE_EQ(m_nelec % 2, spin_odd, "incompatible values of ms2 and nelec");
        size_t nalpha = m_nelec / 2 + (std::abs(m_ms2_restrict)) / 2 + spin_odd;
        return m_ms2_restrict >= 0 ? nalpha : m_nelec - nalpha;
    }

public:
    FrmHilbertData(size_t nelec, int ms2_restrict=0):
        m_nelec(nelec), m_ms2_restrict(ms2_restrict), m_nalpha(nalpha()), m_nbeta(m_nelec-m_nalpha){
    }

    size_t nci(size_t nsite, bool spin_restrict) const {
        if (spin_restrict)
            return integer_utils::combinatorial(nsite, m_nalpha)*integer_utils::combinatorial(nsite, m_nbeta);
        return integer_utils::combinatorial(2*nsite, m_nelec);
    }
};

/**
 * information relating to the many-body attributes of the bosonic sector
 */
struct BosHilbertData {
    /**
     * number of bosons in the system (only respected if the Hamiltonian is boson-number conserving, else it is only
     * referenced in setting initial or reference BosonONVs
     */
    const size_t m_nboson;
    /**
     *  number of bosons permitted to occupy any given mode
     */
    const size_t m_nboson_max;
    BosHilbertData(size_t nboson, size_t nboson_max=defs::max_bos_occ): m_nboson(nboson), m_nboson_max(nboson_max){}

    size_t nci(size_t nmode, bool number_conserve) const {
        if (number_conserve) return integer_utils::combinatorial_with_repetition(nmode, m_nboson);
        else return std::pow(m_nboson_max + 1, nmode);
    }
};

/**
 * combines information relating to many-body attributes of the fermionic and bosonic sectors
 */
struct HilbertData {
    /**
     * fermion sector attributes
     */
    FrmHilbertData m_frm;
    /**
     * boson sector attributes
     */
    BosHilbertData m_bos;
    //HilbertData(FrmHilbertData&& frm, BosHilbertData&& bos): m_frm(frm), m_bos(bos){}
};


#endif //M7_HILBERTDATA_H
