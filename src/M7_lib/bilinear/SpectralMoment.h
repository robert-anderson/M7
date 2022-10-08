//
// Created by Robert J. Anderson on 11/08/2021.
//

#ifndef M7_SPECTRALMOMENT_H
#define M7_SPECTRALMOMENT_H

#include <M7_lib/mae/MaeTable.h>
#include <M7_lib/io/Archivable.h>
#include <M7_lib/field/Fields.h>
#include <M7_lib/util/Exsig.h>


//Communicator<MaeRow<wf_t>, MaeRow<wf_t>, true>
struct SpectralMoment {
    const uint_t m_exsig, m_order;
    SpectralMoment(uint_t exsig, uint_t order): m_exsig(exsig), m_order(order){
        REQUIRE_EQ(order, 1ul, "Spectral moment quantities are currently only implemented for n=1")
        REQUIRE_TRUE(exsig::is_pure_frm(exsig),
                     "Spectral moment excitations must refer to purely fermionic perturbations");
    }
};

class SpecMoms : public Archivable {

public:

    SpecMoms(const conf::SpecMoms& opts): Archivable("spec_moms", opts.m_archivable){}

    operator bool() const {
        return false;//!m_active_ranksigs.empty();
    }

    bool all_stores_empty() const {
//        for (auto& ranksig: m_active_ranksigs)
//            if (!m_rdms[ranksig]->m_store.is_freed())
//                return false;
            return true;
    }

protected:
    void load_fn(const hdf5::NodeReader& /*parent*/) override {

    }

    void save_fn(const hdf5::NodeWriter& /*parent*/) override {

    }

};

#endif //M7_SPECTRALMOMENT_H
