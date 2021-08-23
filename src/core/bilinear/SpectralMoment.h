//
// Created by rja on 11/08/2021.
//

#ifndef M7_SPECTRALMOMENT_H
#define M7_SPECTRALMOMENT_H

#include <src/core/mae/MaeTable.h>
#include <src/core/io/Archivable.h>
#include "src/core/field/Fields.h"


//Communicator<MaeRow<defs::wf_t>, MaeRow<defs::wf_t>, true>
struct SpectralMoment {
    const size_t m_exsig, m_order;
    SpectralMoment(size_t exsig, size_t order): m_exsig(exsig), m_order(order){
        REQUIRE_EQ(order, 1ul, "Spectral moment quantities are currently only implemented for n=1")
        REQUIRE_TRUE(exsig_utils::is_pure_frm(exsig),
                     "Spectral moment excitations must refer to purely fermionic perturbations");
    }
};

class SpecMoms : public Archivable {

public:

    SpecMoms(const fciqmc_config::SpecMoms& opts): Archivable("spec_moms", opts.m_archivable){}

    operator bool() const {
        return false;//!m_active_ranksigs.empty();
    }

    bool all_stores_empty() const {
//        for (auto& ranksig: m_active_ranksigs)
//            if (!m_rdms[ranksig]->m_store.is_cleared())
//                return false;
            return true;
    }

protected:
    void load_fn(hdf5::GroupReader &parent) override {

    }

    void save_fn(hdf5::GroupWriter &parent) override {

    }

};

#endif //M7_SPECTRALMOMENT_H
