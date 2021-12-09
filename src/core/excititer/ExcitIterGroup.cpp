//
// Created by rja on 25/08/2021.
//

#include <src/core/hamiltonian/HubbardHamiltonian.h>
#include "ExcitIterGroup.h"
#include "Hubbard1dSingles.h"
#include "LadderPure.h"
#include "LadderPureHolstein.h"
#include "LadderHopping.h"

void ExcitIterGroup::init() {
    for (size_t exsig=0ul; exsig<defs::nexsig; ++exsig){
        auto ptr = m_excit_iters[exsig].get();
        if (!ptr) continue;
        if (exsig_utils::is_pure_frm(exsig)) m_frm_inds.push_back(m_active_exsigs.size());
        m_active_exsigs.push_back(exsig);
    }
    log_breakdown();
}

void ExcitIterGroup::add(std::unique_ptr<ExcitIter> &&excit_iter) {
    auto exsig = excit_iter->m_exsig;
    REQUIRE_TRUE(m_excit_iters[exsig] == nullptr,
                 "can't specify more than one excitation iterator for the same exsig in an ExcitIterGroup");
    m_excit_iters[exsig] = std::move(excit_iter);
    REQUIRE_LT(exsig, defs::nexsig, "exsig OOB");
}

void ExcitIterGroup::log_breakdown() const {
    std::vector<std::string> strs;
    for (auto exsig: m_active_exsigs) strs.push_back(exsig_utils::to_string(exsig));
    log::info("Excitation iterator generating exsigs: {}", utils::to_string(strs));
}

ExcitIterGroup::ExcitIterGroup(const Hamiltonian &ham) {
    if (ham.m_frm) {
        bool any_singles =
                ham.m_frm->m_contribs_1100.is_nonzero(ex_single) || ham.m_frm->m_contribs_2200.is_nonzero(ex_single);
        if (any_singles) {
            bool is_hubbard = dynamic_cast<const HubbardHamiltonian *>(ham.m_frm.get());
            if (is_hubbard) {
                add(std::unique_ptr<ExcitIter>(new excititers::Hubbard1dSingles(ham)));
            } else {
                add(std::unique_ptr<ExcitIter>(new excititers::Frm(ham, ex_single)));
            }
        }
        if (ham.m_frm->m_contribs_2200.is_nonzero(ex_double)) {
            add(std::unique_ptr<ExcitIter>(new excititers::Frm(ham, ex_double)));
        }
    }

    if (ham.m_ladder) {
        bool any_pures;
        any_pures = ham.m_ladder->m_contribs_0010.is_nonzero(exsig_utils::ex_0010)
                    || ham.m_ladder->m_contribs_1110.is_nonzero(exsig_utils::ex_0010);
        if (any_pures) {
            if (ham.m_ladder->is_holstein()) {
                add(std::unique_ptr<ExcitIter>(new excititers::LadderPureHolstein(ham, exsig_utils::ex_0010)));
            }
            else {
                add(std::unique_ptr<ExcitIter>(new excititers::LadderPure(ham, exsig_utils::ex_0010)));
            }
        }
        any_pures = ham.m_ladder->m_contribs_0001.is_nonzero(exsig_utils::ex_0001)
                    || ham.m_ladder->m_contribs_1101.is_nonzero(exsig_utils::ex_0001);
        if (any_pures) {
            if (ham.m_ladder->is_holstein())
                add(std::unique_ptr<ExcitIter>(new excititers::LadderPureHolstein(ham, exsig_utils::ex_0001)));
            else
                add(std::unique_ptr<ExcitIter>(new excititers::LadderPure(ham, exsig_utils::ex_0001)));
        }

        bool any_hopping;
        any_hopping = ham.m_ladder->m_contribs_1110.is_nonzero(exsig_utils::ex_1110);
        if (any_hopping)
            add(std::unique_ptr<ExcitIter>(new excititers::LadderHopping(ham, exsig_utils::ex_1110)));
        any_hopping = ham.m_ladder->m_contribs_1101.is_nonzero(exsig_utils::ex_1101);
        if (any_hopping)
            add(std::unique_ptr<ExcitIter>(new excititers::LadderHopping(ham, exsig_utils::ex_1101)));
    }

    if (ham.m_bos){
        if(ham.m_bos->m_contribs_0022.is_nonzero(exsig_utils::ex_0022)){
            add(std::unique_ptr<ExcitIter>(new excititers::Bos(ham, exsig_utils::ex_0022)));
        }
    }

    init();
}
