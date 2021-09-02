//
// Created by rja on 25/08/2021.
//

#include "ExcitIterGroup.h"

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
