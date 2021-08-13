//
// Created by rja on 11/08/2021.
//

#ifndef M7_RDM_H
#define M7_RDM_H

#include "src/core/field/Fields.h"

using namespace conn_utils;

struct Rdm {

    const size_t m_ranksig;

    Rdm(size_t ranksig): m_ranksig(ranksig){

    }
};

class Rdms {
    std::array<std::unique_ptr<Rdm>, defs::nexsig> m_rdms;
    const defs::inds m_active_ranksigs;
    const std::array<defs::inds, defs::nexsig> m_exsig_ranks;

    std::array<defs::inds, defs::nexsig> make_exsig_ranks() const {
        std::array<defs::inds, defs::nexsig> exsig_ranks{};
        for (auto& ranksig: m_active_ranksigs){
            auto nfrm_cre = extract_ncref(ranksig);
            auto nfrm_ann = extract_ncref(ranksig);
            while(nfrm_cre!=~0ul && nfrm_ann!=~0ul){
                auto nbos_cre = extract_ncref(ranksig);
                auto nbos_ann = extract_ncref(ranksig);
                while(nbos_cre!=~0ul && nbos_ann!=~0ul){
                    exsig_ranks[exsig(nfrm_cre, nfrm_ann, nbos_cre, nbos_ann)].push_back(ranksig);
                    --nbos_cre;
                    --nbos_ann;
                }
                --nfrm_cre;
                --nfrm_ann;
            }
        }
        return exsig_ranks;
    }

public:
    Rdms(defs::inds ranksigs) : m_active_ranksigs(std::move(ranksigs)), m_exsig_ranks(make_exsig_ranks()){
        for (const auto& ranksig: m_active_ranksigs) {
            REQUIRE_TRUE(ranksig, "multidimensional estimators require a nonzero number of SQ operator indices");
            m_rdms[ranksig] = std::unique_ptr<Rdm>(new Rdm(ranksig));
        }
    }
};

#endif //M7_RDM_H
