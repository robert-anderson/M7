//
// Created by rja on 11/08/2021.
//

#ifndef M7_RDM_H
#define M7_RDM_H

#include "src/core/mae/MaeTable.h"
#include "src/core/field/Fields.h"
#include "src/core/table/Communicator.h"
#include "src/core/io/Archivable.h"

using namespace conn_utils;

class Rdm : Communicator<MaeRow, MaeRow, true> {

    static size_t nrow_estimate(size_t nfrm_cre, size_t nfrm_ann, size_t nbos_cre, size_t nbos_ann, size_t nsite) {
        double nrow = 1.0;
        nrow *= integer_utils::combinatorial(2 * nsite, nfrm_cre);
        nrow *= integer_utils::combinatorial(2 * nsite, nfrm_ann);
        nrow *= integer_utils::combinatorial(nsite, nbos_cre);
        nrow *= integer_utils::combinatorial(nsite, nbos_ann);
        nrow /= integer_utils::factorial(nfrm_cre + nfrm_ann);
        nrow /= integer_utils::factorial(nbos_cre + nbos_ann);
        return nrow;
    }

    static size_t nrow_estimate(size_t exsig, size_t nsite) {
        return nrow_estimate(decode_nfrm_cre(exsig), decode_nfrm_ann(exsig),
                             decode_nbos_cre(exsig), decode_nbos_ann(exsig), nsite);
    }

    const size_t m_ranksig;
public:
    Rdm(const fciqmc_config::Bilinear &opts, size_t ranksig, size_t nsite, size_t nvalue) :
            Communicator<MaeRow, MaeRow, true>(
                    "rdm_" + to_string(ranksig),
                    nrow_estimate(ranksig, nsite),
                    nrow_estimate(ranksig, nsite),
                    opts.m_buffers, opts.m_load_balancing,
                    {{ranksig, nvalue}}, {{ranksig, nvalue}}
            ),
            m_ranksig(ranksig) {

//        m_promoters.reserve(opts.m_rank + 1);
//        for (size_t nins = 0ul; nins <= opts.m_rank; ++nins) m_promoters.emplace_back(nelec + nins - opts.m_rank, nins);

    }

    void save(hdf5::GroupWriter& gw) const {
        m_store.save(gw, to_string(m_ranksig));
    }
};

class Rdms : Archivable {
    std::array<std::unique_ptr<Rdm>, defs::nexsig> m_rdms;
    const defs::inds m_active_ranksigs;
    const std::array<defs::inds, defs::nexsig> m_exsig_ranks;

    std::array<defs::inds, defs::nexsig> make_exsig_ranks() const {
        std::array<defs::inds, defs::nexsig> exsig_ranks{};
        for (auto &ranksig: m_active_ranksigs) {
            auto nfrm_cre = decode_nfrm_cre(ranksig);
            auto nfrm_ann = decode_nfrm_cre(ranksig);
            while (nfrm_cre != ~0ul && nfrm_ann != ~0ul) {
                auto nbos_cre = decode_nfrm_cre(ranksig);
                auto nbos_ann = decode_nfrm_cre(ranksig);
                while (nbos_cre != ~0ul && nbos_ann != ~0ul) {
                    exsig_ranks[encode_exsig(nfrm_cre, nfrm_ann, nbos_cre, nbos_ann)].push_back(ranksig);
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
    Rdms(const fciqmc_config::Bilinear &opts, defs::inds ranksigs, size_t nsite) :
            Archivable("rdms", opts.m_archivable),
            m_active_ranksigs(std::move(ranksigs)), m_exsig_ranks(make_exsig_ranks()) {
        for (const auto &ranksig: m_active_ranksigs) {
            REQUIRE_TRUE(ranksig, "multidimensional estimators require a nonzero number of SQ operator indices");
            m_rdms[ranksig] = std::unique_ptr<Rdm>(new Rdm(opts, ranksig, nsite, 1ul));
        }
    }

private:
    void load_fn(hdf5::GroupReader &parent) override {

    }

    void save_fn(hdf5::GroupWriter &parent) override {
        hdf5::GroupWriter gw("rdms", parent);
        for (const auto& i: m_active_ranksigs) {
            DEBUG_ASSERT_TRUE(m_rdms[i].get(), "active ranksig was not allocated!");
            m_rdms[i]->save(gw);
        }
    }
};

#endif //M7_RDM_H
