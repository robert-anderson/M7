//
// Created by rja on 11/08/2021.
//

#ifndef M7_RDM_H
#define M7_RDM_H

#include "src/core/mae/MaeTable.h"
#include "src/core/field/Fields.h"
#include "src/core/table/Communicator.h"
#include "src/core/io/Archivable.h"
#include "FermionPromoter.h"

using namespace conn_utils;

class Rdm : Communicator<MaeRow, MaeRow, true> {
    const size_t m_ranksig;
    const size_t m_rank, m_nfrm_cre, m_nfrm_ann, m_nbos_cre, m_nbos_ann;
    std::vector<FermionPromoter> m_frm_promoters;
    buffered::MaeInds m_lookup_inds;

    static size_t nrow_estimate(size_t nfrm_cre, size_t nfrm_ann, size_t nbos_cre, size_t nbos_ann, size_t nsite);

    static size_t nrow_estimate(size_t exsig, size_t nsite);

public:
    Rdm(const fciqmc_config::Bilinear &opts, size_t ranksig, size_t nsite, size_t nelec, size_t nvalue);

    void make_contribs(const field::FrmOnv &src_onv, const conn::FrmOnv &conn,
                       const FrmOps &com, const defs::wf_t &contrib);

    void make_contribs(const field::FrmBosOnv &src_onv, const conn::FrmBosOnv &conn,
                       const FrmOps &com, const defs::wf_t &contrib);

    void end_cycle();

    void save(hdf5::GroupWriter& gw) const;
};

class Rdms : Archivable {
    std::array<std::unique_ptr<Rdm>, defs::nexsig> m_rdms;
    const defs::inds m_active_ranksigs;
    const std::array<defs::inds, defs::nexsig> m_exsig_ranks;

    std::array<defs::inds, defs::nexsig> make_exsig_ranks() const;

public:
    Rdms(const fciqmc_config::Bilinear &opts, defs::inds ranksigs, size_t nsite, size_t nelec);

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
