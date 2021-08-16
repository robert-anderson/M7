//
// Created by rja on 11/08/2021.
//

#ifndef M7_RDM_H
#define M7_RDM_H

#include <src/core/basis/Suites.h>
#include "src/core/mae/MaeTable.h"
#include "src/core/field/Fields.h"
#include "src/core/table/Communicator.h"
#include "src/core/io/Archivable.h"
#include "FermionPromoter.h"

using namespace conn_utils;

class Rdm : public Communicator<MaeRow, MaeRow, true> {
    const size_t m_ranksig;
    const size_t m_rank, m_nfrm_cre, m_nfrm_ann, m_nbos_cre, m_nbos_ann;
    std::vector<FermionPromoter> m_frm_promoters;
    buffered::MaeInds m_lookup_inds;

    static size_t nrow_estimate(size_t nfrm_cre, size_t nfrm_ann, size_t nbos_cre, size_t nbos_ann, size_t nsite);

    static size_t nrow_estimate(size_t exsig, size_t nsite);

public:
    Rdm(const fciqmc_config::Bilinears &opts, size_t ranksig, size_t nsite, size_t nelec, size_t nvalue);

    void make_contribs(const field::FrmOnv &src_onv, const conn::FrmOnv &conn,
                       const FrmOps &com, const defs::wf_t &contrib);

    void make_contribs(const field::FrmBosOnv &src_onv, const conn::FrmBosOnv &conn,
                       const FrmOps &com, const defs::wf_t &contrib);

    void end_cycle();

    void save(hdf5::GroupWriter& gw) const;
};

class Rdms : public Archivable {
    std::array<std::unique_ptr<Rdm>, defs::nexsig> m_rdms;
    const defs::inds m_active_ranksigs;
    const std::array<defs::inds, defs::nexsig> m_exsig_ranks;

    suite::Conns m_work_conns;
    FrmOps m_work_com_ops;

    std::array<defs::inds, defs::nexsig> make_exsig_ranks() const;

public:
    const Epoch& m_accum_epoch;
    Rdms(const fciqmc_config::Bilinears &opts, defs::inds ranksigs, size_t nsite, size_t nelec, const Epoch& accum_epoch);

    operator bool() const {
        return !m_active_ranksigs.empty();
    }

    bool takes_contribs_from(const size_t& exsig) const {
        return !m_exsig_ranks[exsig].empty();
    }

    void make_contribs(const field::FrmOnv &src_onv, const conn::FrmOnv &conn,
                       const FrmOps &com, const defs::wf_t &contrib){
        auto exsig = conn.exsig();
        for (auto ranksig: m_exsig_ranks[exsig]) m_rdms[ranksig]->make_contribs(src_onv, conn, com, contrib);
    }

    void make_contribs(const field::FrmBosOnv &src_onv, const conn::FrmBosOnv &conn,
                       const FrmOps &com, const defs::wf_t &contrib){
        auto exsig = conn.exsig();
        for (auto ranksig: m_exsig_ranks[exsig]) m_rdms[ranksig]->make_contribs(src_onv, conn, com, contrib);
    }

    void make_contribs(const field::FrmOnv &src_onv, const field::FrmOnv &dst_onv, const defs::wf_t &contrib){
        m_work_conns[src_onv].connect(src_onv, dst_onv, m_work_com_ops);
        make_contribs(src_onv, m_work_conns[src_onv], m_work_com_ops, contrib);
    }

    void make_contribs(const field::FrmBosOnv &src_onv, const field::FrmBosOnv &dst_onv, const defs::wf_t &contrib){
        m_work_conns[src_onv].connect(src_onv, dst_onv, m_work_com_ops);
        make_contribs(src_onv, m_work_conns[src_onv], m_work_com_ops, contrib);
    }

    bool all_stores_empty() const {
        for (auto& ranksig: m_active_ranksigs)
            if (!m_rdms[ranksig]->m_store.is_cleared())
                return false;
        return true;
    }

    void end_cycle() {
        for (auto& ranksig: m_active_ranksigs) m_rdms[ranksig]->end_cycle();
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
