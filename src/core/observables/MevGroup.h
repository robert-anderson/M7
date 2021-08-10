//
// Created by rja on 01/04/2021.
//

#ifndef M7_MEVGROUP_H
#define M7_MEVGROUP_H

#include <src/core/io/Options.h>
#include <src/core/enumerator/CombinationEnumerator.h>
#include <src/core/hamiltonian/FermionHamiltonian.h>
#include "src/core/table/Communicator.h"
#include "RefExcits.h"


struct MevGroup {
    const fciqmc_config::AvEsts &m_opts;
    Epoch m_accum_epoch;
    std::unique_ptr<FermionRdm> m_fermion_rdm;
    /**
     * optionally accumulate averaged walker occupations of excitations of the reference
     */
    std::unique_ptr<RefExcits> m_ref_excits;
    const size_t m_period;
    const bool m_explicit_hf_conns;
    size_t m_icycle_period_start = ~0ul;

    MevGroup(const fciqmc_config::AvEsts &opts, size_t nsite, size_t nelec, bool explicit_hf_conns=true) :
            m_opts(opts), m_accum_epoch("MEV accumulation"),
            m_fermion_rdm(opts.m_fermion_rdm.m_rank ? new FermionRdm(opts.m_fermion_rdm, nsite, nelec) : nullptr),
            m_ref_excits(opts.m_ref_excits.m_max_exlvl ? new RefExcits(opts.m_ref_excits, nsite) : nullptr),
            m_period(opts.m_stats_period), m_explicit_hf_conns(explicit_hf_conns) {
    }


    size_t iperiod(size_t icycle) {
        if (!m_accum_epoch || m_icycle_period_start == ~0ul) return ~0ul;
        return (icycle-m_icycle_period_start)/m_period;
    }

    bool is_period_cycle(size_t icycle) {
        if (!m_accum_epoch) return false;
        if (!m_period) return false;
        if (m_icycle_period_start == ~0ul || m_icycle_period_start == icycle) {
            m_icycle_period_start = icycle;
            return false;
        }
        return !((icycle - m_icycle_period_start) % m_period);
    }

    /**
     * If we are using mixed estimator, then set the ket weight to the trial WF expectation
     * @param weight
     *  average or instantaneous ONV weight
     * @return
     *  weight unchanged if not using the mixed estimator, else signed unit
     */
    defs::wf_t get_ket_weight(const defs::wf_t& weight) const {
        if (m_fermion_rdm->m_mixed_estimator)
            return consts::real(weight) > 0.0 ? 1.0 : -1.0;
        return weight;
    };

    operator bool() const {
        return m_fermion_rdm || m_ref_excits;
    }

    bool is_bilinear() const {
        if (m_fermion_rdm) return (*this) && !m_fermion_rdm->m_mixed_estimator;
        return false;
    }

    void save(hdf5::GroupWriter& parent) const {
        //if (m_fermion_rdm) m_fermion_rdm->save(fw);
        if (m_ref_excits) m_ref_excits->save(parent);
    }

    /*
    void save(size_t icycle){
        auto fname = log::format(m_opts.m_periodic_output.m_path, iperiod(icycle));
        hdf5::FileWriter fw(fname);
        auto c0 = mpi::all_sum(m_ref_excits->m_av_ref[0]);
        log::info("average reference weight at period {}: {}", iperiod(icycle), c0);
        //if (m_fermion_rdm) m_fermion_rdm->save(fw);
        if (m_ref_excits) m_ref_excits->save(fw);
    }

    void save() const {
        hdf5::FileWriter fw(m_opts.m_io.m_save_path);
        //if (m_fermion_rdm) m_fermion_rdm->save(fw);
        if (m_ref_excits) m_ref_excits->save(fw);
    }
     */
};

/*
 *
//    hdf5::FileWriter fw(std::to_string(m_mevs.iperiod(m_icycle)) + "." + m_opts.write_hdf5_fname);
//    hdf5::GroupWriter gw("solver", fw);
//    if (m_mevs.m_fermion_rdm) {
//        hdf5::GroupWriter gw2("rdm", gw);
//        m_mevs.m_fermion_rdm->save(gw2);
//    }
 */

#if 0
struct BilinearMevGroup {
    static constexpr size_t c_max_rank = 6;
    const size_t m_nsite;
    const size_t m_max_rank = 6;
    typedef BufferedTable<MevRow<defs::wf_t>, 1> mev_table_t;
    std::array<std::unique_ptr<mev_table_t>, c_max_rank> m_rdms;
    conn::Antisym<> m_conn;
    buffered::FermionMevInds m_lookup_inds;

    BilinearMevGroup(size_t nsite, size_t max_rank):
        m_nsite(nsite), m_max_rank(max_rank), m_conn(nsite), m_lookup(){
        for (size_t rank=1; rank<=max_rank; ++rank)
            m_rdms[rank] = std::unique_ptr<mev_table_t>(
                    new mev_table_t(std::to_string(rank)+"-body RDM", {{{rank, rank}, 1}, 100}));
    }

    operator bool() const {
        return m_max_rank;
    }

    void make_contribs_spf_ket(const fields::Onv<>& src_onv, const defs::wf_t& src_weight,
                               const fields::Onv<>& dst_onv) {
        m_conn.connect(src_onv, src_onv);
        const size_t rank = 1;
        auto& rdm = *m_rdms[1];
        for (const auto &com: m_conn.com()) {

            m_lookup.m_row.m_inds[0] = com;
            m_lookup.m_row.m_inds[1] = com;

            size_t irow = *m_rdms[rank]->operator[](m_lookup.m_row.m_inds);
            if (irow==~0ul) irow = m_rdms[rank]->insert(m_lookup.m_row.m_inds);
            rdm.m_row.jump(irow);
            rdm.m_row.m_values[0]+=std::abs(src_weight);
        }
    }

    void make_contribs(const fields::Onv<>& src_onv, const defs::wf_t& src_weight,
                       const fields::Onv<>& dst_onv, const defs::wf_t& dst_weight){
        m_conn.connect(src_onv, dst_onv);
        const auto exlvl = m_conn.nexcit();
        if (exlvl>m_max_rank) return;

        for (size_t rank=1ul; rank<=m_max_rank; ++rank) {
            auto &rdm = *m_rdms[rank].get();
            if (exlvl == 0) {
                ASSERT(m_conn.ncom()==m_nsite);
                for (const auto &com: m_conn.com()) {
                    //m_lookup.m_row.set(0, {com});
                    //m_lookup.m_row.set(1, {com});
                    m_lookup.m_row.m_inds[0] = com;
                    m_lookup.m_row.m_inds[1] = com;

                    size_t irow = *m_rdms[rank]->operator[](m_lookup.m_row.m_inds);
                    if (irow==~0ul) irow = m_rdms[rank]->insert(m_lookup.m_row.m_inds);
                    rdm.m_row.jump(irow);
                    rdm.m_row.m_values[0]+=src_weight*dst_weight;
                }
            }
        }
    }
};


#endif //M7_MEVGROUP_H
#endif //M7_MEVGROUP_H
