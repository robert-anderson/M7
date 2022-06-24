//
// Created by Robert J. Anderson on 10/08/2021.
//

#include "Rdm.h"

#if 0
const uint_t &FermionRdm::nop() const {
    return m_ncre;
}

uint_t FermionRdm::nrow_estimate(uint_t nann, uint_t ncre, uint_t nsite) {
    double nrow = 1.0;
    nrow *= utils::integer::combinatorial(2 * nsite, nann);
    nrow *= utils::integer::combinatorial(2 * nsite, ncre);
    nrow /= utils::integer::factorial(nann + ncre);
    return nrow;
}

FermionRdm::FermionRdm(const conf::FermionRdm &opts, uint_t nrow_crude_est, uint_t nsite, uint_t nelec) :
base_t(
        log::format("{}-body RDM", opts.m_rank),
        nrow_crude_est / mpi::nrank(),
        nrow_crude_est / mpi::nrank(),
        opts.m_buffers, opts.m_load_balancing,
        {
            {opts.m_rank, 1},
            MappedTableBase::nbucket_guess(
                    nrow_crude_est / mpi::nrank(), opts.m_hash_mapping.m_remap_ratio),
                    opts.m_hash_mapping.m_remap_nlookup,
                    opts.m_hash_mapping.m_remap_ratio
                    },
                    {
            {opts.m_rank, 1},
            MappedTableBase::nbucket_guess(
                    nrow_crude_est / mpi::nrank(), opts.m_hash_mapping.m_remap_ratio),
                    opts.m_hash_mapping.m_remap_nlookup,
                    opts.m_hash_mapping.m_remap_ratio
                    }),
                    Archivable("rdm", opts.m_archivable),
                    m_nann(opts.m_rank), m_ncre(opts.m_rank), m_nelec(nelec), m_lookup_inds(opts.m_rank),
                    m_conn(nsite), m_mixed_estimator(opts.m_mixed_estimator), m_com(nsite) {
    m_promoters.reserve(opts.m_rank + 1);
    for (uint_t nins = 0ul; nins <= opts.m_rank; ++nins) m_promoters.emplace_back(nelec + nins - opts.m_rank, nins);
}

void FermionRdm::make_contribs(const fields::FrmOnv &src_onv, const conn::FrmOnv &conn, const FrmOps &com,
                               const defs::wf_t &src_weight, const defs::wf_t &dst_weight) {
    const auto exlvl = conn.m_cre.size();
    if (conn.m_ann.size() > m_nann && conn.m_cre.size() > m_ncre) return;
    const auto nins = nop() - exlvl;
    ASSERT(nins <= nop());

    const auto &promoter = m_promoters[nins];
    for (uint_t icomb = 0ul; icomb < promoter.m_ncomb; ++icomb) {
        auto phase = promoter.apply(icomb, conn, com, m_lookup_inds);
        if (!exlvl || !nins) {ASSERT(!phase); }
        else {
            if (m_lookup_inds.m_ann[0] == m_lookup_inds.m_cre[0]) ASSERT(!phase);
            if (m_lookup_inds.m_ann[0] == m_lookup_inds.m_cre[1]) ASSERT(phase);
            if (m_lookup_inds.m_ann[1] == m_lookup_inds.m_cre[0]) ASSERT(phase);
            if (m_lookup_inds.m_ann[1] == m_lookup_inds.m_cre[1]) ASSERT(!phase);
        }
        auto irank_send = m_ra.get_rank(m_lookup_inds);
        ASSERT(m_lookup_inds.is_ordered());
        auto &send_table = send(irank_send);
        uint_t irow = *send_table[m_lookup_inds];
        if (irow == ~0ul) irow = send_table.insert(m_lookup_inds);

        send_table.m_row.jump(irow);
        /*
         * include the Fermi phase of the excitation
         */
        phase = phase == conn.phase(src_onv);
        auto contrib = (phase ? -1.0 : 1.0) * src_weight * dst_weight;
        send_table.m_row.m_values[0] += contrib;
    }
}

void FermionRdm::load_fn(hdf5::GroupReader &parent) {

}

void FermionRdm::save_fn(hdf5::GroupWriter &parent) {
    m_store.save(parent, std::to_string(nop()), h5_field_names());
}
#endif