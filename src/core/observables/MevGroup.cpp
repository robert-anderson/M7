//
// Created by rja on 01/04/2021.
//

#include <src/core/util/Foreach.h>
#include "MevGroup.h"

FermionPromoter::FermionPromoter(size_t ncom, size_t nop_insert) :
        m_nop_insert(nop_insert),
        m_ncomb(integer_utils::combinatorial(ncom, nop_insert)),
        m_all_combs(nop_insert * m_ncomb) {
    if (!nop_insert) return;

    foreach::rtnd::Ordered<> foreach_comb(ncom, nop_insert);
    size_t icomb = 0ul;
    auto fn = [&](){
        for (size_t i = 0ul; i < nop_insert; ++i) {
            auto j = icomb * nop_insert + i;
            ASSERT(j < m_all_combs.size());
            m_all_combs[j] = foreach_comb[i];
        }
        ++icomb;
    };
    foreach_comb(fn);
}

const defs::mev_ind_t *FermionPromoter::begin(const size_t &icomb) const {
    ASSERT(icomb < m_ncomb);
    return m_all_combs.data() + icomb * m_nop_insert;
}

bool FermionPromoter::apply(const size_t &icomb, const conn::Antisym<0> &conn, fields::FermionMevInds &inds) const {
    auto comb_begin = begin(icomb);
    inds.zero();
    size_t ann_passed = 0ul;
    size_t cre_passed = 0ul;
    for (size_t iins = 0ul; iins < m_nop_insert; ++iins) {
        auto ins = conn.com(comb_begin[iins]);
        while (ann_passed < conn.nann() && conn.ann(ann_passed) < ins) {
            inds.m_ann[ann_passed + iins] = conn.ann(ann_passed);
            ++ann_passed;
        }
        inds.m_ann[ann_passed + iins] = ins;

        while (cre_passed < conn.ncre() && conn.cre(cre_passed) < ins) {
            inds.m_cre[cre_passed + iins] = conn.cre(cre_passed);
            ++cre_passed;
        }
        inds.m_cre[cre_passed + iins] = ins;
    }
    auto phase = (ann_passed + cre_passed) & 1ul;

    // the rest of the promoted connection is the same as the connection
    while (ann_passed < conn.nann()) {
        inds.m_ann[ann_passed + m_nop_insert] = conn.ann(ann_passed);
        ++ann_passed;
    }
    while (cre_passed < conn.ncre()) {
        inds.m_cre[cre_passed + m_nop_insert] = conn.cre(cre_passed);
        ++cre_passed;
    }
    return phase;
}

const size_t &FermionRdm::nop() const {
    return m_ncre;
}

size_t FermionRdm::nrow_guess(size_t nann, size_t ncre, size_t nsite) {
    double nrow = 1.0;
    nrow *= integer_utils::combinatorial(2 * nsite, nann);
    nrow *= integer_utils::combinatorial(2 * nsite, ncre);
    nrow /= integer_utils::factorial(nann + ncre);
    return nrow;
}

FermionRdm::FermionRdm(const fciqmc_config::FermionRdm &opts, size_t nsite, size_t nelec) :
        base_t(
                log::format("{}-body RDM", opts.m_rank),
                opts.m_buffers, opts.m_load_balancing,
                {
                        {opts.m_rank, 1},
                        MappedTableBase::nbucket_guess(
                                nrow_guess(opts.m_rank, opts.m_rank, nsite) / mpi::nrank(), 3)
                },
                {
                        {opts.m_rank, 1},
                        MappedTableBase::nbucket_guess(
                                nrow_guess(opts.m_rank, opts.m_rank, nsite) / mpi::nrank(), 3)
                }),
        m_nann(opts.m_rank), m_ncre(opts.m_rank), m_nelec(nelec), m_lookup_inds(opts.m_rank),
        m_conn(nsite), m_mixed_estimator(opts.m_mixed_estimator) {
    m_store.resize(100);
    m_comm.resize(100);
    m_promoters.reserve(opts.m_rank+1);
    for (size_t nins=0ul; nins<=opts.m_rank; ++nins) m_promoters.emplace_back(nelec+nins-opts.m_rank, nins);
}

void FermionRdm::make_contribs(const conn::Antisym<0> &conn, const defs::wf_t &src_weight,
                               const defs::wf_t &dst_weight) {
    const auto exlvl = conn.nexcit();
    if (conn.nann() > m_nann && conn.ncre() > m_ncre) return;
    const auto nins = nop()-exlvl;
    ASSERT(nins<=nop());

    const auto& promoter = m_promoters[nins];
    for (size_t icomb=0ul; icomb<promoter.m_ncomb; ++icomb){
        auto phase = promoter.apply(icomb, conn, m_lookup_inds);
        if (!exlvl || !nins) {ASSERT(!phase);}
        else {
            if (m_lookup_inds.m_ann[0]==m_lookup_inds.m_cre[0]) ASSERT(!phase);
            if (m_lookup_inds.m_ann[0]==m_lookup_inds.m_cre[1]) ASSERT(phase);
            if (m_lookup_inds.m_ann[1]==m_lookup_inds.m_cre[0]) ASSERT(phase);
            if (m_lookup_inds.m_ann[1]==m_lookup_inds.m_cre[1]) ASSERT(!phase);
        }
        auto irank_send = m_ra.get_rank(m_lookup_inds);
        ASSERT(m_lookup_inds.is_ordered());
        auto& send_table = send(irank_send);
        size_t irow = *send_table[m_lookup_inds];
        if (irow == ~0ul) irow = send_table.insert(m_lookup_inds);

        send_table.m_row.jump(irow);
        /*
         * include the Fermi phase of the excitation
         */
        phase = phase==conn.phase();
        auto contrib = (phase ? -1.0 : 1.0) * src_weight * dst_weight;
        send_table.m_row.m_values[0] += contrib;
    }
}