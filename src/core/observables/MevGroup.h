//
// Created by rja on 01/04/2021.
//

#ifndef M7_MEVGROUP_H
#define M7_MEVGROUP_H

#include <src/core/io/Options.h>
#include <src/core/enumerator/CombinationEnumerator.h>
#include "src/core/table/Communicator.h"
#include "AverageCoefficients.h"


/**
 * given a Connection, iterate over common operators
 *
 * Fermionic and bosonic promotion require sorting of indices,
 * fermionic additionally requires the number of swaps involved in the sort to be
 * recorded in order to compute the phase of the promotion
 *
 * e.g.
 *     01234 56789
 *    (01101,01011)
 *    (00111,11001)
 *    ( x   ,   x )
 *    ann: [1, 8]
 *    (   x ,x    )
 *    cre: [3, 5]
 *    (  x x, x  x)
 *    com: [2, 4, 6, 9]
 *    phase: -1
 *
 *    promotion to 3-body RDM element:
 *    enumerate common indices i and count number of set indices between i and the beginning of
 *    the string.
 */
struct FermionPromoter {
    /**
     * number of common creation-annihilation operator pairs to insert into connection
     */
    const size_t m_nop_insert;
    /**
     * total number of possible combinations
     */
    const size_t m_ncomb;
    /**
     * enumeration of all possible combinations
     */
    std::vector<defs::mev_ind_t> m_all_combs;
    FermionPromoter(size_t ncom, size_t nop_insert);

private:
    /**
     * @param icomb
     *  combination index
     * @return
     *  const pointer to beginning of combination in m_all_combs
     */
    const defs::mev_ind_t* begin(const size_t& icomb) const;

public:
    /**
     * apply the icomb-th promotion to the connection given and store the result inds
     * @param icomb
     *  combination index
     * @param conn
     *  connection (no repeated SQ operator indices between ann and cre vectors)
     * @param inds
     *  MEV index field
     * @return
     *  antisymmetric phase associated with sorting both ann and cre to ascending order
     */
    bool apply(const size_t& icomb, const conn::Antisym<>& conn, fields::FermionMevInds& inds) const;
};

struct FermionRdm : Communicator<MevRow<defs::wf_t>, MevRow<defs::wf_t>, true> {
    typedef Communicator<MevRow<defs::wf_t>, MevRow<defs::wf_t>, true> base_t;
    const size_t m_nann, m_ncre, m_nelec;
    /**
     * working indices for building promotions and looking up the MEV tables
     */
    buffered::FermionMevInds m_lookup_inds;
    /**
     * all promoters required for an RDM of this rank. Index refers to the number of SQ ops
     * being inserted
     */
     std::vector<FermionPromoter> m_promoters;

     const size_t& nop() const {
         return m_ncre;
     }

    static size_t nrow_guess(size_t nann, size_t ncre, size_t nsite) {
        double nrow = 1.0;
        nrow *= integer_utils::combinatorial(2 * nsite, nann);
        nrow *= integer_utils::combinatorial(2 * nsite, ncre);
        nrow /= integer_utils::factorial(nann + ncre);
        return nrow;
    }

    FermionRdm(const Options &opts, size_t nop, size_t nsite, size_t nelec) :
            base_t(
                    log::format("{}-body RDM", nop),
                    opts.mev_buffer_expansion_factor,
                    opts.nload_balance_block_per_rank * mpi::nrank(),
                    opts.load_balance_period,
                    {
                            {nop, 1},
                            MappedTableBase::nbucket_guess(
                                    nrow_guess(nop, nop, nsite) / mpi::nrank(), 3)
                    },
                    {
                            {nop, 1},
                            MappedTableBase::nbucket_guess(
                                    nrow_guess(nop, nop, nsite) / mpi::nrank(), 3)
                    },
                    opts.acceptable_load_imbalance),
            m_nann(nop), m_ncre(nop), m_nelec(nelec), m_lookup_inds(nop) {
        m_promoters.reserve(nop+1);
        for (size_t nins=0ul; nins<=nop; ++nins) m_promoters.emplace_back(nelec+nins-nop, nins);
    }

#if 0
    void make_contribs_spf_ket(const fields::Onv<> &src_onv, const defs::wf_t &src_weight,
                               const fields::Onv<> &dst_onv) {
        m_conn.connect(src_onv, src_onv);
        const size_t rank = 1;
        auto &rdm = *m_rdms[1];
        for (const auto &com: m_conn.com()) {

            m_lookup.m_row.m_inds[0] = com;
            m_lookup.m_row.m_inds[1] = com;

            size_t irow = *m_rdms[rank]->operator[](m_lookup.m_row.m_inds);
            if (irow == ~0ul) irow = m_rdms[rank]->insert(m_lookup.m_row.m_inds);
            rdm.m_row.jump(irow);
            rdm.m_row.m_values[0] += std::abs(src_weight);
        }
    }
#endif

    void make_contribs(const conn::Antisym<> &conn, const defs::wf_t &src_weight, const defs::wf_t &dst_weight) {
        const auto exlvl = conn.nexcit();
        if (conn.nann() > m_nann && conn.ncre() > m_ncre) return;
        const auto nins = exlvl-nop();
        ASSERT(nins<=exlvl);

        const auto& promoter = m_promoters[nins];
        for (size_t icomb=0ul; icomb<promoter.m_ncomb; ++icomb){
            auto phase = promoter.apply(icomb, conn, m_lookup_inds);
            auto& send = m_comm.send(m_ra.get_rank(m_lookup_inds));
            size_t irow = send[m_lookup_inds];
            if (irow == ~0ul) irow = send.insert(m_lookup_inds);
            send.m_row.jump(irow);
            auto contrib = src_weight * dst_weight;
            if (!phase) contrib = -contrib;
            send.m_row.m_values[0] += contrib;
        }
    }

    void make_contribs_spf_ket(const conn::Antisym<> &conn, const defs::wf_t &src_weight){
        make_contribs(conn, std::abs(src_weight), 1);
     }
};

struct BilinearMevGroup {

};

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
