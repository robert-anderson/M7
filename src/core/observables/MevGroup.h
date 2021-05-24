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
    const defs::mev_ind_t *begin(const size_t &icomb) const;

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
    bool apply(const size_t &icomb, const conn::Antisym<0> &conn, fields::FermionMevInds &inds) const;
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

    conn::Antisym<0> m_conn;

    const size_t &nop() const;

    static size_t nrow_guess(size_t nann, size_t ncre, size_t nsite);

    FermionRdm(const Options &opts, size_t nop, size_t nsite, size_t nelec);

    void make_contribs(const conn::Antisym<0> &conn, const defs::wf_t &src_weight, const defs::wf_t &dst_weight);

    void make_contribs(const fields::FermionOnv &src_onv, const defs::wf_t &src_weight,
                       const fields::FermionOnv &dst_onv, const defs::wf_t &dst_weight) {
        m_conn.connect(src_onv, dst_onv);
        make_contribs(m_conn, src_weight, dst_weight);
    }

    void make_contribs(const fields::Onv<1> &src_onv, const defs::wf_t &src_weight,
                       const fields::Onv<1> &dst_onv, const defs::wf_t &dst_weight) {
        make_contribs(src_onv.m_frm, src_weight, dst_onv.m_frm, dst_weight);
    }

    void make_contribs(const fields::FermionOnv &src_onv, const defs::wf_t &src_weight,
                       const fields::FermionOnv &dst_onv, const defs::wf_t &dst_weight, const size_t &nop_conn) {
        m_conn.connect(src_onv, dst_onv);
        if (m_conn.nexcit() != nop_conn) return;
        make_contribs(m_conn, src_weight, dst_weight);
    }

    void make_contribs_spf_ket(const conn::Antisym<0> &conn, const defs::wf_t &src_weight);

    void make_contribs_spf_ket(const fields::FermionOnv &src_onv, const defs::wf_t &src_weight,
                               const fields::FermionOnv &dst_onv) {
        m_conn.connect(src_onv, dst_onv);
        make_contribs_spf_ket(m_conn, src_weight);
    }

    void end_cycle() {
        if (!send().buffer_dsize()) return;
        communicate();
        auto &row = m_comm.recv().m_row;
        if (!m_comm.recv().m_hwm) return;
        for (row.restart(); row.in_range(); row.step()) {
            auto irow_store = *m_store[row.m_inds];
            if (irow_store == ~0ul) irow_store = m_store.insert(row.m_inds);
            m_store.m_row.jump(irow_store);
            m_store.m_row.m_values += row.m_values;
        }
        m_comm.recv().clear();
    }

    void h5_read(hdf5::GroupReader &parent) {
        m_store.clear();
        BufferedTable<MevRow<defs::wf_t>> m_buffer("", {{m_nann, m_ncre}});
        m_buffer.push_back();
        RowHdf5Reader<MevRow<defs::wf_t>> row_reader(m_buffer.m_row, parent, std::to_string(nop()), h5_field_names());

        row_reader.restart();
        for (size_t iitem = 0ul; iitem < row_reader.m_nitem; ++iitem) {
            row_reader.read(iitem);
            auto &send_table = send(m_ra.get_rank(row_reader.m_inds));
            // should never read in the same inds twice
            ASSERT(!send_table[row_reader.m_inds]);
            auto irow = send_table.insert(row_reader.m_inds);
            send_table.m_row.jump(irow);
            send_table.m_row.m_values = row_reader.m_values;
        }
    }

    void h5_write(hdf5::GroupWriter &parent) {
        m_store.write(parent, std::to_string(nop()), h5_field_names());
    }

    std::vector<std::string> h5_field_names() {
        return {m_store.m_row.m_inds.m_ann.m_name,
                m_store.m_row.m_inds.m_cre.m_name,
                m_store.m_row.m_values.m_name};
    }
};

struct MevGroup {
    Epoch m_accum_epoch;
    std::unique_ptr<FermionRdm> m_fermion_rdm;
    const size_t m_period;
    const bool m_explicit_hf_conns, m_output_periodically;
    size_t m_icycle_period_start = ~0ul;

    MevGroup(const Options &opts, size_t nsite, size_t nelec) :
            m_accum_epoch("MEV accumulation"),
            m_fermion_rdm(opts.rdm_rank ? new FermionRdm(opts, opts.rdm_rank, nsite, nelec) : nullptr),
            m_period(opts.ncycle_mev_period), m_explicit_hf_conns(opts.explicit_hf_conn_mevs),
            m_output_periodically(opts.output_mevs_periodically) {}


    size_t iperiod(size_t icycle) {
        if (!m_accum_epoch || m_icycle_period_start == ~0ul) return ~0ul;
        return (icycle-m_icycle_period_start)/m_period;
    }

    bool is_period_cycle(size_t icycle) {
        if (!m_accum_epoch) return false;
        if (m_icycle_period_start == ~0ul || m_icycle_period_start == icycle) {
            m_icycle_period_start = icycle;
            return false;
        }
        return !((icycle - m_icycle_period_start) % m_period);
    }
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
