//
// Created by rja on 01/04/2021.
//

#ifndef M7_MEVGROUP_H
#define M7_MEVGROUP_H

#include "MevTable.h"
#include "AverageCoefficients.h"

struct BilinearMevGroup {
    static constexpr size_t c_max_rank = 6;
    const size_t m_nsite;
    const size_t m_max_rank = 6;
    typedef BufferedTable<MevRow<defs::wf_t>, 1> mev_table_t;
    typedef BufferedTable<MevRow<defs::wf_t>, 0> mev_lookup_t;
    std::array<std::unique_ptr<mev_table_t>, c_max_rank> m_rdms;
    conn::Antisym<> m_conn;
    mev_lookup_t m_lookup;


    BilinearMevGroup(size_t nsite, size_t max_rank):
        m_nsite(nsite), m_max_rank(max_rank), m_conn(nsite), m_lookup("", {{{1, 1}, 1}}){
        m_lookup.push_back(1);
        m_lookup.m_row.restart();
        for (size_t rank=1; rank<=max_rank; ++rank)
            m_rdms[rank] = std::unique_ptr<mev_table_t>(
                    new mev_table_t(std::to_string(rank)+"-body RDM", {{{rank, rank}, 1}, 100}));
    }

    operator bool() const {
        return m_max_rank;
    }

    void make_contribs(const fields::Onv<>& src_onv, const defs::wf_t& src_weight,
                       const fields::Onv<>& dst_onv, const defs::wf_t& dst_weight){

        m_conn.connect(src_onv, dst_onv);
        const auto exlvl = m_conn.nexcit();
        if (exlvl>m_max_rank) return;

//        utils::print(m_conn.ann());
//        utils::print(m_conn.cre());
//        utils::print(m_conn.com());

        for (size_t rank=1ul; rank<=m_max_rank; ++rank) {
            auto& rdm = *m_rdms[rank];
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
