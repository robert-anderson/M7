//
// Created by rja on 12/11/2020.
//

#include <src/core/table/BufferedTable.h>
#include <src/core/field/Elements.h>
#include "gtest/gtest.h"
#include "src/core/table/MappedTable.h"
#include "src/core/excitgen/BosonCouplingSamplers.h"

struct ExcitGenTestTable : MappedTable<fields::FermiBosOnv> {
    fields::FermiBosOnv m_onv;
    fields::Number<size_t> m_frequency;
    fields::Number<defs::prob_t> m_weight;
    ExcitGenTestTable(size_t nsite):
    MappedTable<fields::FermiBosOnv>(m_onv, 1000),
    m_onv(this, {nsite, nsite}, "occupation number vector"),
    m_frequency(this, "number of times the ONV was drawn"),
    m_weight(this, "cumulative reciprocal probability")
    {}
};

TEST(BosonCouplingSamplers, SingleOnvTest){
    const size_t nsite = 6;
    const size_t nboson_max = 3;
    BufferedTable<ExcitGenTestTable> bt(nsite);

    const size_t ndraw = 10000;

    elements::FermiBosOnv src_onv(nsite, nsite);
    elements::FermiBosOnv dst_onv(nsite, nsite);
    src_onv = {{0, 4, 6, 11}, {1, 0, 0, 1, 3, 2}};

    std::cout << src_onv.to_string() << std::endl;

    PRNG prng = PRNG(18, 1e4);
    BosonCouplings bc(nsite, nboson_max, 1.0, 0.5);
    BosonCouplingSamplers sampler(bc, nboson_max, prng);
    OccupiedOrbitals occ_orbs(src_onv.m_fonv);
    conn::AsFermiBosOnv aconn(src_onv);
    defs::prob_t prob;
    defs::ham_t helem;
    for (size_t idraw=0ul; idraw<ndraw; ++idraw){
        sampler.draw_single(src_onv, dst_onv, occ_orbs, prob, helem, aconn);
        auto irow = *bt[dst_onv];
        if (irow==~0ul){
            bt.expand(1);
            irow = bt.insert(dst_onv);
        }
        bt.m_frequency(irow)++;
        bt.m_weight(irow)+=1.0/prob;
    }
    /*
     * TODO: James
     *  ASSERT that the m_weights are approximately the same,
     *  would be better to refactor the above into a function that can be called
     *  multiple times with different ndraw values, and check that the variance
     *  decreases as expected
     *  ASSERT that all expected excitations are generated at least once (by counting,
     *  i.e. check the m_hwm member of bt is correct)
     */
    for(size_t irow = 0ul; irow < bt.m_hwm; ++irow){
        std::cout << bt.m_weight(irow) << std::endl;
    }
}