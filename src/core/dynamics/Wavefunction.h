//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_WAVEFUNCTION_H
#define M7_WAVEFUNCTION_H


#include <src/core/fermion/Determinant.h>
#include <src/core/io/Options.h>
#include "WalkerList.h"
#include "SpawnList.h"
#include "Propagator.h"
#include "src/core/parallel/Distributed.h"
#include "DeterministicSubspace.h"
#include "Reference.h"

class FciqmcCalculation;

class Wavefunction {
    FciqmcCalculation *m_fciqmc;
    const Options &m_input;
    Reference m_reference;
    const std::unique_ptr<Propagator> &m_prop;
    std::unique_ptr<DeterministicSubspace> m_detsub = nullptr;

    WalkerList m_data;
    SpawnList m_send, m_recv;

    Hybrid<defs::wf_t> m_aborted_weight;
    Hybrid<int64_t> m_ninitiator;

    bool m_in_semistochastic_epoch = false;

public:

    /*
     * Square norm is sum_i(|w_i|^2)
     */
    Hybrid<defs::wf_comp_t> m_square_norm;
    Hybrid<defs::wf_comp_t> m_d_square_norm;
    /*
     * Walker number is sum_i(|w_i|)
     */
    Distributed<defs::wf_comp_t> m_nwalker;
    Hybrid<defs::wf_comp_t> m_d_nwalker;
    defs::wf_comp_t m_nwalker_growth_rate;

    /*
     * number of occupied determinants, the principal variable is distributed, since
     * rank-resolved data could inform load balancing
     */
    Distributed<size_t> m_nocc_det;
    Hybrid<int> m_d_nocc_det;


    Wavefunction(FciqmcCalculation *fciqmc);

    ~Wavefunction(){
        std::cout << "# initiators: " << m_data.verify_ninitiator(m_input.nadd_initiator)<< std::endl;
    }

    /**
     * Effect arithmetic updates on member variables (no communication)
     */
    void update(const size_t& icycle);

    void propagate();

    void communicate();

    void consolidate_incoming_weight() {
        // TODO
    }

    void annihilate();

    /**
     * Perform the necessary thread and MPI reductions on members
     */
    void synchronize();

    void write_iter_stats(FciqmcStatsFile* stats_file);

private:
    void annihilate_row(const size_t &irow_recv);
};


#endif //M7_WAVEFUNCTION_H
