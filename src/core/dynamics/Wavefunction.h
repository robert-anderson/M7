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

class FciqmcCalculation;

class Wavefunction {
    FciqmcCalculation *m_fciqmc;
    const Options &m_input;
    Determinant& m_reference;
    const std::unique_ptr<Propagator> &m_prop;

    WalkerList m_data;
    SpawnList m_send, m_recv;

    defs::ham_t m_ref_proj_energy_num;
    defs::wf_t m_reference_weight;
    size_t m_reference_row;

    defs::wf_t m_aborted_weight;
    int m_ninitiator = 0;

public:

    /*
     * Square norm is sum_i(|w_i|^2)
     */
    defs::wf_comp_t m_square_norm;
    defs::wf_comp_t m_delta_square_norm;
    /*
     * Walker number is sum_i(|w_i|)
     */
    defs::wf_comp_t m_nw;
    defs::wf_comp_t m_delta_nw;
    defs::wf_comp_t m_nw_growth_rate = 0;

    size_t m_noccupied_determinant;


    Wavefunction(FciqmcCalculation *fciqmc);

    ~Wavefunction(){
        std::cout << "# initiators: " << m_data.verify_ninitiator(m_input.nadd_initiator)<< std::endl;
    }

    void propagate();

    void communicate();

    void consolidate_incoming_weight() {
        // TODO
    }

    void annihilate();

    void write_iter_stats(FciqmcStatsFile &stats_file);

private:
    void annihilate_row(const size_t &irow_recv, defs::wf_comp_t &aborted_weight,
                        defs::wf_comp_t &delta_square_norm, defs::wf_comp_t &delta_nw);
};


#endif //M7_WAVEFUNCTION_H
