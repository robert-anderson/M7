//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_WAVEFUNCTION_H
#define M7_WAVEFUNCTION_H


#include <src/core/fermion/Determinant.h>
#include <src/core/io/InputOptions.h>
#include "WalkerList.h"
#include "SpawnList.h"
#include "Propagator.h"

class Wavefunction {
    WalkerList m_data;
    SpawnList m_send, m_recv;

    Determinant m_reference;
    defs::ham_t m_ref_proj_energy_num;
    defs::wf_t m_reference_weight;
    size_t m_reference_row;

    defs::wf_t m_aborted_weight;
    int m_ninitiator = 0;

public:

    const InputOptions &m_input;
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


    Wavefunction(const InputOptions &input, const std::unique_ptr<Propagator> &propagator,
                 const Determinant &reference);

    void propagate(std::unique_ptr<Propagator> &propagator);

    void communicate();

    void consolidate_incoming_weight() {
        // TODO
    }

    void annihilate(const std::unique_ptr<Propagator> &propagator);

    void write_iter_stats(FciqmcStatsFile &stats_file);

    defs::ham_comp_t norm() const;

private:
    void annihilate_row(const size_t& irow_recv, const std::unique_ptr<Propagator> &propagator,
                        Connection& connection, defs::wf_comp_t& aborted_weight, defs::wf_comp_t &delta_square_norm);
};


#endif //M7_WAVEFUNCTION_H
