//
// Created by rja on 14/07/2021.
//

#ifndef M7_ANNIHILATOR_H
#define M7_ANNIHILATOR_H


#include <src/core/wavefunction/Wavefunction.h>
#include <src/core/wavefunction/Reference.h>

struct Annihilator {
    Wavefunction &m_wf;
    const Hamiltonian& m_ham;
    const References& m_refs;
    const defs::wf_comp_t m_nadd;
    const size_t& m_icycle;

private:
    SpawnTableRow m_work_row1;
    SpawnTableRow m_work_row2;
    const bool m_send_parents;
    const comparators::index_cmp_fn_t m_sort_cmp_fn;

    comparators::index_cmp_fn_t make_sort_cmp_fn();

public:

    Annihilator(Wavefunction &wf, const Hamiltonian& ham, const References& refs, const size_t& icycle, defs::wf_comp_t nadd);

    void sort_recv();

    void annihilate_row(const size_t &dst_ipart, const field::Mbf &dst_mbf, const defs::wf_t &delta_weight,
                        bool allow_initiation, const size_t &irow_store);

    void handle_dst_mbf_block(SpawnTableRow &block_start, SpawnTableRow &current,
                              const defs::wf_t &total_delta, const size_t& irow_store);

    void loop_over_dst_mbfs();

};


#endif //M7_ANNIHILATOR_H
