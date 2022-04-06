//
// Created by rja on 06/04/2022.
//

#ifndef M7_EXCITGENGROUP2_H
#define M7_EXCITGENGROUP2_H

#include "M7_lib/hamiltonian/Hamiltonian.h"

using namespace exsig_utils;

class ExcitGenGroup2 {
    /**
     * forward list containing all excitation generators for all hamiltonian terms
     */
    ExcitGen2::excit_gen_list_t m_list;
    /**
     * one dynamically-constructed excitation generator can be used to generate excitations of many different exsigs
     */
    std::array<std::vector<ExcitGen2*>, defs::nexsig> m_exsigs_to_exgen_lists = {};

public:
    ExcitGenGroup2(const Hamiltonian &ham, const fciqmc_config::Propagator &opts, PRNG &prng) {
        ExcitGen2::excit_gen_list_t list;
        list = ham.m_frm->make_excit_gens(prng, opts);
        m_list.merge(list);
        list = ham.m_frmbos->make_excit_gens(prng, opts);
        m_list.merge(list);
        list = ham.m_bos->make_excit_gens(prng, opts);
        m_list.merge(list);
        /*
         * exsigs have a many-to-many relationship with excitation generators.
         */
        for (const auto& excit_gen: m_list) {
            for (auto& exsig: excit_gen->m_exsigs){
                m_exsigs_to_exgen_lists[exsig].push_back(excit_gen.get());
            }
        }
    }

};


#endif //M7_EXCITGENGROUP2_H
