//
// Created by rja on 25/06/2021.
//

#ifndef M7_FCIQMCCONFIG_H
#define M7_FCIQMCCONFIG_H


#include "Parameters.h"

namespace fciqmc_config {

    struct Buffers : config::Section {
        config::Param<double> m_store_fac_init;
        config::Param<double> m_store_exp_fac;
        config::Param<double> m_comm_fac_init;
        config::Param<double> m_comm_exp_fac;
        Buffers(config::Group* parent);
    };

    struct Serialization : config::Section {
        config::Param<std::string> m_save_path;
        config::Param<std::string> m_load_path;
        Serialization(config::Group* parent) :
        config::Section(parent, "serialization", "options relating to filesystem save and load of structures in an M7 calculation"),
        m_save_path(this, "save_path", {}, "path to which the HDF5 file containing the structure should be saved"),
        m_load_path(this, "load_path", {}, "path from which the HDF5 file containing the structure should be loaded")
        {}
    };

    struct Prng : config::Section {
        config::Param<size_t> m_seed;
        config::Param<size_t> m_ngen_block;
        Prng(config::Group* parent);
    };

    struct LoadBalancing : config::Section {
        config::Param<size_t> m_nblock_per_rank;
        config::Param<size_t> m_period;
        config::Param<double> m_acceptable_imbalance;
        config::Param<size_t> m_nnull_updates_deactivate;
        LoadBalancing(config::Group* parent);
    };

    struct Wavefunction : config::Section {
        config::Param<double> m_nw_init;
        config::Param<size_t> m_nroot;
        config::Param<bool> m_replicate;
        config::Param<long> m_spin_restrict;
        config::Param<std::vector<std::string>> m_init_reference_onv;
        Serialization m_serialization;
        LoadBalancing m_load_balancing;
        Wavefunction(config::Group* parent);
    };

    struct Reweight : config::Section {
        config::Param<size_t> m_ncycle;
        config::Param<size_t> m_delay;
        Reweight(config::Group* parent);
    };

    struct Shift : config::Section {
        config::Param<double> m_init;
        config::Param<double> m_damp;
        config::Param<size_t> m_period;
        config::Param<size_t> m_ncycle_av;
        config::Param<bool> m_jump;
        Reweight m_reweight;
        Shift(config::Group* parent);

    };

    struct Propagator : config::Section {
        config::Param<bool> m_exact;
        config::Param<double> m_nw_target;
        config::Param<double> m_max_bloom;
        config::Param<double> m_nadd;
        config::Param<double> m_tau_init;
        config::Param<bool> m_static_tau;
        config::Param<double> m_min_spawn_mag;
        config::Param<double> m_min_death_mag;
        config::Param<bool> m_consolidate_spawns;

        Propagator(config::Group* parent);
    };

    struct Semistochastic : config::Section {

    };

    struct Document : config::Document {
        Prng m_prng;
        Wavefunction m_wavefunction;
        Shift m_shift;
        Propagator m_propagator;
        Document(const yaml::File* file);
    };
}


#endif //M7_FCIQMCCONFIG_H
