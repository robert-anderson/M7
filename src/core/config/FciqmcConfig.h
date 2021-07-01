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

        explicit Buffers(config::Group *parent);
    };

    struct Io : config::Section {
        config::Param<std::string> m_save_path;
        config::Param<std::string> m_load_path;

        explicit Io(config::Group *parent, std::string path_default = "");
    };

    struct Prng : config::Section {
        config::Param<size_t> m_seed;
        config::Param<size_t> m_ngen_block;

        explicit Prng(config::Group *parent);
    };

    struct LoadBalancing : config::Section {
        config::Param<size_t> m_nblock_per_rank;
        config::Param<size_t> m_period;
        config::Param<double> m_acceptable_imbalance;
        config::Param<size_t> m_nnull_updates_deactivate;

        explicit LoadBalancing(config::Group *parent);
    };

    struct Reference : config::Section {
        config::Param<std::vector<std::string>> m_init_onv;
        config::Param<double> m_redef_thresh;

        explicit Reference(config::Group *parent);
    };

    struct Wavefunction : config::Section {
        config::Param<double> m_nw_init;
        config::Param<size_t> m_nroot;
        config::Param<bool> m_replicate;
        config::Param<long> m_spin_restrict;
        Buffers m_buffers;
        Io m_io;
        LoadBalancing m_load_balancing;

        explicit Wavefunction(config::Group *parent);
    };

    struct Reweight : config::Section {
        config::Param<size_t> m_ncycle;
        config::Param<size_t> m_delay;

        explicit Reweight(config::Group *parent);
    };

    struct Shift : config::Section {
        config::Param<double> m_init;
        config::Param<double> m_damp;
        config::Param<size_t> m_period;
        config::Param<size_t> m_ncycle_av;
        config::Param<bool> m_jump;
        Reweight m_reweight;

        explicit Shift(config::Group *parent);

    };

    struct Semistochastic : config::Section {
        config::Param<size_t> m_size;
        config::Param<size_t> m_delay;

        explicit Semistochastic(config::Group *parent);
    };

    struct Fcidump : config::Section {
        config::Param<std::string> m_path;
        config::Param<bool> m_spin_major;

        explicit Fcidump(config::Group *parent);
    };

    struct Stats : config::Section {
        config::Param<std::string> m_path;
        config::Param<bool> m_parallel;

        explicit Stats(config::Group *parent);
    };

    struct SpfWeightedTwf : config::Section {
        config::Param<double> m_fermion_fac;
        config::Param<double> m_boson_fac;

        explicit SpfWeightedTwf(config::Group *parent);

        void verify() override;
    };

    struct FermionRdm : config::Section {
        config::Param<size_t> m_rank;
        config::Param<bool> m_mixed_estimator;
        Buffers m_buffers;
        LoadBalancing m_load_balancing;

        explicit FermionRdm(config::Group *parent);
    };

    struct InstEsts : config::Section {
        config::Param<bool> m_spf_uniform_twf;
        SpfWeightedTwf m_spf_weighted_twf;

        explicit InstEsts(config::Group *parent);
    };

    struct RefExcits : config::Section {
        config::Param<size_t> m_max_exlvl;
        Buffers m_buffers;

        explicit RefExcits(config::Group *parent);
    };

    struct PeriodicOutput : config::Section {
        config::Param<size_t> m_period;
        config::Param<std::string> m_path;
        config::Param<bool> m_clobber;
        explicit PeriodicOutput(config::Group *parent);

        void verify() override;
    };

    struct AvEsts : config::Section {
        config::Param<size_t> m_delay;
        config::Param<size_t> m_ncycle;
        FermionRdm m_fermion_rdm;
        RefExcits m_ref_excits;
        Io m_io;
        PeriodicOutput m_periodic_output;

        explicit AvEsts(config::Group *parent);
    };

    struct Hamiltonian : config::Section {
        Fcidump m_fcidump;
        config::Param<defs::ham_t> m_boson_frequency;
        config::Param<defs::ham_t> m_boson_coupling;
        config::Param<defs::ham_t> m_nboson_max;

        Hamiltonian(config::Group *parent);

        void verify() override;
    };

    struct Propagator : config::Section {
        config::Param<size_t> m_ncycle;
        config::Param<bool> m_exact;
        config::Param<std::string> m_excit_gen;
        config::Param<double> m_nw_target;
        config::Param<double> m_max_bloom;
        config::Param<double> m_nadd;
        config::Param<double> m_tau_init;
        config::Param<bool> m_static_tau;
        config::Param<double> m_min_spawn_mag;
        config::Param<double> m_min_death_mag;
        config::Param<double> m_min_excit_class_prob;
        config::Param<double> m_psingle_init;
        config::Param<size_t> m_nenough_spawns_for_dynamic_tau;
        config::Param<bool> m_consolidate_spawns;
        Semistochastic m_semistochastic;

        explicit Propagator(config::Group *parent);

        void verify();
    };

    struct Document : config::Document {
        Prng m_prng;
        Wavefunction m_wavefunction;
        Reference m_reference;
        Shift m_shift;
        Propagator m_propagator;
        Hamiltonian m_hamiltonian;
        Stats m_stats;
        InstEsts m_inst_ests;
        AvEsts m_av_ests;

        explicit Document(const yaml::File *file = nullptr);

        void verify();
    };
}


#endif //M7_FCIQMCCONFIG_H
