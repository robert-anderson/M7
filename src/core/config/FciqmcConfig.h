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

    struct Serialization : config::Section {
        config::Param<std::string> m_save_path;
        config::Param<std::string> m_load_path;

        explicit Serialization(config::Group *parent, std::string path_default="");
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
        Serialization m_serialization;
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

    struct AvCoeffs : config::Section {
        config::Param<size_t> m_max_exlvl;
        Buffers m_buffers;

        explicit AvCoeffs(config::Group *parent) :
                config::Section(parent, "av_coeffs",
                                "options relating to the average coefficients of excitations of the reference"),
                m_max_exlvl(this, "max_exlvl", 0ul,
                            "maximum excitation level from the reference for which to accumulate average coefficients"),
                m_buffers(this) {}
    };

    struct FermionRdm : config::Section {
        config::Param<size_t> m_rank;
        config::Param<bool> m_mixed_estimator;
        Buffers m_buffers;
        Serialization m_serialization;
        LoadBalancing m_load_balancing;

        explicit FermionRdm(config::Group *parent) :
                config::Section(parent, "fermion_rdm",
                                "options relating to the accumulation and sampling of fermion RDM elements"),
                m_rank(this, "rank", 0ul, "Rank of fermion RDM to accumulate"),
                m_mixed_estimator(this, "mixed_estimator", false, "replace one instance of the wavefunction in the bilinear RDM definition with an SPF TWF"),
                m_buffers(this), m_serialization(this), m_load_balancing(this)
                {}
    };

    struct Observables : config::Section {
        config::Param<size_t> m_output_period;
        config::Param<size_t> m_delay;
        config::Param<size_t> m_ncycle;
        config::Param<bool> m_spf_uniform_twf;
        SpfWeightedTwf m_spf_weighted_twf;
        AvCoeffs m_av_coeffs;
        FermionRdm m_fermion_rdm;

        explicit Observables(config::Group *parent);
    };

    struct Hamiltonian : config::Section {
        Fcidump m_fcidump;
        config::Param<defs::ham_t> m_boson_frequency;
        config::Param<defs::ham_t> m_boson_coupling;
        config::Param<defs::ham_t> m_nboson_max;

        Hamiltonian(config::Group *parent) :
                config::Section(parent, "hamiltonian", "options relating to the Hamiltonian operator"),
                m_fcidump(this),
                m_boson_frequency(this, "boson_frequency", 0.0,
                                  "frequency of onsite boson modes for Hubbard-Holstein model"),
                m_boson_coupling(this, "boson_coupling", 0.0,
                                 "coupling of onsite boson modes for Hubbard-Holstein model"),
                m_nboson_max(this, "nboson_max", 0ul, "maximum allowed occupation of bosonic modes") {}

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
        Observables m_observables;

        explicit Document(const yaml::File *file= nullptr);

        void verify();
    };
}


#endif //M7_FCIQMCCONFIG_H
