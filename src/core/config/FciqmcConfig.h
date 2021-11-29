//
// Created by rja on 25/06/2021.
//

#ifndef M7_FCIQMCCONFIG_H
#define M7_FCIQMCCONFIG_H


#include "Parameters.h"

namespace fciqmc_config {

    struct HashMapping : config::Section {
        config::Param<double> m_remap_ratio;
        config::Param<size_t> m_remap_nlookup;

        explicit HashMapping(config::Group *parent);
    };

    struct Buffers : config::Section {
        config::Param<double> m_store_fac_init;
        config::Param<double> m_store_exp_fac;
        config::Param<double> m_comm_fac_init;
        config::Param<double> m_comm_exp_fac;

        explicit Buffers(config::Group *parent);
    };

    struct Archive : config::Section {
        config::Param<std::string> m_load_path;
        config::Param<std::string> m_save_path;
        config::Param<std::string> m_chkpt_path;
        /**
         * the period options are not mutually exclusive, the Archive class will make sure that a checkpoint file does
         * not get dumped twice on the same MC cycle
         */
        config::Param<size_t> m_period;
        config::Param<size_t> m_period_mins;

        explicit Archive(config::Group *parent);

        void verify() override;

        bool do_load() const {
            return !m_load_path.get().empty();
        }

        bool do_save() const {
            return !m_save_path.get().empty();
        }

        bool do_chkpts() const {
            return !m_chkpt_path.get().empty() && (m_period_mins || m_period);
        }
    };

    struct Archivable : config::Section {
        config::Param<bool> m_load;
        config::Param<bool> m_save;
        config::Param<bool> m_chkpt;

        explicit Archivable(config::Group *parent);
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

    struct MbfDef : config::Section {
        config::Param<std::vector<defs::inds>> m_frm;
        config::Param<std::vector<defs::inds>> m_bos;
        //config::Param<std::vector<defs::inds>> m_csf;

        explicit MbfDef(config::Group *parent, std::string name);
    };

    struct Reference : config::Section {
        config::Param<std::string> m_frm_onv_init;
        config::Param<defs::inds> m_bos_onv_init;
        MbfDef m_mbf_init;
        config::Param<bool> m_init_mbf_neel;
        config::Param<double> m_redef_thresh;

        explicit Reference(config::Group *parent);
    };

    struct Wavefunction : config::Section {
        config::Param<double> m_nw_init;
        config::Param<size_t> m_nroot;
        config::Param<long> m_spin_restrict;
        Buffers m_buffers;
        HashMapping m_hash_mapping;
        Archivable m_archivable;
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
        config::Param<double> m_l1_fraction_cutoff;
        config::Param<size_t> m_delay;

        explicit Semistochastic(config::Group *parent);

        void verify() override;
    };

    struct Fcidump : config::Section {
        config::Param<std::string> m_path;
        config::Param<std::string> m_eb_path;
        config::Param<std::string> m_bos_path;
        config::Param<bool> m_spin_major;

        explicit Fcidump(config::Group *parent);
    };

    struct Stats : config::Section {
        config::Param<std::string> m_path;
        config::Param<size_t> m_period;
        config::Param<bool> m_parallel;

        explicit Stats(config::Group *parent);

        void verify() override;
    };

    struct SpfWeightedTwf : config::Section {
        config::Param<double> m_fermion_fac;
        config::Param<double> m_boson_fac;

        explicit SpfWeightedTwf(config::Group *parent);

        void verify() override;
    };

    struct Bilinears : config::Section {
        config::Param<std::vector<std::string>> m_ranks;
        Buffers m_buffers;
        HashMapping m_hash_mapping;
        LoadBalancing m_load_balancing;
        Archivable m_archivable;

        explicit Bilinears(config::Group *parent, std::string name, std::string description);

        void verify() override;
    };

    struct Rdms : Bilinears {
        config::Param<bool> m_explicit_ref_conns;

        explicit Rdms(config::Group *parent, std::string name, std::string description);
    };

    struct SpecMoms : Bilinears {
        config::Param<bool> m_stochastic;
        config::Param<double> m_nattempt_per_walker;

        explicit SpecMoms(config::Group *parent, std::string name, std::string description);
    };

    struct InstEsts : config::Section {
        config::Param<bool> m_spf_uniform_twf;
        SpfWeightedTwf m_spf_weighted_twf;

        explicit InstEsts(config::Group *parent);
    };

    struct RefExcits : config::Section {
        config::Param<size_t> m_max_exlvl;
        Buffers m_buffers;
        Archivable m_archivable;

        explicit RefExcits(config::Group *parent);
    };

    struct AvEsts : config::Section {
        config::Param<size_t> m_delay;
        config::Param<size_t> m_ncycle;
        config::Param<size_t> m_stats_period;
        config::Param<std::string> m_stats_path;
        Rdms m_rdm;
        SpecMoms m_spec_mom;
        RefExcits m_ref_excits;

        explicit AvEsts(config::Group *parent);

        bool any_bilinears() const {
            return !(m_rdm.m_ranks.get().empty() && m_spec_mom.m_ranks.get().empty());
        }
    };

    struct Hamiltonian : config::Section {
        Fcidump m_fcidump;
        config::Param<int> m_charge;
        config::Param<bool> m_elecs;
        config::Param<size_t> m_nboson_max;

        Hamiltonian(config::Group *parent);

        void verify() override;
    };

    struct Propagator : config::Section {
        config::Param<size_t> m_ncycle;
        config::Param<bool> m_stochastic;
        config::Param<std::string> m_excit_gen;
        config::Param<double> m_nw_target;
        config::Param<defs::wf_comp_t> m_max_bloom;
        config::Param<defs::wf_comp_t> m_nadd;
        config::Param<double> m_tau_init;
        config::Param<double> m_tau_min;
        config::Param<double> m_tau_max;
        config::Param<bool> m_static_tau;
        config::Param<bool> m_static_probs;
        config::Param<double> m_min_spawn_mag;
        config::Param<double> m_min_death_mag;
        config::Param<double> m_min_exlvl_prob;
        config::Param<std::vector<double>> m_exlvl_probs_init;
        config::Param<size_t> m_ndraw_min_for_dynamic;
        config::Param<size_t> m_period;
        Semistochastic m_semistochastic;

        explicit Propagator(config::Group *parent);

        void verify();
    };

    struct Document : config::Document {
        Prng m_prng;
        Archive m_archive;
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
