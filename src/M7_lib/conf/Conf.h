//
// Created by Robert J. Anderson on 25/06/2021.
//

#ifndef M7_CONF_H
#define M7_CONF_H


#include "ConfComponents.h"
#include "HamiltonianConf.h"

namespace conf {
    
    using namespace conf_components;

    struct HashMapping : Section {
        Param<double> m_remap_ratio;
        Param<size_t> m_remap_nlookup;

        explicit HashMapping(Group *parent);
    };

    struct Buffers : Section {
        Param<double> m_store_fac_init;
        Param<double> m_store_exp_fac;
        Param<double> m_comm_fac_init;
        Param<double> m_comm_exp_fac;

        explicit Buffers(Group *parent);
    };

    struct Archive : Section {
        Param<std::string> m_load_path;
        Param<std::string> m_save_path;
        Param<std::string> m_chkpt_path;
        /**
         * the period options are not mutually exclusive, the Archive class will make sure that a checkpoint file does
         * not get dumped twice on the same MC cycle
         */
        Param<size_t> m_period;
        Param<size_t> m_period_mins;

        explicit Archive(Group *parent);

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

    struct Archivable : Section {
        Param<bool> m_load;
        Param<bool> m_save;
        Param<bool> m_chkpt;

        explicit Archivable(Group *parent);
    };

    struct Prng : Section {
        Param<size_t> m_seed;
        Param<size_t> m_ngen_block;

        explicit Prng(Group *parent);
    };

    struct LoadBalancing : Section {
        Param<size_t> m_nblock_per_rank;
        Param<size_t> m_period;
        Param<double> m_acceptable_imbalance;
        Param<size_t> m_nnull_updates_deactivate;

        explicit LoadBalancing(Group *parent);
    };

    struct MbfDef : Section {
        Param<std::vector<defs::inds>> m_frm;
        Param<std::vector<defs::inds>> m_bos;
        Param<bool> m_neel;
        //Param<std::vector<defs::inds>> m_csf;

        explicit MbfDef(Group *parent, std::string name);

        bool enabled() const override;
    };

    struct Reference : Section {
        MbfDef m_mbf_init;
        Param<double> m_redef_thresh;

        explicit Reference(Group *parent);
    };

    struct Basis : Section {
        Param<size_t> m_bos_occ_cutoff;
        explicit Basis(Group *parent):
        Section(parent, "basis", "options relating to the single particle basis functions and subsets thereof"),
        m_bos_occ_cutoff(this, "bos_occ_cutoff", defs::max_bos_occ, "maximum allowed occupation of each boson mode"){}

        void verify() override {
            REQUIRE_LE(m_bos_occ_cutoff, defs::max_bos_occ, log::format("given nboson_max exceeds limit of {}", defs::max_bos_occ));
        }
    };

    struct Particles : Section {
        Param<size_t> m_nelec;
        Param<int> m_ms2;
        Param<size_t> m_nboson;

        explicit Particles(Group *parent):
            Section(parent, "particles", "options relating to the particle number sector"),
            m_nelec(this, "nelec", 0ul, "number of electrons in the system (conserved)"),
            m_ms2(this, "ms2", defs::undefined_ms2, "2*Ms sector in which the system is to be restricted (taken as reference hint if H does not conserve Sz"),
            m_nboson(this, "nboson", 0ul, "number of bosons in the system (taken as reference hint if H does not conserve boson number"){}
    };

    struct Wavefunction : Section {
        Param<double> m_nw_init;
        Param<size_t> m_nroot;
        Buffers m_buffers;
        HashMapping m_hash_mapping;
        Archivable m_archivable;
        LoadBalancing m_load_balancing;

        explicit Wavefunction(Group *parent);
    };

    struct Reweight : Section {
        Param<size_t> m_ncycle;
        Param<size_t> m_delay;

        explicit Reweight(Group *parent);
    };

    struct Shift : Section {
        Param<double> m_init;
        Param<double> m_damp;
        Param<size_t> m_period;
        Param<size_t> m_ncycle_av;
        Param<bool> m_jump;
        Reweight m_reweight;

        explicit Shift(Group *parent);

    };

    struct Semistochastic : Section {
        Param<size_t> m_size;
        Param<double> m_l1_fraction_cutoff;
        Param<size_t> m_delay;

        explicit Semistochastic(Group *parent);

        void verify() override;
    };

    struct Stats : Section {
        Param<std::string> m_path;
        Param<size_t> m_period;
        Param<bool> m_parallel;

        explicit Stats(Group *parent);

        void verify() override;
    };

    struct SpfWeightedTwf : Section {
        Param<double> m_fermion_fac;
        Param<double> m_boson_fac;

        explicit SpfWeightedTwf(Group *parent);

        void verify() override;
    };

    struct Bilinears : Section {
        Param<std::vector<std::string>> m_ranks;
        Buffers m_buffers;
        HashMapping m_hash_mapping;
        LoadBalancing m_load_balancing;
        Archivable m_archivable;

        explicit Bilinears(Group *parent, std::string name, std::string description);

        void verify() override;
    };

    struct Rdms : Bilinears {
        Param<bool> m_explicit_ref_conns;

        explicit Rdms(Group *parent, std::string name, std::string description);
    };

    struct SpecMoms : Bilinears {
        Param<bool> m_stochastic;
        Param<double> m_nattempt_per_walker;

        explicit SpecMoms(Group *parent, std::string name, std::string description);
    };

    struct InstEsts : Section {
        Param<bool> m_spf_uniform_twf;
        SpfWeightedTwf m_spf_weighted_twf;

        explicit InstEsts(Group *parent);
    };

    struct RefExcits : Section {
        Param<size_t> m_max_exlvl;
        Buffers m_buffers;
        Archivable m_archivable;

        explicit RefExcits(Group *parent);
    };

    struct AvEsts : Section {
        Param<size_t> m_delay;
        Param<size_t> m_ncycle;
        Param<size_t> m_stats_period;
        Param<std::string> m_stats_path;
        Rdms m_rdm;
        SpecMoms m_spec_mom;
        RefExcits m_ref_excits;

        explicit AvEsts(Group *parent);

        bool any_bilinears() const {
            return !(m_rdm.m_ranks.get().empty() && m_spec_mom.m_ranks.get().empty());
        }
    };

    struct Propagator : Section {
        Param<size_t> m_ncycle;
        Param<bool> m_stochastic;
        Param<std::string> m_excit_gen;
        Param<double> m_nw_target;
        Param<defs::wf_comp_t> m_max_bloom;
        Param<defs::wf_comp_t> m_nadd;
        Param<double> m_tau_init;
        Param<double> m_tau_min;
        Param<double> m_tau_max;
        Param<bool> m_static_tau;
        Param<bool> m_static_probs;
        Param<double> m_min_spawn_mag;
        Param<double> m_min_death_mag;
        Param<double> m_min_exlvl_prob;
        Param<std::vector<double>> m_exlvl_probs_init;
        Param<size_t> m_ndraw_min_for_dynamic;
        Param<size_t> m_period;
        Param<double> m_imp_samp_exp;
        Semistochastic m_semistochastic;

        explicit Propagator(Group *parent);

        void verify();
    };

    struct Document : conf_components::Document {
        Prng m_prng;
        Archive m_archive;
        Basis m_basis;
        Particles m_particles;
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


#endif //M7_CONF_H
