//
// Created by Robert J. Anderson on 25/06/2021.
//

#ifndef M7_CONF_H
#define M7_CONF_H


#include "ConfComponents.h"
#include "HamiltonianConf.h"
#include "M7_lib/basis/BasisData.h"
#include "ExcitgenConf.h"

namespace conf {
    
    using namespace conf_components;

    struct HashMapping : Section {
        static constexpr double c_default_remap_ratio = 2.0;
        static constexpr uint_t c_default_remap_nlookup = 500ul;
        Param<double> m_remap_ratio;
        Param<uint_t> m_remap_nlookup;

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
        Param<str_t> m_load_path;
        Param<str_t> m_save_path;
        Param<str_t> m_chkpt_path;
        /**
         * the period options are not mutually exclusive, the Archive class will make sure that a checkpoint file does
         * not get dumped twice on the same MC cycle
         */
        Param<uint_t> m_period;
        Param<uint_t> m_period_mins;

        explicit Archive(Group *parent);

        bool do_load() const {
            return !m_load_path.m_value.empty();
        }

        bool do_save() const {
            return !m_save_path.m_value.empty();
        }

        bool do_chkpts() const {
            return !m_chkpt_path.m_value.empty() && (m_period_mins || m_period);
        }

    protected:
        void validate_node_contents() override;
    };

    struct Archivable : Section {
        Param<bool> m_load;
        Param<bool> m_save;
        Param<bool> m_chkpt;

        explicit Archivable(Group *parent);
    };

    struct Prng : Section {
        Param<uint_t> m_seed;
        Param<uint_t> m_ngen_block;

        explicit Prng(Group *parent);
    };

    struct Distribution : Section {
        static constexpr uint_t c_default_nblock_per_rank = 6ul;
        static constexpr uint_t c_default_period = 10ul;
        static constexpr double c_default_imbalance_thresh = 0.05;
        Param<uint_t> m_nblock_per_rank;
        Param<uint_t> m_period;
        Param<double> m_imbalance_thresh;

        explicit Distribution(Group *parent);
    };

    struct MbfDef : Section {
        Param<v_t<uintv_t>> m_frm;
        Param<v_t<uintv_t>> m_bos;
        Param<bool> m_neel;
        //Param<v_t<uintv_t>> m_csf;

        explicit MbfDef(Group *parent, str_t name);
    };

    struct Reference : Section {
        MbfDef m_mbf_init;
        Param<double> m_redef_thresh;

        explicit Reference(Group *parent);
    };

    struct Basis : Section {
        Param<uint_t> m_bos_occ_cutoff;
        explicit Basis(Group *parent):
        Section(parent, "basis", "options relating to the single particle basis functions and subsets thereof"),
        m_bos_occ_cutoff(this, "bos_occ_cutoff", sys::bos::c_max_occ,
                         "maximum allowed occupation of each boson mode"){}

    protected:
        void validate_node_contents() override {
            REQUIRE_LE(m_bos_occ_cutoff, sys::bos::c_max_occ,
                       logging::format("given nboson_max exceeds limit of {}", sys::bos::c_max_occ));
        }
    };

    struct Particles : Section {
        Param<uint_t> m_nelec;
        Param<int> m_ms2;
        Param<uint_t> m_nboson;

        explicit Particles(Group *parent):
            Section(parent, "particles", "options relating to the particle number sector"),
            m_nelec(this, "nelec", 0ul, "number of electrons in the system (conserved)"),
            m_ms2(this, "ms2", sys::frm::c_undefined_ms2, "2*Ms sector in which the system is to be restricted (taken as reference hint if H does not conserve Sz"),
            m_nboson(this, "nboson", 0ul, "number of bosons in the system (taken as reference hint if H does not conserve boson number"){}
    };

    struct Wavefunction : Section {
        Param<double> m_nw_init;
        Param<uint_t> m_nroot;
        Param<bool> m_fci_init;
        Buffers m_buffers;
        HashMapping m_hash_mapping;
        Archivable m_archivable;
        Distribution m_distribution;

        explicit Wavefunction(Group *parent);
    };

    struct Shift : Section {
        Param<double> m_init;
        Param<double> m_damp;
        Param<double> m_target_damp;
        Param<uint_t> m_period;
        Param<uint_t> m_ncycle_av;
        Param<bool> m_jump;

        explicit Shift(Group *parent);

    };

    struct Semistochastic : Section {
        Param<uint_t> m_size;
        Param<double> m_l1_fraction_cutoff;
        Param<uint_t> m_delay;

        explicit Semistochastic(Group *parent);

    protected:
        void validate_node_contents() override;
    };

    struct Stats : Section {
        Param<str_t> m_path;
        Param<uint_t> m_period;
        Param<bool> m_parallel;

        explicit Stats(Group *parent);

    protected:
        void validate_node_contents() override;
    };

    struct SpfWeightedTwf : Section {
        Param<double> m_fermion_fac;
        Param<double> m_boson_fac;

        explicit SpfWeightedTwf(Group *parent);

    protected:
        void validate_node_contents() override;
    };

    struct Bilinears : Section {
        Param<strv_t> m_ranks;
        Buffers m_buffers;
        HashMapping m_hash_mapping;
        Distribution m_distribution;
        Archivable m_archivable;

        explicit Bilinears(Group *parent, str_t name, str_t description);

    protected:
        void validate_node_contents() override;
    };

    struct Fock4rdm : Section {
        Param<str_t> m_fock_path;
        explicit Fock4rdm(Group* parent);
    };

    struct Rdms : Bilinears {
        Param<bool> m_spinfree;
        Fock4rdm m_fock_4rdm;

        explicit Rdms(Group *parent, str_t name, str_t description);
    };

    struct SpecMoms : Bilinears {
        Param<bool> m_stochastic;
        Param<double> m_nattempt_per_walker;

        explicit SpecMoms(Group *parent, str_t name, str_t description);
    };

    struct InstEsts : Section {
        Param<bool> m_spf_uniform_twf;
        SpfWeightedTwf m_spf_weighted_twf;

        explicit InstEsts(Group *parent);
    };

    struct HfExcits : Section {
        Param<uint_t> m_max_exlvl;
        Buffers m_buffers;
        Archivable m_archivable;

        explicit HfExcits(Group *parent);
    };

    struct Mae : Section {
        Param<uint_t> m_delay;
        Param<uint_t> m_ncycle;
        Param<uint_t> m_stats_period;
        Param<str_t> m_stats_path;
        Rdms m_rdm;
        SpecMoms m_spec_mom;
        HfExcits m_hf_excits;

        explicit Mae(Group *parent);

        bool any_bilinears() const {
            return !(m_rdm.m_ranks.m_value.empty() && m_spec_mom.m_ranks.m_value.empty());
        }
    };


//    struct FrmExcitGens : Section {
//        SingleChoice<str_t> m_general_doubles;
//        SingleChoice<str_t> m_general_singles;
//        SingleChoice<str_t> m_hubbard;
//        SingleChoice<str_t> m_heisenberg;
//        FrmExcitGens(Group *parent): Section(parent, "excitgen",
//             "options relating to fermion hamiltonian excitation generation"),
//             m_general_doubles(this, "general_doubles", {"pchb"},
//                               "approach for doubles generation in the case that fermion hamiltonian is defined by an FCIDUMP"),
//                                     m_general_singles(this, "general_singles", {"uniform"},
//                                                       "approach for singles generation in the case that fermion hamiltonian is defined by an FCIDUMP"),
//                                     m_hubbard(this, "hubbard", {"uniform", ""},
//                                                       "approach for singles generation in the case that fermion hamiltonian is defined by an FCIDUMP"),
//                                     m_general_singles(this, "general_singles", {"uniform"},
//                                                       "approach for singles generation in the case that fermion hamiltonian is defined by an FCIDUMP"),
//
//
//    };


//    struct ExcitGens : Section {
//        SingleChoice<str_t> m_frm_doubles;
//        ExcitGens(Group *parent): Section(parent, "excitgen", "options relating to excitation generation")
//        m_frm_doubles(){}
//    };

    struct GuideWavefunction : Selection {
        struct ExpFac : Section {
            Param<double> m_fac;
            ExpFac(GuideWavefunction* parent, str_t name, str_t description);

        protected:
            void validate_node_contents() override;
        };

        struct GutzwillerLike : ExpFac {
            GutzwillerLike(GuideWavefunction* parent);
        };
        struct SuppressMultiOcc : ExpFac {
            SuppressMultiOcc(GuideWavefunction* parent);
        };

        GutzwillerLike m_gutzwiller_like;
        SuppressMultiOcc m_suppress_multi_occ;
        GuideWavefunction(Group *parent, const str_t& name);
    };


    struct Propagator : Section {
        Param<uint_t> m_ncycle;
        Param<bool> m_stochastic;
        conf::ExcitGen m_excit_gen;
        Param<double> m_nw_target;
        Param<wf_comp_t> m_max_bloom;
        Param<wf_comp_t> m_nadd;
        Param<double> m_tau_init;
        Param<double> m_tau_min;
        Param<double> m_tau_max;
        Param<bool> m_static_tau;
        Param<bool> m_static_probs;
        Param<double> m_min_spawn_mag;
        Param<double> m_min_death_mag;
        Param<double> m_min_exlvl_prob;
        Param<v_t<double>> m_exlvl_probs_init;
        Param<uint_t> m_ndraw_min_for_dynamic;
        Param<uint_t> m_period;
        GuideWavefunction m_imp_samp_guide;
        Semistochastic m_semistochastic;

        explicit Propagator(Group *parent);

    protected:
        void validate_node_contents() override;
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
        Mae m_av_ests;

        explicit Document(const str_t& fname="");

    protected:
        void validate_node_contents() override;
    };
}


#endif //M7_CONF_H
