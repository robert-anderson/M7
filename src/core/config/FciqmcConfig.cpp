//
// Created by rja on 25/06/2021.
//

#include "FciqmcConfig.h"


fciqmc_config::HashMapping::HashMapping(config::Group *parent) :
        config::Section(parent, "hash_mapping", "options relating to the behavior of hash-mapped tables"),
        m_remap_ratio(this, "remap_ratio", 2.0,
                      "ratio of bucket-searching skips to total lookups required to trigger remapping with a larger number of buckets"),
        m_remap_nlookup(this, "remap_nlookup", 500ul,
                        "number of lookups required before remapping based on the skips/lookups ratio is considered") {}

fciqmc_config::Buffers::Buffers(config::Group *parent) :
        config::Section(parent, "buffers",
                        "options relating to the allocation and reallocation behavior of a Communicator"),
        m_store_fac_init(this, "store_fac_init", 1.0,
                         "a crude estimate for the number of rows ultimately required by the store buffer is computed, then that estimate is multiplied by this parameter to determine the initial row allocation"),
        m_store_exp_fac(this, "store_expand_fac", 0.5,
                        "additional number of rows that should be added to the store buffer's capacity as a fraction of the required number of additional rows"),
        m_comm_fac_init(this, "comm_fac_init", 1.0,
                        "a crude estimate for the number of rows per rank ultimately required by the send buffer is computed, then that estimate is multiplied by this parameter to determine the initial row allocation"),
        m_comm_exp_fac(this, "comm_expand_fac", 0.5,
                       "additional number of rows that should be added to the communicating buffers' capacities as a fraction of the required number of additional rows") {}


fciqmc_config::Archive::Archive(config::Group *parent) :
        config::Section(parent, "archive",
                        "options relating to archives: HDF5 binary storage of calculation data"),
        m_load_path(this, "load_path", "",
                    "path to the HDF5 file from which the calculation is to be restarted"),
        m_save_path(this, "save_path", "",
                    "path to the HDF5 file into which calculation data is to be dumped at the end of the calculation"),
        m_chkpt_path(this, "chkpt_path", "",
                     "path to the HDF5 file (or file name format) into which calculation data is to be dumped periodically during the calculation"),
        m_period(this, "period", 0ul, "number of MC cycles between checkpoint dumps"),
        m_period_mins(this, "period_mins", 0ul, "time in minutes between checkpoint dumps") {}

void fciqmc_config::Archive::verify() {
    if (static_cast<bool>(m_period) && static_cast<bool>(m_period_mins))
        log::warn("both cycle number and time periods are defined for checkpointing");

    auto &str = m_chkpt_path.get();
    size_t token_count = std::count(str.cbegin(), str.cend(), '{');
    REQUIRE_LE_ALL(token_count, 1ul, "checkpoint paths can have at most one {} token");
    if (token_count) {
        auto it_open = std::find(str.cbegin(), str.cend(), '{');
        auto it_close = std::find(str.cbegin(), str.cend(), '}');
        REQUIRE_EQ_ALL(std::distance(it_open, it_close), 1l,
                       "checkpoint path for periodic output should contain at most one {} formatting point");
        log::info(
                "formatting token found in path, successive checkpoints will not overwrite previous checkpoints from the same run");
    } else
        log::info(
                "formatting token not found in path, successive checkpoints will overwrite previous checkpoints from the same run");
}

fciqmc_config::Archivable::Archivable(config::Group *parent) :
        config::Section(parent, "archive",
                        "options relating to archiving behavior"),
        m_load(this, "load", false,
               "attempt to load the object from the archive at the beginning of the calculation"),
        m_save(this, "save", false,
               "save the object to the archive at the end of the calculation"),
        m_chkpt(this, "chkpt", false,
                "attempt to save the object to the checkpoint archive at each checkpointing period") {}


fciqmc_config::Prng::Prng(config::Group *parent) :
        config::Section(parent, "prng",
                        "options relating to the random number generator used in stochastic calculations"),
        m_seed(this, "seed", 123ul, "value with which to seed the mt19937 PRNG"),
        m_ngen_block(this, "ngen_block", 10000ul,
                     "size of the block of PRNGs generated each time the buffer is depleted") {}

fciqmc_config::LoadBalancing::LoadBalancing(config::Group *parent) :
        config::Section(parent, "load_balancing",
                        "options relating to the allocation of records among MPI ranks so as to share the workload more equally"),
        m_nblock_per_rank(this, "nblock_per_rank", 6ul, "number of rank allocation blocks to create per MPI rank"),
        m_period(this, "period", 10ul, "number of MC cycles between load-balancing block transactions"),
        m_acceptable_imbalance(this, "acceptable_imbalance", 0.05,
                               "fractional difference in the work figure of the busiest and laziest ranks at which the load balancing is deactivated after a given number of consecutive failed attempts"),
        m_nnull_updates_deactivate(this, "nnull_updates_deactivate", 20ul,
                                   "number of consecutive attempted updates which do not exceed the maximum acceptable imbalance required to meet the deactivation criterion") {}

fciqmc_config::Reference::Reference(config::Group *parent) :
        config::Section(parent, "reference", "options relating to the reference MBF"),
        m_init_mbf(this, "init_reference_mbf", {},
                   "string representations of the MBFs to use as the init references for each root. If replication is used, the init state will be used for both replicas of a root population."),
        m_redef_thresh(this, "redef_thresh", 10.0,
                       "when the highest-weighted non-reference MBF (the candidate) reaches this multiple of the weight on the reference, the candidate will be adopted as the new reference") {}

fciqmc_config::Wavefunction::Wavefunction(config::Group *parent) :
        config::Section(parent, "wavefunction",
                        "options relating to the storage and update of a distributed many-body wavefunction"),
        m_nw_init(this, "nw_init", 1ul, "L1 norm of the initial wavefunction"),
        m_nroot(this, "nroot", 1ul, "number of the lowest-lying eigenvectors of the hamiltonian to target"),
        m_spin_restrict(this, "spin_restrict", 0ul,
                        "2Ms value in which to restrict the fermion sector if the Hamiltonian conserves secondary spin quantum number"),
        m_buffers(this), m_hash_mapping(this), m_archivable(this), m_load_balancing(this) {}

fciqmc_config::Reweight::Reweight(config::Group *parent) :
        config::Section(parent, "reweight", "options relating to the on-the-fly correction of population control bias"),
        m_ncycle(this, "ncycle", 1500ul, "number of past cycles to accumulate the reweighting product over"),
        m_delay(this, "delay", 1000ul,
                "number of MC cycles to wait before beginning to reweight coefficients based on historical shift data") {}

fciqmc_config::Shift::Shift(config::Group *parent) :
        config::Section(parent, "shift",
                        "options relating to the diagonal shift parameter and the manner in which it is varied"),
        m_init(this, "init", 0.0, "initial shift relative to the energy of the initial reference MBF"),
        m_damp(this, "damp", 0.05, "damping factor in shift update"),
        m_period(this, "period", 5, "number of MC cycles between shift updates"),
        m_ncycle_av(this, "ncycle_av", 100ul, "number of cycles over which to maintain a rolling average"),
        m_jump(this, "jump", false,
               "ignore growth data in the shift update, and instead use a projected energy estimator"),
        m_reweight(this) {}

fciqmc_config::Semistochastic::Semistochastic(config::Group *parent) :
        config::Section(parent, "semistochastic", "options related to semi-stochastic propagation"),
        m_size(this, "size", 0ul, "number of MBFs selected to comprise the semi-stochastic space"),
        m_l1_fraction_cutoff(this, "l1_fraction_cutoff", 1.0,
                             "requisite fraction of the total number of walkers required to reside on an MBF for inclusion in the semistochastic space"),
        m_delay(this, "delay", 0ul,
                "number of MC cycles to wait after the onset of variable shift mode before initializing the semi-stochastic space") {}

void fciqmc_config::Semistochastic::verify() {
    REQUIRE_NE_ALL(bool(m_size), m_l1_fraction_cutoff == 1.0, "incompatible methods of subspace selection specified");
    REQUIRE_LE_ALL(m_l1_fraction_cutoff, 1.0, "cutoff cannot exceed 1.0");
}

fciqmc_config::Fcidump::Fcidump(config::Group *parent) :
        config::Section(parent, "fcidump", "options relating to the FCIDUMP file"),
        m_path(this, "path", "FCIDUMP", "path to file defining fermionic Hamiltonian"),
        m_eb_path(this, "eb_path", "EBDUMP", "path to file defining fermion-boson coupling"),
        m_bos_path(this, "bos_path", "BOSDUMP", "path to file defining bosonic Hamiltonian"),
        m_spin_major(this, "spin_major", false,
                     "if true, spin-resolved FCIDUMP orders the spin orbitals aaa...bbb..., and ababab... if false.") {}

fciqmc_config::Stats::Stats(config::Group *parent) :
        config::Section(parent, "stats",
                        "options relating to the recording of time series statistics about the calculation"),
        m_path(this, "path", "M7.stats",
               "path to the file to which MC cycle statistics will be output"),
        m_parallel(this, "parallel", false,
                   "output additional stats from each MPI rank") {}

fciqmc_config::SpfWeightedTwf::SpfWeightedTwf(config::Group *parent) :
        config::Section(parent, "spf_weighted_twf",
                        "options related to the weighted trial wavefunction defined for Hamiltonians in which the sign structure of the FCI wavefunction is known"),
        m_fermion_fac(this, "fermion_fac", 0.0,
                      "Exponential constant penalising fermion double-occupancy in weighted TWF"),
        m_boson_fac(this, "boson_fac", 0.0, "Exponential constant penalising boson occupancy in weighted TWF") {}

void fciqmc_config::SpfWeightedTwf::verify() {
    Section::verify();
    if (!defs::enable_bosons) {
        REQUIRE_EQ_ALL(m_boson_fac, 0.0,
                       "Boson exponential parameter is non-zero but bosons are compile time disabled. "
                       "Specify -DENABLE_BOSONS to cmake and recompile");
    }
}

fciqmc_config::Bilinears::Bilinears(config::Group *parent, std::string name, std::string description) :
        config::Section(parent, name, description),
        m_ranks(this, "ranks", {}, "Ranks to accumulate"),
        m_buffers(this), m_hash_mapping(this), m_load_balancing(this), m_archivable(this) {}

void fciqmc_config::Bilinears::verify() {
    for (const auto &rank: m_ranks.get()) {
        REQUIRE_TRUE_ALL(rank.size() == 1 || rank.size() == 4, "invalid rank specifier");
    }
}


fciqmc_config::Rdms::Rdms(config::Group *parent, std::string name, std::string description) :
        Bilinears(parent, name, description),
        m_explicit_ref_conns(this, "explicit_ref_conns", true,
                             "if true, take contributions from reference connections into account exactly") {}

fciqmc_config::SpecMoms::SpecMoms(config::Group *parent, std::string name, std::string description) :
        Bilinears(parent, name, description),
        m_stochastic(this, "stochastic", true,
                     "if false, perform exact evaluation of contributing connections"),
        m_nattempt_per_walker(this, "nattempt_per_walker", 1.0,
                              "number of attempts to generate contributions per integerized walker on the source MBF") {}

fciqmc_config::InstEsts::InstEsts(config::Group *parent) :
        config::Section(parent, "inst_ests",
                        "options relating to instantaneous (MC cycle-resolved) estimators"),

        m_spf_uniform_twf(this, "spf_uniform_twf", false,
                          "switch on estimation of energy by uniform TWF (applicable only in sign problem-free systems)"),
        m_spf_weighted_twf(this) {}

fciqmc_config::RefExcits::RefExcits(config::Group *parent) :
        config::Section(parent, "ref_excits",
                        "options relating to averaged amplitudes of reference-connected MBFs"),
        m_max_exlvl(this, "max_exlvl", 0ul,
                    "maximum excitation level from the reference for which to accumulate average amplitudes"),
        m_buffers(this), m_archivable(this) {}

fciqmc_config::AvEsts::AvEsts(config::Group *parent) :
        config::Section(parent, "av_ests",
                        "options related to quantities estimated from the many-body wavefunction(s) and averaged on-the-fly over a number of MC cycles"),
        m_delay(this, "delay", 1000ul,
                "number of MC cycles to wait after the onset of variable shift mode before beginning to accumulate MEVs"),
        m_ncycle(this, "ncycle", ~0ul,
                 "number of MC cycles for which to accumulate MEVs before terminating the calculation"),
        m_stats_period(this, "stats_period", 100ul,
                       "number of MC cycles between computation and output of all contracted values computed from the averaged estimators"),
        m_stats_path(this, "stats_path", "M7.av_ests.stats", "output path for contracted value statistics"),
        m_rdm(this, "rdm", "options relating to the accumulation and sampling of RDM elements"),
        m_spec_mom(this, "spec_mom", "options relating to the accumulation and sampling of spectral moments"),
        m_ref_excits(this) {}

fciqmc_config::Hamiltonian::Hamiltonian(config::Group *parent) :
        config::Section(parent, "hamiltonian", "options relating to the Hamiltonian operator"),
        m_fcidump(this),
        m_charge(this, "charge", 0,
                          "electron deficit relative to the number given in the FCIDUMP header (positive value to remove elecs)"),
        m_nboson_max(this, "nboson_max", 0ul, "maximum allowed occupation of bosonic modes") {}

void fciqmc_config::Hamiltonian::verify() {
    Section::verify();
    if (!defs::enable_bosons) {
        REQUIRE_EQ_ALL(m_nboson_max, 0ul,
                       "Maximum boson number per mode is non-zero but bosons are compile time disabled. "
                       "Specify -DENABLE_BOSONS to cmake and recompile");
    }

}

fciqmc_config::Propagator::Propagator(config::Group *parent) :
        config::Section(parent, "propagator",
                        "options relating to the propagation of the wavefunction from one MC cycle to the next"),
        m_ncycle(this, "ncycle", ~0ul, "number of MC cycles for which to iterate the solver"),
        m_stochastic(this, "stochastic", true,
                     "if false, perform exact projective FCI (only practical for debugging in small systems)"),
        m_excit_gen(this, "excit_gen", "pchb", "excitation generation algorithm to use"),
        m_nw_target(this, "nw_target", 0ul, "the L1 norm of the wavefunction at which the shift should begin to vary"),
        m_max_bloom(this, "max_bloom", 0.0,
                    "the maximum acceptable magnitude for an off-diagonal propagated contribution. If tau is dynamic, it is updated to keep spawned contributions below this magnitude"),
        m_nadd(this, "nadd", 3.0, "MBFs with weight above this value are granted initiator status"),
        m_tau_init(this, "tau_init", 1e-3, "initial value for the timestep"),
        m_tau_min(this, "tau_min", 1e-5, "minimum allowed value for the dynamic timestep"),
        m_tau_max(this, "tau_max", 1e-2, "maximum allowed value for the dynamic timestep"),
        m_static_tau(this, "static_tau", true, "keep tau value fixed"),
        m_static_probs(this, "static_probs", true, "keep excitation level probabilities fixed"),
        m_min_spawn_mag(this, "min_spawn_mag", 0.4,
                        "spawn magnitude threshold - smaller spawns are stochastically rounded about this value"),
        m_min_death_mag(this, "min_death_mag", 0.0, "death magnitude threshold"),
        m_min_exlvl_prob(this, "min_exlvl_prob", 1e-3,
                         "prevent the probability of drawing an excitation class falling below this threshold"),
        m_exlvl_probs_init(this, "exlvl_probs_init", {},
                           "initial probabilities with which to attempt to draw excitations"),
        m_ndraw_min_for_dynamic(this, "ndraw_min_for_dynamic", 1000ul,
                                "number of spawns logged for excitation type magnitudes to be used in tau and probability update"),
        m_period(this, "period", 10ul,
                 "number of MC cycles between updates of tau and probabilities if requested"),
        m_semistochastic(this) {}

void fciqmc_config::Propagator::verify() {
    if (consts::float_is_zero(m_min_death_mag.get())) {
        m_min_death_mag = m_min_spawn_mag.get();
        log::warn("{} was zero, defaulting to the specified value of {}",
                  m_min_death_mag.m_path.to_string(), m_min_spawn_mag.m_path.to_string());
    }
    if (consts::float_is_zero(m_max_bloom.get())) {
        m_max_bloom = m_nadd.get();
        log::warn("{} was zero, defaulting to the specified value of {}",
                  m_max_bloom.m_path.to_string(), m_nadd.m_path.to_string());
    }
}

fciqmc_config::Document::Document(const yaml::File *file) :
        config::Document(file, "FCIQMC options",
                         "Configuration document prescribing the behavior of an FCIQMC calculation in M7"),
        m_prng(this), m_archive(this), m_wavefunction(this), m_reference(this),
        m_shift(this), m_propagator(this),
        m_hamiltonian(this), m_stats(this), m_inst_ests(this), m_av_ests(this) {}

void fciqmc_config::Document::verify() {
    config::Document::verify();
    REQUIRE_LT_ALL(m_wavefunction.m_nw_init, m_propagator.m_nw_target,
                   "initial number of walkers must not exceed the target population");
    if (m_wavefunction.m_nw_init < m_propagator.m_nadd) {
        m_wavefunction.m_nw_init = m_propagator.m_nadd.get();
        log::warn("initial number of walkers must be at least the initiator threshold");
    }
}