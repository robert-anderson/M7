//
// Created by rja on 25/06/2021.
//

#include "FciqmcConfig.h"

fciqmc_config::Buffers::Buffers(config::Group *parent) :
        config::Section(parent, "buffers",
                        "options relating to the allocation and reallocation behavior of a Communicator"),
        m_store_fac_init(this, "store_fac_init", 0.5,
                         "the number of rows ultimately required by the store buffer is estimated, then that estimate is multiplied by this parameter to determine the initial row allocation"),
        m_store_exp_fac(this, "store_expand_fac", 2.0,
                        "additional number of rows that should be added to the store buffer's capacity as a fraction of the required number of additional rows"),
        m_comm_fac_init(this, "comm_fac_init", 0.5,
                        "the numbers of rows ultimately required by the communicating pair of buffers (send/recv) are estimated, then that estimate is multiplied by this parameter to determine the initial row allocations"),
        m_comm_exp_fac(this, "comm_expand_fac", 2.0,
                       "additional number of rows that should be added to the communicating buffers' capacities as a fraction of the required number of additional rows") {}

fciqmc_config::Io::Io(config::Group *parent, std::string path_default) :
        config::Section(parent, "io",
                        "options relating to filesystem save and load of structures in an M7 calculation"),
        m_save_path(this, "save_path", path_default,
                    "path to which the HDF5 file containing the structure should be saved"),
        m_load_path(this, "load_path", path_default,
                    "path from which the HDF5 file containing the structure should be loaded") {}

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
        config::Section(parent, "reference", "options relating to the reference ONV"),
        m_init_onv(this, "init_reference_onv", {},
                   "string representations of the ONVs to use as the init references for each root. If replication is used, the init state will be used for both replicas of a root population."),
        m_redef_thresh(this, "redef_thresh", 10.0,
                       "when the highest-weighted non-reference ONV (the candidate) reaches this multiple of the weight on the reference, the candidate will be adopted as the new reference") {}

fciqmc_config::Wavefunction::Wavefunction(config::Group *parent) :
        config::Section(parent, "wavefunction",
                        "options relating to the storage and update of a distributed many-body wavefunction"),
        m_nw_init(this, "nw_init", 1ul, "L1 norm of the init wavefunction"),
        m_nroot(this, "nroot", 1ul, "number of the lowest-lying eigenvectors of the hamiltonian to target"),
        m_replicate(this, "replicate", false, "evolve a statistically-independent replica of each walker population"),
        m_spin_restrict(this, "spin_restrict", 0ul,
                        "2Ms value in which to restrict the fermion sector if the Hamiltonian conserves magnetic spin numbers"),
        m_buffers(this), m_io(this), m_load_balancing(this) {}

fciqmc_config::Reweight::Reweight(config::Group *parent) :
        config::Section(parent, "reweight", "options relating to the on-the-fly correction of population control bias"),
        m_ncycle(this, "ncycle", 1500ul, "number of past cycles to accumulate the reweighting product over"),
        m_delay(this, "delay", 1000ul,
                "number of MC cycles to wait before beginning to reweight coefficients based on historical shift data") {}

fciqmc_config::Shift::Shift(config::Group *parent) :
        config::Section(parent, "shift",
                        "options relating to the diagonal shift parameter and the manner in which it is varied"),
        m_init(this, "init", 0.0, "initial shift relative to the energy of the initial reference ONV"),
        m_damp(this, "damp", 0.05, "damping factor in shift update"),
        m_period(this, "period", 5, "number of MC cycles between shift updates"),
        m_ncycle_av(this, "ncycle_av", 100ul, "number of cycles over which to maintain a rolling average"),
        m_jump(this, "jump", false,
               "ignore growth data in the shift update, and instead use a projected energy estimator"),
        m_reweight(this) {}

fciqmc_config::Semistochastic::Semistochastic(config::Group *parent) :
        config::Section(parent, "semistochastic", "options related to semi-stochastic propagation"),
        m_size(this, "size", 0ul, "number of ONVs selected to comprise the semi-stochastic space"),
        m_delay(this, "delay", 0ul,
                "number of MC cycles to wait after the onset of variable shift mode before initializing the semi-stochastic space") {}

fciqmc_config::Fcidump::Fcidump(config::Group *parent) :
        config::Section(parent, "fcidump", "options relating to the FCIDUMP file"),
        m_path(this, "path", "FCIDUMP", "path to file defining fermionic Hamiltonian"),
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
    if (!defs::enable_bosons)
        REQUIRE_EQ_ALL(m_boson_fac, 0.0,
                       "Boson exponential parameter is non-zero but bosons are compile time disabled");
}

fciqmc_config::FermionRdm::FermionRdm(config::Group *parent) :
        config::Section(parent, "fermion_rdm",
                        "options relating to the accumulation and sampling of fermion RDM elements"),
        m_rank(this, "rank", 0ul, "Rank of fermion RDM to accumulate"),
        m_mixed_estimator(this, "mixed_estimator", false,
                          "replace one instance of the wavefunction in the bilinear RDM definition with an SPF TWF"),
        m_buffers(this), m_load_balancing(this) {}

fciqmc_config::InstEsts::InstEsts(config::Group *parent) :
        config::Section(parent, "inst_ests",
                        "options relating to instantaneous (MC cycle-resolved) estimators"),

        m_spf_uniform_twf(this, "spf_uniform_twf", false,
                          "switch on estimation of energy by uniform TWF (applicable only in sign problem-free systems)"),
        m_spf_weighted_twf(this) {}

fciqmc_config::RefExcits::RefExcits(config::Group *parent) :
        config::Section(parent, "ref_excits",
                        "options relating to averaged amplitudes of reference-connected ONVs"),
        m_max_exlvl(this, "max_exlvl", 0ul,
                    "maximum excitation level from the reference for which to accumulate average amplitudes"),
        m_buffers(this) {}


fciqmc_config::PeriodicOutput::PeriodicOutput(config::Group *parent) :
        config::Section(parent, "periodic_output", "options relating to structures which are output periodically to disk"),
        m_period(this, "period", 0ul, "number of MC cycles between outputs"),
        m_path(this, "path", "", "path to which the structure is output - requires {} token when not clobbering"),
        m_clobber(this, "clobber", false, "overwrite the same file with subsequent data"){}

void fciqmc_config::PeriodicOutput::verify() {
    const auto &str = m_path.get();
    if (str.empty()){
        REQUIRE_EQ_ALL(m_period, 0ul, "non-zero period assigned but no path specified for output");
    }
    else {
        auto token_count = std::count(str.cbegin(), str.cend(), '{');
        REQUIRE_EQ_ALL(token_count, !m_clobber, "if clobbering, no token required, else one token is required");
        if (token_count) {
            auto it_open = std::find(str.cbegin(), str.cend(), '{');
            auto it_close = std::find(str.cbegin(), str.cend(), '}');
            REQUIRE_EQ_ALL(std::distance(it_open, it_close), 1l,
                           "path for periodic output should contain only one {} formatting point");
        }
    }
}

fciqmc_config::AvEsts::AvEsts(config::Group *parent) :
        config::Section(parent, "av_ests",
                        "options related to quantities estimated from the many-body wavefunction(s) and averaged on-the-fly over a number of MC cycles"),
        m_delay(this, "delay", 0ul,
                "number of MC cycles to wait after the onset of variable shift mode before beginning to accumulate MEVs"),
        m_ncycle(this, "ncycle", ~0ul,
                 "number of MC cycles for which to accumulate MEVs before terminating the calculation"),
        m_fermion_rdm(this), m_ref_excits(this), m_io(this, "av_ests.h5"),
        m_periodic_output(this){}

fciqmc_config::Hamiltonian::Hamiltonian(config::Group *parent) :
        config::Section(parent, "hamiltonian", "options relating to the Hamiltonian operator"),
        m_fcidump(this),
        m_boson_frequency(this, "boson_frequency", 0.0,
                          "frequency of onsite boson modes for Hubbard-Holstein model"),
        m_boson_coupling(this, "boson_coupling", 0.0,
                         "coupling of onsite boson modes for Hubbard-Holstein model"),
        m_nboson_max(this, "nboson_max", 0ul, "maximum allowed occupation of bosonic modes") {}

void fciqmc_config::Hamiltonian::verify() {
    Section::verify();
    if (!defs::enable_bosons) {
        REQUIRE_EQ_ALL(m_boson_coupling, 0.0,
                       "Boson coupling parameter is non-zero but bosons are compile time disabled");
        REQUIRE_EQ_ALL(m_boson_frequency, 0.0,
                       "Boson frequency parameter is non-zero but bosons are compile time disabled");
        REQUIRE_EQ_ALL(m_nboson_max, 0ul,
                       "Maximum boson number per mode is non-zero but bosons are compile time disabled");
    }
}

fciqmc_config::Propagator::Propagator(config::Group *parent) :
        config::Section(parent, "propagator",
                        "options relating to the propagation of the wavefunction from one MC cycle to the next"),
        m_ncycle(this, "ncycle", 10000ul, "number of MC cycles for which to iterate the solver"),
        m_exact(this, "exact", false, "perform exact projective FCI (only practical for debugging in small systems)"),
        m_excit_gen(this, "excit_gen", "pchb", "excitation generation algorithm to use"),
        m_nw_target(this, "nw_target", 0ul, "the L1 norm of the wavefunction at which the shift should begin to vary"),
        m_max_bloom(this, "max_bloom", 0.0,
                    "the maximum acceptable magnitude for an off-diagonal propagated contribution. If tau is dynamic, it is updated to keep spawned contributions below this magnitude"),
        m_nadd(this, "nadd", 3.0, "ONVs with weight above this value are granted initiator status"),
        m_tau_init(this, "tau_init", 0.001, "initial value for the timestep"),
        m_static_tau(this, "static_tau", true, "keep tau value fixed"),
        m_min_spawn_mag(this, "min_spawn_mag", 0.4,
                        "spawn magnitude threshold - smaller spawns are stochastically rounded about this value"),
        m_min_death_mag(this, "min_death_mag", 0.0, "death magnitude threshold"),
        m_min_excit_class_prob(this, "min_excit_class_prob", 1e-3,
                               "prevent the probability of drawing an excitatation class falling below this threshold"),
        m_psingle_init(this, "psingle_init", 0.0,
                       "initial probability with which to attempt to draw single excitations"),
        m_nenough_spawns_for_dynamic_tau(this, "nenough_spawns_for_dynamic_tau", 1000ul,
                                         "number of spawns logged for excitation type magnitudes to be used in tau update"),
        m_consolidate_spawns(this, "consolidate_spawns", false,
                             "sort and consolidate received spawns so that there is at most one update to any ONV weight in an annihilation loop"),
        m_semistochastic(this) {}

void fciqmc_config::Propagator::verify() {
    if (consts::float_is_zero(m_min_death_mag.get())) {
        m_min_death_mag = m_min_spawn_mag.get();
        log::warn("{} was zero, defaulting to the specified value of {}",
                  m_min_death_mag.m_path.to_string(), m_min_spawn_mag.m_path.to_string());
    }
}

fciqmc_config::Document::Document(const yaml::File *file) :
        config::Document(file, "FCIQMC options",
                         "Configuration document prescribing the behavior of an FCIQMC calculation in M7"),
        m_prng(this), m_wavefunction(this), m_reference(this), m_shift(this), m_propagator(this),
        m_hamiltonian(this), m_stats(this), m_inst_ests(this), m_av_ests(this) {}

void fciqmc_config::Document::verify() {
    config::Group::verify();
    REQUIRE_LT_ALL(m_wavefunction.m_nw_init, m_propagator.m_nw_target,
                   "initial number of walkers must not exceed the target population");
}
