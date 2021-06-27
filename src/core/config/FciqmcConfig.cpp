//
// Created by rja on 25/06/2021.
//

#include "FciqmcConfig.h"

fciqmc_config::Prng::Prng(config::Group *parent) :
        config::Section(parent, "prng",
                        "options relating to the random number generator used in stochastic calculations"),
        m_seed(this, "seed", 123ul, "value with which to seed the PRNG"),
        m_ngen_block(this, "ngen_block", 10000ul,
                     "size of the block of PRNGs generated each time the buffer is depleted") {}

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

fciqmc_config::Serialization::Serialization(config::Group *parent) :
        config::Section(parent, "serialization",
                        "options relating to filesystem save and load of structures in an M7 calculation"),
        m_save_path(this, "save_path", {}, "path to which the HDF5 file containing the structure should be saved"),
        m_load_path(this, "load_path", {}, "path from which the HDF5 file containing the structure should be loaded") {}

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
        m_serialization(this), m_load_balancing(this), m_reference(this) {}

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

fciqmc_config::Propagator::Propagator(config::Group *parent) :
        config::Section(parent, "propagator",
                        "options relating to the propagation of the wavefunction from one MC cycle to the next"),
        m_exact(this, "exact", false, "perform exact projective FCI (only practical for debugging in small systems)"),
        m_nw_target(this, "nw_target", 0ul, "the L1 norm of the wavefunction at which the shift should begin to vary"),
        m_max_bloom(this, "max_bloom", 0.0,
                    "the maximum acceptable magnitude for an off-diagonal propagated contribution. If tau is dynamic, it is updated to keep spawned contributions below this magnitude"),
        m_nadd(this, "nadd", 3.0, "ONVs with weight above this value are granted initiator status"),
        m_tau_init(this, "tau_init", 0.001, "initial value for the timestep"),
        m_static_tau(this, "static_tau", true, "keep tau value fixed"),
        m_min_spawn_mag(this, "min_spawn_mag", 0.4,
                        "spawn magnitude threshold - smaller spawns are stochastically rounded about this value"),
        m_min_death_mag(this, "min_death_mag", 0.0, "death magnitude threshold"),
        m_consolidate_spawns(this, "consolidate_spawns", false,
                             "sort and consolidate received spawns so that there is at most one update to any ONV weight in an annihilation loop") {}

fciqmc_config::Document::Document(const yaml::File *file) :
        config::Document(file, "FCIQMC options",
                         "Configuration document prescribing the behavior of an FCIQMC calculation in M7"),
        m_prng(this), m_wavefunction(this), m_shift(this), m_propagator(this) {}

fciqmc_config::Semistochastic::Semistochastic(config::Group *parent) :
        config::Section(parent, "semistochastic", "options related to semi-stochastic propagation"),
        m_size(this, "size", 0ul, "number of ONVs selected to comprise the semi-stochastic space"){}
