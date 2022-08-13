//
// Created by Robert J. Anderson on 25/06/2021.
//

#include "Conf.h"

using namespace conf_components;

conf::HashMapping::HashMapping(Group *parent) :
        Section(parent, "hash_mapping", "options relating to the behavior of hash-mapped tables"),
        m_remap_ratio(this, "remap_ratio", 2.0,
                      "ratio of bucket-searching skips to total lookups required to trigger remapping with a larger number of buckets"),
        m_remap_nlookup(this, "remap_nlookup", 500ul,
                        "number of lookups required before remapping based on the skips/lookups ratio is considered") {}

conf::Buffers::Buffers(Group *parent) :
        Section(parent, "buffers",
                        "options relating to the allocation and reallocation behavior of a Communicator"),
        m_store_fac_init(this, "store_fac_init", 1.0,
                         "a crude estimate for the number of rows ultimately required by the store buffer is computed, then that estimate is multiplied by this parameter to determine the initial row allocation"),
        m_store_exp_fac(this, "store_expand_fac", 0.5,
                        "additional number of rows that should be added to the store buffer's capacity as a fraction of the required number of additional rows"),
        m_comm_fac_init(this, "comm_fac_init", 1.0,
                        "a crude estimate for the number of rows per rank ultimately required by the send buffer is computed, then that estimate is multiplied by this parameter to determine the initial row allocation"),
        m_comm_exp_fac(this, "comm_expand_fac", 0.5,
                       "additional number of rows that should be added to the communicating buffers' capacities as a fraction of the required number of additional rows") {}


conf::Archive::Archive(Group *parent) :
        Section(parent, "archive",
                        "options relating to archives: HDF5 binary storage of calculation data"),
        m_load_path(this, "load_path", "",
                    "path to the HDF5 file from which the calculation is to be restarted"),
        m_save_path(this, "save_path", "",
                    "path to the HDF5 file into which calculation data is to be dumped at the end of the calculation"),
        m_chkpt_path(this, "chkpt_path", "",
                     "path to the HDF5 file (or file name format) into which calculation data is to be dumped periodically during the calculation"),
        m_period(this, "period", 0ul, "number of MC cycles between checkpoint dumps"),
        m_period_mins(this, "period_mins", 0ul, "time in minutes between checkpoint dumps") {}
        
void conf::Archive::validate_node_contents() {
    if (bool(m_period) && bool(m_period_mins))
        logging::warn("both cycle number and time periods are defined for checkpointing");

    auto &str = m_chkpt_path.m_value;
    uint_t token_count = std::count(str.cbegin(), str.cend(), '{');
    REQUIRE_LE(token_count, 1ul, "checkpoint paths can have at most one {} token");
    if (token_count) {
        auto it_open = std::find(str.cbegin(), str.cend(), '{');
        auto it_close = std::find(str.cbegin(), str.cend(), '}');
        REQUIRE_EQ(std::distance(it_open, it_close), 1l,
                   "checkpoint path for periodic output should contain at most one {} formatting point");
        logging::info("formatting token found in path, "
                  "successive checkpoints will not overwrite previous checkpoints from the same run");
    } else
        logging::info("formatting token not found in path, "
                  "successive checkpoints will overwrite previous checkpoints from the same run");
}

conf::Archivable::Archivable(Group *parent) :
        Section(parent, "archive",
                        "options relating to archiving behavior"),
        m_load(this, "load", false,
               "attempt to load the object from the archive at the beginning of the calculation"),
        m_save(this, "save", false,
               "save the object to the archive at the end of the calculation"),
        m_chkpt(this, "chkpt", false,
                "attempt to save the object to the checkpoint archive at each checkpointing period") {}


conf::Prng::Prng(Group *parent) :
        Section(parent, "prng",
                        "options relating to the random number generator used in stochastic calculations"),
        m_seed(this, "seed", 123ul, "value with which to seed the mt19937 PRNG"),
        m_ngen_block(this, "ngen_block", 10000ul,
                     "size of the block of PRNGs generated each time the buffer is depleted") {}

conf::LoadBalancing::LoadBalancing(Group *parent) :
        Section(parent, "load_balancing",
                        "options relating to the allocation of records among MPI ranks so as to share the workload more equally"),
        m_nblock_per_rank(this, "nblock_per_rank", 6ul, "number of rank allocation blocks to create per MPI rank"),
        m_period(this, "period", 10ul, "number of MC cycles between load-balancing block transactions"),
        m_acceptable_imbalance(this, "acceptable_imbalance", 0.05,
                               "fractional difference in the work figure of the busiest and laziest ranks at which the load balancing is deactivated after a given number of consecutive failed attempts"),
        m_nnull_updates_deactivate(this, "nnull_updates_deactivate", 20ul,
                                   "number of consecutive attempted updates which do not exceed the maximum acceptable imbalance required to meet the deactivation criterion") {}

conf::Reference::Reference(Group *parent) :
        Section(parent, "reference", "options relating to the reference MBF"),
        m_mbf_init(this, "mbf_init"),
        m_redef_thresh(this, "redef_thresh", 2.0,
            "when the highest-weighted non-reference MBF (the candidate) reaches this multiple of the weight on the "
            "reference, the candidate will be adopted as the new reference. set to exactly 0.0 to disable") {}

conf::Wavefunction::Wavefunction(Group *parent) :
        Section(parent, "wavefunction",
                        "options relating to the storage and update of a distributed many-body wavefunction"),
        m_nw_init(this, "nw_init", 1ul, "L1 norm of the initial wavefunction"),
        m_nroot(this, "nroot", 1ul, "number of the lowest-lying eigenvectors of the hamiltonian to target"),
        m_fci_init(this, "fci_init", false, "call the ARPACK interface to initialize the required roots to their exact values"),
        m_buffers(this), m_hash_mapping(this), m_archivable(this), m_load_balancing(this) {}

conf::Shift::Shift(Group *parent) :
        Section(parent, "shift",
                        "options relating to the diagonal shift parameter and the manner in which it is varied"),
        m_init(this, "init", 0.0, "initial shift relative to the energy of the initial reference MBF"),
        m_damp(this, "damp", 0.05, "walker growth-related damping factor in shift update"),
        m_target_damp(this, "target_damp", 0.0, "damping factor in shift update related to growth relative to the target walker population"),
        m_period(this, "period", 5, "number of MC cycles between shift updates"),
        m_ncycle_av(this, "ncycle_av", 100ul, "number of cycles over which to maintain a rolling average"),
        m_jump(this, "jump", false,
               "ignore growth data in the shift update, and instead use a projected energy estimator")
               {}

conf::Semistochastic::Semistochastic(Group *parent) :
        Section(parent, "semistochastic", "options related to semi-stochastic propagation", Explicit),
        m_size(this, "size", 0ul, "number of MBFs selected to comprise the semi-stochastic space"),
        m_l1_fraction_cutoff(this, "l1_fraction_cutoff", 1.0,
                             "requisite fraction of the total number of walkers required to reside on an MBF for inclusion in the semistochastic space"),
        m_delay(this, "delay", 0ul,
                "number of MC cycles to wait after the onset of variable shift mode before initializing the semi-stochastic space") {}

void conf::Semistochastic::validate_node_contents() {
    REQUIRE_NE(bool(m_size), m_l1_fraction_cutoff == 1.0, "incompatible methods of subspace selection specified");
    REQUIRE_LE(m_l1_fraction_cutoff, 1.0, "cutoff must not exceed 1.0");
}

conf::Stats::Stats(Group *parent) :
        Section(parent, "stats",
                        "options relating to the recording of time series statistics about the calculation"),
        m_path(this, "path", "M7.stats",
               "path to the file to which MC cycle statistics will be output"),
        m_period(this, "period", 1,
               "number of MC cycles between output of averaged stats"),
        m_parallel(this, "parallel", false,
                   "output additional stats from each MPI rank") {}

void conf::Stats::validate_node_contents() {
    REQUIRE_TRUE(m_period.m_value, "Stats output period must be non-zero");
}

conf::SpfWeightedTwf::SpfWeightedTwf(Group *parent) :
        Section(parent, "spf_weighted_twf",
                        "options related to the weighted trial wavefunction defined for Hamiltonians in which the "
                        "sign structure of the FCI wavefunction is known", Explicit),
        m_fermion_fac(this, "fermion_fac", 0.0,
                      "Exponential constant penalising fermion double-occupancy in weighted TWF"),
        m_boson_fac(this, "boson_fac", 0.0, "Exponential constant penalising boson occupancy in weighted TWF") {}


void conf::SpfWeightedTwf::validate_node_contents() {
    if (!c_enable_bosons)
        REQUIRE_EQ(m_boson_fac.m_value, 0.0, "Boson exponential parameter is non-zero but bosons are compile time "
                                             "disabled. Set CMake variable -DMBF_TYPE_IND to 1 or 2 and recompile");
}

conf::Bilinears::Bilinears(Group *parent, str_t name, str_t description) :
        Section(parent, name, description, Explicit),
        m_ranks(this, "ranks", {}, "Ranks to accumulate"),
        m_buffers(this), m_hash_mapping(this), m_load_balancing(this), m_archivable(this) {}

void conf::Bilinears::validate_node_contents() {
    for (const auto &rank: m_ranks.m_value)
        REQUIRE_TRUE(rank.size() == 1 || rank.size() == 4, "invalid rank specifier");
}

conf::Rdms::Rdms(Group *parent, str_t name, str_t description) :
        Bilinears(parent, name, description),
        m_explicit_ref_conns(this, "explicit_ref_conns", true,
                             "if true, take contributions from reference connections into account exactly") {}

conf::SpecMoms::SpecMoms(Group *parent, str_t name, str_t description) :
        Bilinears(parent, name, description),
        m_stochastic(this, "stochastic", true,
                     "if false, perform exact evaluation of contributing connections"),
        m_nattempt_per_walker(this, "nattempt_per_walker", 1.0,
                              "number of attempts to generate contributions per integerized walker on the source MBF") {}

conf::InstEsts::InstEsts(Group *parent) :
        Section(parent, "inst_ests",
                        "options relating to instantaneous (MC cycle-resolved) estimators", Explicit),

        m_spf_uniform_twf(this, "spf_uniform_twf", false,
                          "switch on estimation of energy by uniform TWF (applicable only in sign problem-free systems)"),
        m_spf_weighted_twf(this) {}

conf::RefExcits::RefExcits(Group *parent) :
        Section(parent, "ref_excits",
                        "options relating to averaged amplitudes of reference-connected MBFs", Explicit),
        m_max_exlvl(this, "max_exlvl", 0ul,
                    "maximum excitation level from the reference for which to accumulate average amplitudes"),
        m_buffers(this), m_archivable(this) {}

conf::AvEsts::AvEsts(Group *parent) :
        Section(parent, "av_ests",
                        "options related to quantities estimated from the many-body wavefunction(s) and averaged "
                        "on-the-fly over a number of MC cycles", Explicit),
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


conf::GuideWavefunction::ExpFac::ExpFac(conf::GuideWavefunction* parent, str_t name, str_t description) :
        Section(parent, name, description, conf_components::Explicit),
        m_fac(this, "fac", 0.0, "exponential scaling factor"){}

void conf::GuideWavefunction::ExpFac::validate_node_contents() {
    if (m_fac.m_value < 0.0) logging::warn("importance sampling exponent factor is negative");
}

conf::GuideWavefunction::GutzwillerLike::GutzwillerLike(conf::GuideWavefunction* parent) :
        ExpFac(parent, "gutzwiller_like",
               "amplitudes are given by the exponential of the basis function energy"){}

conf::GuideWavefunction::SuppressMultiOcc::SuppressMultiOcc(conf::GuideWavefunction* parent) :
        ExpFac(parent, "suppress_multi_occ",
               "amplitudes are given by the exponential of the number of double occupancies (for fermions) or the "
               "sum of the pair number operator values (for bosons)"){}

conf::GuideWavefunction::GuideWavefunction(Group* parent, const str_t& name) :
    Selection(parent, name, "options relating to a wavefunction used in guiding propagation",
        Exactly, 1ul, conf_components::Explicit),
        m_gutzwiller_like(this), m_suppress_multi_occ(this){}


conf::Propagator::Propagator(Group *parent) :
        Section(parent, "propagator",
                        "options relating to the propagation of the wavefunction from one MC cycle to the next"),
        m_ncycle(this, "ncycle", ~0ul, "number of MC cycles for which to iterate the solver"),
        m_stochastic(this, "stochastic", true,
                     "if false, perform exact projective FCI (only practical for debugging in small systems)"),
        m_excit_gen(this),
        m_nw_target(this, "nw_target", 10000ul, "the L1 norm of the wavefunction at which the shift should begin to vary"),
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
        m_imp_samp_guide(this, "imp_samp_guide"), m_semistochastic(this) {}

void conf::Propagator::validate_node_contents() {
    if (m_min_death_mag.m_value==0.0) {
        m_min_death_mag.m_value = m_min_spawn_mag.m_value;
        logging::warn("{} was zero, defaulting to the specified value of {}",
                  m_min_death_mag.m_path.m_string, m_min_spawn_mag.m_path.m_string);
    }
    if (m_max_bloom.m_value==0.0) {
        m_max_bloom.m_value = m_nadd.m_value;
        logging::warn("{} was zero, defaulting to the specified value of {}",
                  m_max_bloom.m_path.m_string, m_nadd.m_path.m_string);
    }
}

conf::Document::Document(const str_t& fname) :
        conf_components::Document(fname, "a calculation in M7"),
        m_prng(this), m_archive(this), m_basis(this), m_particles(this),
        m_wavefunction(this), m_reference(this),
        m_shift(this), m_propagator(this),
        m_hamiltonian(this), m_stats(this), m_inst_ests(this), m_av_ests(this) {}

void conf::Document::validate_node_contents() {
    REQUIRE_LE(m_wavefunction.m_nw_init, m_propagator.m_nw_target,
               "initial number of walkers must not exceed the target population");
    if (m_wavefunction.m_nw_init < m_propagator.m_nadd) {
        m_wavefunction.m_nw_init.m_value = m_propagator.m_nadd.m_value;
        logging::warn("initial number of walkers must be at least the initiator threshold");
    }
    REQUIRE_LE(m_particles.m_nboson, m_basis.m_bos_occ_cutoff,
               "number of bosons in a number-conserving system mustn't exceed the maximum occupation cutoff");
}

conf::MbfDef::MbfDef(Group *parent, str_t name) :
        Section(parent, name, "definition of a vector of many-body basis functions", Explicit),
        m_frm(this, "fermion", {},
              "fermion sector occupation of the MBF (each element can be a boolean "
              "array of occupation status of all spin-orbitals or list of occupied spin-obrital indices)"),
        m_bos(this, "boson", {}, "boson sector occupation of the MBF as an arrays of occupation "
                                   "levels of all modes"),
        m_neel(this, "neel", false,
                        "initialize the MBF to a Neel state rather than assuming Aufbau principle"){}
