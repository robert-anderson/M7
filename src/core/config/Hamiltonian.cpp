//
// Created by anderson on 12/8/21.
//

#include "Hamiltonian.h"

fciqmc_config::Fcidump::Fcidump(config::Group *parent) :
        config::Section(parent, "fcidump", "options relating to the FCIDUMP file"),
        m_path(this, "path", defs::enable_fermions ? "FCIDUMP" : "", "path to file defining fermionic Hamiltonian"),
        m_spin_major(this, "spin_major", false,
                     "if true, spin-resolved FCIDUMP orders the spin orbitals aaa...bbb..., and ababab... if false.") {}

fciqmc_config::Bosdump::Bosdump(config::Group *parent) :
        config::Section(parent, "bosdump",
                        "options relating to 4-indexed text file defining arbitrary number-conserving boson interactions"),
        m_path(this, "path", defs::enable_bosons ? "BOSDUMP" : "", "path to BOSDUMP format file") {}

fciqmc_config::Ebdump::Ebdump(config::Group *parent) :
        config::Section(parent, "bosdump",
                        "options relating to 3-indexed text file defining arbitrary couplings between electron hopping (and one-electron density) and boson (de-)excitations"),
        m_path(this, "path", defs::enable_bosons ? "EBDUMP" : "", "path to EBDUMP format file") {}

fciqmc_config::Hubbard::Hubbard(config::Group *parent) :
        config::Section(parent, "hubbard",
                        "parameters of the arbitrarily dimensioned Hubbard model in the site basis. half-filling is the default, and doping can be achieved by modifying the charge parameter of the hamiltonian section"),
        m_site_shape(this, "site_shape", {},
                     "dimensionality of the orthogonally-coordinated N-dimensional Hubbard lattice"),
        m_boundary_conds(this, "boundary_conds", {},
                         "boundary conditions for each dimension of the Hubbard lattice (-1: anti-periodic, 0: open, 1: periodic)"),
        m_repulsion(this, "repulsion", 0ul, "on-site repulsion coefficient \"U\"") {}

void fciqmc_config::Hubbard::verify() {
    REQUIRE_EQ(m_site_shape.get().size(), m_boundary_conds.get().size(),
               "boundary conditions must be defined for each element of the lattice shape");
}

fciqmc_config::FermionHamiltonian::FermionHamiltonian(config::Group *parent) :
        config::Section(parent, "fermion", "options relating to the fermion hamiltonian terms"),
        m_fcidump(this), m_hubbard(this),
        m_charge(this, "charge", 0,
                 "electron deficit relative to the default value in the FCIDUMP file or that assumed by the model system (positive value to remove elecs)"),
        m_elecs(this, "elecs", true, "include fermionic operators in the Hamiltonian"),
        m_ms2_restrict(this, "ms2_restrict", 0ul,
                       "2Ms value in which to restrict the fermion sector if the Hamiltonian conserves z-axis projection of spin quantum number") {}

fciqmc_config::LadderHamiltonian::LadderHamiltonian(config::Group *parent) :
        config::Section(parent, "ladder",
                        "options relating to the \"ladder\" hamiltonian term, which couples number-conserving electronic single excitations to boson (de-)excitations"),
        m_ebdump(this),
        m_holstein_coupling(this, "holstein_coupling", 0.0,
                            "constant coupling between the electronic density and the boson (de-)excitation operators"),
        m_nboson_max(this, "nboson_max", 0ul,
                     "maximum allowed occupation of bosonic modes. Disregard bosonic operators entirely in the Hamiltonian if set to 0") {}

void fciqmc_config::LadderHamiltonian::verify() {
    if (!defs::enable_bosons) {
        REQUIRE_EQ_ALL(m_nboson_max, 0ul,
                       "Maximum boson number per mode is non-zero but bosons are compile time disabled. "
                       "Set CMake variable -DMBF_TYPE to \"fermion-boson\" (or 1) and recompile");
    }
    REQUIRE_LE_ALL(m_nboson_max, defs::max_bos_occ,
                   log::format("Maximum boson number mustn't exceed the capacity of the integer container ({})",
                               defs::max_bos_occ));
}

fciqmc_config::BosonHamiltonian::BosonHamiltonian(config::Group *parent) :
        config::Section(parent, "boson", "options relating to the number-conserving boson hamiltonian terms"),
        m_bosdump(this),
        m_holstein_omega(this, "holstein_omega", 0.0, "constant frequency of the boson modes in the Holstein model") {}

fciqmc_config::Hamiltonian::Hamiltonian(config::Group *parent) :
        config::Section(parent, "hamiltonian", "options relating to the Hamiltonian operator terms"),
        m_fermion(this), m_ladder(this), m_boson(this) {}
