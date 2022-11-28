//
// Created by Robert J. Anderson on 12/8/21.
//

#include "HamiltonianConf.h"

conf::Fcidump::Fcidump(Group *parent) :
        Section(parent, "fcidump", "options relating to the FCIDUMP file", Explicit),
        m_path(this, "path", "FCIDUMP", "path to file defining fermionic Hamiltonian"),
        m_unrestrict_style(this, "unrestrict_style",
           {
               {"minor",  "spin orbitals in the ududud... ordering"},
               {"major",  "spin orbitals in the uuu...ddd... ordering"},
               {"blocks", "spatial orbitals in Molpro-style blocks (uu|uu), (uu|dd), (dd|dd), h_uu, h_dd"}
           }, "minor", "protocol assumed when reading integrals expressed in an unrestricted basis"){}

conf::Bosdump::Bosdump(Group *parent) :
        Section(parent, "bosdump", "options relating to 4-indexed text file defining arbitrary "
                                   "number-conserving boson interactions", Explicit),
        m_path(this, "path", "BOSDUMP", "path to BOSDUMP format file") {}

conf::Ebdump::Ebdump(Group *parent) :
        Section(parent, "ebdump","options relating to 3-indexed text file defining arbitrary couplings "
                                 "between electron hopping (and one-electron density) and boson (de-)excitations", Explicit),
        m_path(this, "path", "EBDUMP", "path to EBDUMP format file"),
        m_spin_major(this, "spin_major", false,
                     "if true, spin-resolved EBDUMP orders the spin orbitals aaa...bbb..., and ababab... otherwise.") {}

conf::LatticeModel::LatticeModel(Group *parent, str_t name, str_t description) :
        Section(parent, name, description, Explicit),
        m_topology(this, "topology", "ortho",
                   "geometric layout of the N-dimensional lattice"),
        m_site_shape(this, "site_shape", {},
                     "dimensionality of the N-dimensional lattice"),
        m_boundary_conds(this, "boundary_conds", {{-1, "anti-periodic"}, {0, "open"}, {1, "periodic"}}, {},
                         "boundary conditions for each dimension of the lattice", true){}

void conf::LatticeModel::validate_node_contents() {
    REQUIRE_EQ(m_site_shape.m_value.size(), m_boundary_conds.m_value.size(),
        "boundary conditions must be defined for each element of the lattice shape");
}

conf::Hubbard::Hubbard(Group *parent) :
        LatticeModel(parent, "hubbard",
                     "parameters of the arbitrarily dimensioned Hubbard model in the site basis. half-filling is the default, and doping can be achieved by modifying the charge parameter of the hamiltonian section"),
        m_repulsion(this, "repulsion", 0.0, "on-site repulsion coefficient \"U\" in units of the hopping") {}

conf::Heisenberg::Heisenberg(Group *parent) :
        LatticeModel(parent, "heisenberg",
                     "parameters of the arbitrarily dimensioned Heisenberg spin model in the site basis."),
        m_coupling(this, "coupling", 1.0,
                   "interaction coefficient \"J\": just scales the energy unless there are other terms in the Hamiltonian") {}


conf::FrmHam::FrmHam(Group *parent) :
        Section(parent, "fermion", "options relating to the fermion hamiltonian terms", Implicit),
        m_fcidump(this), m_hubbard(this), m_heisenberg(this),
        m_spin_penalty_j(this, "spin_penalty_j", 0.0, "scalar multiple of the total spin operator used in the spin penalty fermion Hamiltonian modification"){}

void conf::FrmHam::validate_node_contents() {
    uint_t ndefined = m_fcidump.m_enabled + m_hubbard.m_enabled + m_heisenberg.m_enabled;
    REQUIRE_LE(ndefined, 1ul, "conflicting hamiltonian definitions are defined");
}

conf::FrmBosHam::FrmBosHam(Group *parent) :
        Section(parent, "ladder",
                        "options relating to the fermion-boson product term, which couples number-conserving "
                        "electronic single excitations to boson (de-)excitations", Explicit),
        m_ebdump(this),
        m_holstein_coupling(this, "holstein_coupling", 0.0,
                            "constant coupling between the electronic density and the boson (de-)excitation operators"){}

conf::InteractingBoseGas::InteractingBoseGas(Group *parent) :
        Section(parent, "interacting_bose_gas",
                "options relating to the N-dimensional interacting boson gas", Explicit),
        m_ndim(this, "ndim", 1ul, "number of dimensions"),
        m_nwave(this, "nwave", 0ul, "number of plane waves per dimension (number of degrees of freedom is 2*nwave + 1)"),
        m_ek_scale(this, "ek_scale", 0.0, "scale of the kinetic energy in terms of the scattering interaction strength"){}


conf::BosHam::BosHam(Group *parent) :
        Section(parent, "boson", "options relating to the number-conserving boson hamiltonian terms", Explicit),
        m_bosdump(this),
        m_num_op_weight(this, "num_op_weight", 0.0, "scalar factor of the bosonic number operator"),
        m_interacting_bose_gas(this), m_hubbard(this){}

conf::Hamiltonian::Hamiltonian(Group *parent) :
        Section(parent, "hamiltonian", "options relating to the Hamiltonian operator terms"),
        m_fermion(this), m_ladder(this), m_boson(this) {}

