//
// Created by rja on 28/07/2021.
//

#include "Hamiltonian.h"
#include "M7_lib/hamiltonian/frm/GeneralFrmHam.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"
#include "M7_lib/hamiltonian/frmbos/InteractingBoseGasBosHam.h"

BasisData Hamiltonian::make_bd() const {
    if (m_frmbos->enabled()) return m_frmbos->m_bd;
    FrmBasisData frm_bd(m_frm->m_nsite, m_frm->m_point_group_map);
    BosBasisData bos_bd(m_bos->m_nmode, m_nboson_max);
    return {frm_bd, bos_bd};
}

std::unique_ptr<FrmHam> Hamiltonian::make_frm(const fciqmc_config::FermionHamiltonian &opts) {
    if (opts.m_hubbard.enabled())
        return make_frm<HubbardFrmHam>(opts);
    else if (opts.m_heisenberg.enabled())
        return make_frm<HeisenbergFrmHam>(opts);
    else if (opts.m_fcidump.enabled())
        return make_frm<GeneralFrmHam>(opts);
    return std::unique_ptr<FrmHam>(new NullFrmHam);
}

std::unique_ptr<FrmBosHam> Hamiltonian::make_ladder(const fciqmc_config::LadderHamiltonian &opts, size_t nsite) {
    if (opts.m_holstein_coupling) {
        auto nboson_max = opts.m_nboson_max.get();
        auto g = opts.m_holstein_coupling.get();
        return std::unique_ptr<FrmBosHam>(new HolsteinLadderHam(nsite, nboson_max, g));
    }
    else if (opts.m_ebdump.enabled()) return std::unique_ptr<FrmBosHam>(new GeneralLadderHam(opts));
    return std::unique_ptr<FrmBosHam>(new NullLadderHam);
}

std::unique_ptr<BosHam> Hamiltonian::make_bos(const fciqmc_config::BosonHamiltonian &opts, size_t nsite) {
    if (opts.m_holstein_omega) {
        auto omega = opts.m_holstein_omega.get();
        return std::unique_ptr<BosHam>(new HolsteinBosHam(nsite, omega));
    }
    else if (opts.m_interacting_bose_gas.enabled())
        return std::unique_ptr<BosHam>(new InteractingBoseGasBosHam(opts));
    else if (opts.m_bosdump.enabled())
        return std::unique_ptr<BosHam>(new GeneralBosHam(opts));
    return std::unique_ptr<BosHam>(new NullBosHam);
}

Hamiltonian::Hamiltonian(const fciqmc_config::Hamiltonian &opts) :
        m_frm(make_frm(opts.m_fermion)),
        m_frmbos(make_ladder(opts.m_ladder, m_frm->m_nsite)),
        m_bos(make_bos(opts.m_boson, m_frm->m_nsite)),
        m_nboson_max(m_frmbos->m_nboson_max),
        m_bd(make_bd()), m_work_conn(m_bd){
    REQUIRE_TRUE(m_bd.m_frm.m_nsite || m_bd.m_bos.m_nmode, "No system defined");
    if (m_frm->disabled()) log::info("Fermion Hamiltonian is disabled");
    if (defs::enable_bosons) {
        if (m_frmbos->disabled()) log::info("Fermion-boson ladder Hamiltonian is disabled");
        if (m_bos->disabled()) log::info("Number-conserving boson Hamiltonian is disabled");
    }
}

size_t Hamiltonian::nci() const {
    return m_frm->nci() * m_bos->nci();
}

size_t Hamiltonian::nelec() const {
    return m_frm->m_nelec;
}

size_t Hamiltonian::nboson() const {
    return m_bos->m_nboson;
}

bool Hamiltonian::complex_valued() const {
    return m_frm->m_complex_valued;
}
