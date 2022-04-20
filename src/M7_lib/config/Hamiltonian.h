//
// Created by anderson on 12/8/21.
//

#ifndef M7_CONFIG_HAMILTONIAN_H
#define M7_CONFIG_HAMILTONIAN_H

#include "Parameters.h"

namespace fciqmc_config {

    struct Fcidump : config::Section {
        config::Param<std::string> m_path;
        config::Param<bool> m_spin_major;

        explicit Fcidump(config::Group *parent);

        bool enabled() const override;
    };

    struct Bosdump : config::Section {
        config::Param<std::string> m_path;

        explicit Bosdump(config::Group *parent);

        bool enabled() const override;
    };

    struct Ebdump : config::Section {
        config::Param<std::string> m_path;
        config::Param<bool> m_spin_major;

        explicit Ebdump(config::Group *parent);

        bool enabled() const override;
    };

    struct LatticeModel : config::Section {
        config::Param<std::string> m_topology;
        config::Param<defs::inds> m_site_shape;
        config::Param<std::vector<int>> m_boundary_conds;

        LatticeModel(config::Group *parent, std::string name, std::string description);

        void verify() override;

        bool enabled() const override;
    };

    struct Hubbard : LatticeModel {
        config::Param<defs::ham_comp_t> m_repulsion;

        explicit Hubbard(config::Group *parent);
    };


    struct Heisenberg : LatticeModel {
        config::Param<defs::ham_comp_t> m_coupling;

        explicit Heisenberg(config::Group *parent);
    };


    struct FermionHamiltonian : config::Section {
        Fcidump m_fcidump;
        Hubbard m_hubbard;
        Heisenberg m_heisenberg;
        config::Param<int> m_charge;
        config::Param<int> m_ms2_restrict;
        config::Param<double> m_spin_penalty_j;

        explicit FermionHamiltonian(config::Group *parent);

        void verify() override;

        bool enabled() const override;
    };

    struct FrmBosHamiltonian : config::Section {
        Ebdump m_ebdump;
        config::Param<defs::ham_t> m_holstein_coupling;

        explicit FrmBosHamiltonian(config::Group *parent);

        bool enabled() const override;
    };

    struct InteractingBoseGas : config::Section {
        config::Param<size_t> m_ndim;
        config::Param<size_t> m_nwave;
        config::Param<defs::ham_t> m_ek_scale;

        explicit InteractingBoseGas(config::Group *parent);

        bool enabled() const override;
    };

    struct BosonHamiltonian : config::Section {
        Bosdump m_bosdump;
        config::Param<size_t> m_nboson;
        config::Param<size_t> m_nboson_max;
        config::Param<defs::ham_comp_t> m_holstein_omega;
        InteractingBoseGas m_interacting_bose_gas;

        explicit BosonHamiltonian(config::Group *parent);

        bool enabled() const override;

        void verify() override {
            REQUIRE_LE(m_nboson_max, defs::max_bos_occ, log::format("given nboson_max exceeds limit of {}", defs::max_bos_occ));
            REQUIRE_LE(m_nboson, m_nboson_max, "number of bosons in a number-conserving system mustn't exceed the maximum occupation cutoff");
        }
    };

    struct Hamiltonian : config::Section {
        FermionHamiltonian m_fermion;
        FrmBosHamiltonian m_ladder;
        BosonHamiltonian m_boson;

        explicit Hamiltonian(config::Group *parent);
    };

}


#endif //M7_CONFIG_HAMILTONIAN_H
