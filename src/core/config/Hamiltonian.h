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

        explicit Ebdump(config::Group *parent);

        bool enabled() const override;
    };

    struct Hubbard : config::Section {
        config::Param<defs::inds> m_site_shape;
        config::Param<std::vector<int>> m_boundary_conds;
        config::Param<defs::ham_comp_t> m_repulsion;

        Hubbard(config::Group *parent);

        void verify() override;

        bool enabled() const override;
    };

    struct FermionHamiltonian : config::Section {
        Fcidump m_fcidump;
        Hubbard m_hubbard;
        config::Param<int> m_charge;
        config::Param<int> m_ms2_restrict;

        FermionHamiltonian(config::Group *parent);

        void verify() override;

        bool enabled() const override;
    };

    struct LadderHamiltonian : config::Section {
        Ebdump m_ebdump;
        config::Param<defs::ham_t> m_holstein_coupling;
        config::Param<size_t> m_nboson_max;

        LadderHamiltonian(config::Group *parent);

        void verify() override;

        bool enabled() const override;
    };

    struct BosonHamiltonian : config::Section {
        Bosdump m_bosdump;
        config::Param<defs::ham_comp_t > m_holstein_omega;
        BosonHamiltonian(config::Group *parent);

        bool enabled() const override;
    };

    struct Hamiltonian : config::Section {
        FermionHamiltonian m_fermion;
        LadderHamiltonian m_ladder;
        BosonHamiltonian m_boson;

        Hamiltonian(config::Group *parent);
    };

}


#endif //M7_CONFIG_HAMILTONIAN_H
