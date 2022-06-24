//
// Created by Robert J. Anderson on 12/8/21.
//

#ifndef M7_HAMILTONIAN_CONF_H
#define M7_HAMILTONIAN_CONF_H

#include "ConfComponents.h"

namespace conf {

    using namespace conf_components;

    struct Fcidump : Section {
        Param<std::string> m_path;
        Param<bool> m_spin_major;

        explicit Fcidump(Group *parent);

        bool enabled() const override;
    };

    struct Bosdump : Section {
        Param<std::string> m_path;

        explicit Bosdump(Group *parent);

        bool enabled() const override;
    };

    struct Ebdump : Section {
        Param<std::string> m_path;
        Param<bool> m_spin_major;

        explicit Ebdump(Group *parent);

        bool enabled() const override;
    };

    struct LatticeModel : Section {
        Param<std::string> m_topology;
        Param<defs::uintv_t> m_site_shape;
        Param<std::vector<int>> m_boundary_conds;

        LatticeModel(Group *parent, std::string name, std::string description);

        void verify() override;

        bool enabled() const override;
    };

    struct Hubbard : LatticeModel {
        Param<defs::ham_comp_t> m_repulsion;

        explicit Hubbard(Group *parent);
    };


    struct Heisenberg : LatticeModel {
        Param<defs::ham_comp_t> m_coupling;

        explicit Heisenberg(Group *parent);
    };


    struct FrmHam : Section {
        Fcidump m_fcidump;
        Hubbard m_hubbard;
        Heisenberg m_heisenberg;
        Param<double> m_spin_penalty_j;

        explicit FrmHam(Group *parent);

        void verify() override;

        bool enabled() const override;
    };

    struct FrmBosHam : Section {
        Ebdump m_ebdump;
        Param<defs::ham_t> m_holstein_coupling;

        explicit FrmBosHam(Group *parent);

        bool enabled() const override;
    };

    struct InteractingBoseGas : Section {
        Param<uint_t> m_ndim;
        Param<uint_t> m_nwave;
        Param<defs::ham_t> m_ek_scale;

        explicit InteractingBoseGas(Group *parent);

        bool enabled() const override;
    };

    struct BosHam : Section {
        Bosdump m_bosdump;
        Param<defs::ham_comp_t> m_num_op_weight;
        InteractingBoseGas m_interacting_bose_gas;

        explicit BosHam(Group *parent);

        bool enabled() const override;
    };

    struct Hamiltonian : Section {
        FrmHam m_fermion;
        FrmBosHam m_ladder;
        BosHam m_boson;

        explicit Hamiltonian(Group *parent);
    };

}

#endif //M7_HAMILTONIAN_CONF_H