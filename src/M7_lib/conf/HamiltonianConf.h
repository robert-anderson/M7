//
// Created by Robert J. Anderson on 12/8/21.
//

#ifndef M7_HAMILTONIAN_CONF_H
#define M7_HAMILTONIAN_CONF_H

#include "ConfComponents.h"

namespace conf {

    using namespace conf_components;

    struct Fcidump : Section {
        Param<str_t> m_path;
        SingleChoice<str_t> m_unrestrict_style;

        explicit Fcidump(Group *parent);
    };

    struct Bosdump : Section {
        Param<str_t> m_path;

        explicit Bosdump(Group *parent);
    };

    struct Ebdump : Section {
        Param<str_t> m_path;
        Param<bool> m_spin_major;

        explicit Ebdump(Group *parent);
    };

    struct LatticeModel : Section {
        Param<str_t> m_topology;
        Param<uintv_t> m_site_shape;
        MultiChoice<int> m_boundary_conds;

        LatticeModel(Group *parent, str_t name, str_t description);

    protected:
        void validate_node_contents() override;
    };

    struct Hubbard : LatticeModel {
        Param<ham_comp_t> m_repulsion;

        explicit Hubbard(Group *parent);
    };


    struct Heisenberg : LatticeModel {
        Param<ham_comp_t> m_coupling;

        explicit Heisenberg(Group *parent);
    };


    struct FrmHam : Section {
        Fcidump m_fcidump;
        Hubbard m_hubbard;
        Heisenberg m_heisenberg;
        Param<double> m_spin_penalty_j;

        explicit FrmHam(Group *parent);

    protected:
        void validate_node_contents() override;
    };

    struct FrmBosHam : Section {
        Ebdump m_ebdump;
        Param<ham_t> m_holstein_coupling;

        explicit FrmBosHam(Group *parent);
    };

    struct InteractingBoseGas : Section {
        Param<uint_t> m_ndim;
        Param<uint_t> m_nwave;
        Param<ham_t> m_ek_scale;

        explicit InteractingBoseGas(Group *parent);
    };

    struct BosHam : Section {
        Bosdump m_bosdump;
        Param<ham_comp_t> m_num_op_weight;
        InteractingBoseGas m_interacting_bose_gas;
        Hubbard m_hubbard;

        explicit BosHam(Group *parent);
    };

    struct Hamiltonian : Section {
        FrmHam m_fermion;
        FrmBosHam m_ladder;
        BosHam m_boson;

        explicit Hamiltonian(Group *parent);
    };

}

#endif //M7_HAMILTONIAN_CONF_H