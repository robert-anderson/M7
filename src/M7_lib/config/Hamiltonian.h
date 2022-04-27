//
// Created by anderson on 12/8/21.
//

#ifndef M7_CONFIG_HAMILTONIAN_H
#define M7_CONFIG_HAMILTONIAN_H

#include "Parameters.h"

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
        Param<defs::inds> m_site_shape;
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
        Param<size_t> m_nelec;
        Param<int> m_ms2;
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
        Param<size_t> m_ndim;
        Param<size_t> m_nwave;
        Param<defs::ham_t> m_ek_scale;

        explicit InteractingBoseGas(Group *parent);

        bool enabled() const override;
    };

    struct BosHam : Section {
        Bosdump m_bosdump;
        Param<size_t> m_nboson;
        Param<size_t> m_bos_occ_cutoff;
        Param<defs::ham_comp_t> m_num_op_weight;
        InteractingBoseGas m_interacting_bose_gas;

        explicit BosHam(Group *parent);

        bool enabled() const override;

        void verify() override {
            REQUIRE_LE(m_bos_occ_cutoff, defs::max_bos_occ, log::format("given nboson_max exceeds limit of {}", defs::max_bos_occ));
            REQUIRE_LE(m_nboson, m_bos_occ_cutoff, "number of bosons in a number-conserving system mustn't exceed the maximum occupation cutoff");
        }
    };

    struct Hamiltonian : Section {
        FrmHam m_fermion;
        FrmBosHam m_ladder;
        BosHam m_boson;
        

        explicit Hamiltonian(Group *parent);
    };

}


#endif //M7_CONFIG_HAMILTONIAN_H
