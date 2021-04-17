//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_WAVEFUNCTION_H
#define M7_WAVEFUNCTION_H


#include <src/core/io/Options.h>
#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/parallel/Reduction.h>
#include "src/core/table/Communicator.h"
#include "src/core/dynamics/WalkerTable.h"
#include "src/core/dynamics/SpawnTable.h"
#include "src/core/parallel/RankAllocator.h"
#include "src/core/parallel/ReductionMember.h"
#include "src/core/field/Fields.h"
#include "src/core/sort/QuickSorter.h"

struct Wavefunction : Communicator<WalkerTableRow, SpawnTableRow> {

    const Options &m_opts;
    const size_t m_nsite;

    NdEnumeration<defs::ndim_wf> m_part_inds;

    ReductionSyndicate m_summables;

    NdReduction<size_t, defs::ndim_wf> m_ninitiator;
    NdReduction<int, defs::ndim_wf> m_delta_ninitiator;
    Reduction<size_t> m_nocc_onv;
    Reduction<int> m_delta_nocc_onv;
    NdReduction<defs::wf_comp_t, defs::ndim_wf> m_nwalker;
    NdReduction<defs::wf_comp_t, defs::ndim_wf> m_delta_nwalker;
    NdReduction<defs::wf_comp_t, defs::ndim_wf> m_l2_norm_square;
    NdReduction<defs::wf_comp_t, defs::ndim_wf> m_delta_l2_norm_square;
    NdReduction<defs::wf_comp_t, defs::ndim_wf> m_nannihilated;

    MappedTable<UniqueOnvRow> m_unique_recvd_onvs;
    MappedTable<OnvRow> m_parent_recvd_onvs;

    Wavefunction(const Options &opts, size_t nsite);

    std::vector<std::string> h5_field_names();

    void h5_write(hdf5::GroupWriter &parent, std::string name = "wavefunction");

    void h5_read(hdf5::GroupReader &parent, const Hamiltonian<> &ham, const fields::Onv<> &ref,
                 std::string name = "wavefunction");


    void begin_cycle();

    void end_cycle();

//    void update(size_t icycle, double work_time) {
//        m_ra.update(icycle, work_time, m_walkers, m_walkers.m_key_field);
//    }


    defs::wf_comp_t square_norm(const size_t& ipart) const;

    defs::wf_comp_t l1_norm(const size_t& ipart) const;

    void grant_initiator_status(const size_t &ipart);

    void revoke_initiator_status(const size_t &ipart);

    void set_weight(const size_t &ipart, const defs::wf_t &new_weight);

    void change_weight(const size_t &ipart, const defs::wf_t &delta);

    void scale_weight(const size_t &ipart, const double &factor);

    void zero_weight(const size_t &ipart);

    size_t create_row(const fields::Onv<> &onv, const defs::ham_comp_t &hdiag);

    void remove_row();

    size_t create_walker_(const size_t& icycle, const size_t& ipart, const fields::Onv<> &onv, const defs::ham_t weight,
                          const defs::ham_comp_t &hdiag, bool refconn);


    TableBase::Loc create_walker(const size_t &icycle, const size_t &ipart, const fields::Onv<> &onv,
                                 const defs::ham_t weight, const defs::ham_comp_t &hdiag, bool refconn);

    // create on all parts
    TableBase::Loc create_walker(const size_t &icycle, const fields::Onv<> &onv,
                                 const defs::ham_t weight, const defs::ham_comp_t &hdiag, bool refconn);

    size_t add_spawn(const fields::Onv<> &dst_onv, const defs::wf_t &delta,
                     bool initiator, bool deterministic, size_t dst_ipart);

    size_t add_spawn(const fields::Onv<> &dst_onv, const defs::wf_t &delta,
                     bool initiator, bool deterministic, size_t dst_ipart,
                     const fields::Onv<> &src_onv, const defs::wf_t &src_weight);

    const size_t& npart(){
        return m_part_inds.nelement();
    }

};

#endif //M7_WAVEFUNCTION_H
