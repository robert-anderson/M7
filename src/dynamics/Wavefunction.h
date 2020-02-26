//
// Created by Robert John Anderson on 2020-02-04.
//

#ifndef M7_WAVEFUNCTION_H
#define M7_WAVEFUNCTION_H

#include <src/enumerators/VectorCombinationEnumerator.h>
#include "src/dynamics/WalkerList.h"
#include "src/sample/HeatBathSampler.h"
#include "SpawnList.h"

class Wavefunction {
    WalkerList m_walker_list;
    SpawnList m_spawn_list;

public:

    defs::ham_comp_t m_square_norm;
    defs::ham_comp_t m_delta_square_norm;

    Wavefunction(const Determinant &ref, size_t nwalker_initial, size_t nrow_walkers, size_t nrow_send,
                 size_t nrow_recv) :
        m_walker_list(ref.nspatorb() * 2, nrow_walkers),
        m_spawn_list(ref.nspatorb() * 2, nrow_send, nrow_recv) {
        auto irow = m_walker_list.m_list->push(ref);
        *m_walker_list.m_list->view<defs::ham_t>(0) = nwalker_initial;
        m_square_norm = std::pow(nwalker_initial, 2);
    }

    /*
    void evolve(const HeatBathSampler &hb){
#pragma omp parallel for
        for (size_t irow=0ul; irow<m_walker_list.high_water_mark(); ++irow){
            if (m_walker_list.is_free(irow)) continue;
            auto weight = m_walker_list.weight(irow);
            auto nattempt = std::ceil(std::abs(weight));
            auto det = m_walker_list.det(irow);
            {
                auto sampler = hb.sample_excitations(det);
                for (size_t iattempt=0ul; iattempt<nattempt; ++iattempt){
                    auto excit = sampler.draw();
                }
            }
        }
    }*/


    void annihilation() {
        for (size_t irow_recv = 0ul; irow_recv < m_spawn_list.nrecv(); ++irow_recv) {
            auto det = m_spawn_list.recv_det(irow_recv);
            auto spawned_weight = m_spawn_list.recv_weight(irow_recv);

            defs::ham_t stored_weight;
            {
                auto mutex = m_walker_list.m_list->key_mutex(det);
                auto irow_store = m_walker_list.m_list->push(mutex, det);
                stored_weight = *m_walker_list.m_list->view<defs::ham_t>(irow_store, m_walker_list.m_iweight);
                *m_walker_list.m_list->view<defs::ham_t>(irow_store, m_walker_list.m_iweight) += spawned_weight;
            }
            m_delta_square_norm += std::pow(abs(stored_weight + spawned_weight), 2);
            m_delta_square_norm -= std::pow(abs(stored_weight), 2);
        }
        m_spawn_list.m_conn->m_recv.zero();
    }

    void death(const AbInitioHamiltonian &h, const defs::ham_comp_t &shift) {
        for (size_t irow = 0ul; irow < m_walker_list.high_water_mark(); ++irow) {
            auto det = m_walker_list.det(irow);
            auto hdiag = h.get_element_0(det);
            auto stored_weight = m_walker_list.weight(irow);
            auto delta_weight = -0.01 * (hdiag - shift) * stored_weight;
            m_delta_square_norm += std::pow(abs(stored_weight + delta_weight), 2);
            m_delta_square_norm -= std::pow(abs(stored_weight), 2);
            *m_walker_list.m_list->view<defs::ham_t>(irow, m_walker_list.m_iweight) += delta_weight;
        }
    }

    void spawn(const AbInitioHamiltonian &h, const defs::ham_t &weight,
               const Determinant &det, const Determinant &excited) {

        defs::ham_t helement = h.get_element(det, excited);
        if (consts::float_is_zero(helement)) return;
        auto idst_rank = 0;//m_det_block_destinations[DeterminantHasher()(ket) % m_det_block_destinations.size()];
        m_spawn_list.create(excited, idst_rank, -0.01 * helement * weight, false);
        //if (det==m_reference) m_hf_connection_buffer_flag->set(send, idst, irow);
    }

    void evolve_exact(const AbInitioHamiltonian &h) {
        for (size_t irow = 0ul; irow < m_walker_list.high_water_mark(); ++irow) {
            if (m_walker_list.is_free(irow)) continue;
            auto weight = m_walker_list.weight(irow);
            auto nattempt = std::ceil(std::abs(weight));
            auto det = m_walker_list.det(irow);

            auto occs = DeterminantSetEnumerator(det).enumerate();
            assert(occs.size());
            auto unoccs = DeterminantClrEnumerator(det).enumerate();
            assert(unoccs.size());

            Determinant excited(det.nspatorb());
            for (auto occ : occs) {
                for (auto unocc :unoccs) {
                    excited = det.get_excited_det(occ, unocc);
                    spawn(h, weight, det, excited);
                }
            }

            VectorCombinationEnumerator occ_enumerator(occs, 2);
            defs::inds occ_inds(2);

            while (occ_enumerator.next(occ_inds)) {
                {
                    VectorCombinationEnumerator unocc_enumerator(unoccs, 2);
                    defs::inds unocc_inds(2);
                    while (unocc_enumerator.next(unocc_inds)) {
                        excited = det.get_excited_det(occ_inds, unocc_inds);
                        spawn(h, weight, det, excited);
                    }
                }
            }
        }
    }
};

/*
#include "DataSystem.h"
#include "src/data/PerforableMappedTable.h"
#include <memory>
#include <src/propagator/Propagator.h>
#include <src/enumerators/BitfieldEnumerator.h>
#include <src/enumerators/VectorCombinationEnumerator.h>
#include "src/consts.h"
#include "src/io/InputOptions.h"


class Flag {
    const size_t m_ispec;
    const size_t m_iflag;
public:
    Flag(const size_t ispec, const size_t iflag):
    m_ispec(ispec), m_iflag(iflag){}
    bool get(Table& table, const size_t &isegment, const size_t irow){
        return table.view<BitfieldNew>(isegment, irow).get(m_iflag);
    }
    void set(Table& table, const size_t &isegment, const size_t irow){
        return table.view<BitfieldNew>(isegment, irow).set(m_iflag);
    }
    void clr(Table& table, const size_t &isegment, const size_t irow){
        return table.view<BitfieldNew>(isegment, irow).clr(m_iflag);
    }
};

class FlagSet {
    Specification &m_spec;
    const size_t m_ispec;
public:
    FlagSet(Specification &spec):
     m_spec(spec), m_ispec(m_spec.add<BitfieldNew>(0)){}

    auto add(){
        return std::make_unique<Flag>(m_ispec, m_spec.m_bitfield_lengths[m_ispec]++);
    }
};
*/

/*
class Wavefunction {
    const InputOptions &m_input;
    std::unique_ptr<PerforableMappedTable<Determinant>> m_store;
    std::unique_ptr<Table> m_send_buffer;
    std::unique_ptr<Table> m_recv_buffer;

    Specification m_store_spec;
    Specification m_buffer_spec;

    size_t m_weight_store_entry;
    size_t m_hdiag_store_entry;
    size_t m_onv_store_entry;
    size_t m_weight_buffer_entry;
    size_t m_onv_buffer_entry;

    std::unique_ptr<Flag> m_hf_connection_store_flag;
    std::unique_ptr<Flag> m_hf_connection_buffer_flag;



    std::vector<size_t> m_det_block_destinations;

    defs::ham_comp_t m_square_norm;
    defs::ham_comp_t m_delta_square_norm;

    defs::ham_t m_component_norm;
    defs::ham_t m_delta_component_norm;

public:
    Determinant m_reference;

    auto get_hdiag(const size_t &irow) {
        return m_store->view<consts::component_t<defs::ham_t>::type>(0, irow, m_hdiag_store_entry);
    }

    auto get_stored_weight(const size_t &irow) {
        return m_store->view<defs::ham_t>(0, irow, m_weight_store_entry);
    }

    Wavefunction(const InputOptions &input, const Determinant reference) :
        m_input(input), m_reference(reference),
        m_det_block_destinations(input.nload_balance_block * mpi::nrank()) {


        m_weight_store_entry = m_store_spec.add<defs::ham_t>(1);
        m_hdiag_store_entry = m_store_spec.add<consts::component_t<defs::ham_t>::type>(1);
        m_onv_store_entry = m_store_spec.add<Determinant>(reference.nspatorb());

        m_weight_buffer_entry = m_buffer_spec.add<defs::ham_t>(1);
        m_onv_buffer_entry = m_buffer_spec.add<Determinant>(reference.nspatorb());


        FlagSet store_flags(m_store_spec);
        FlagSet buffer_flags(m_store_spec);
        m_hf_connection_store_flag = store_flags.add();
        m_hf_connection_buffer_flag = buffer_flags.add();


        m_store = std::make_unique<PerforableMappedTable<Determinant>>(
            m_store_spec,
            (size_t) (input.nwalker_target * input.store_factor_initial)
        );
        m_send_buffer = std::make_unique<Table>(
            m_buffer_spec, (size_t) (input.nwalker_target * input.store_factor_initial)
        );
        m_recv_buffer = std::make_unique<Table>(
            m_buffer_spec, (size_t) (input.nwalker_target * input.store_factor_initial * mpi::nrank())
        );

        //m_store->push_view<defs::ham_t>(0, reference)[m_weight_store_entry] = input.nwalker_initial;
        m_square_norm = std::pow(input.nwalker_initial, 2);
        m_component_norm = input.nwalker_initial;

         // load balancing
        for (size_t i =0ul; i < m_det_block_destinations.size(); ++i) {
            m_det_block_destinations[i] = i % mpi::nrank();
        }
    }

    void fill_send_buffer(Table &send, const Propagator &prop,
                          const defs::ham_t &weight, Determinant &bra, const Determinant &ket) {
        defs::ham_t helement = prop.m_h.get_element(bra, ket);
        if (consts::float_is_zero(helement)) return;
        auto idst = m_det_block_destinations[DeterminantHasher()(ket) % m_det_block_destinations.size()];
        auto irow = send.safe_push(idst, 1);
        send.view<Determinant>(idst, irow, m_onv_buffer_entry) = ket;
        *send.view<defs::ham_t>(idst, irow, m_weight_buffer_entry) = -prop.tau * helement * weight;
        if (bra==m_reference) m_hf_connection_buffer_flag->set(send, idst, irow);
    }

    void spawn_from_determinant(Table &send, const Propagator &prop, const size_t &irow) {
        auto det = m_store->view<Determinant>(0, irow, m_onv_store_entry);
        auto weight = *m_store->view<defs::ham_t>(0, irow, m_weight_store_entry);
        assert(!det.is_zero());
        auto occs = DeterminantSetEnumerator(det).enumerate();
        assert(occs.size());
        auto unoccs = DeterminantClrEnumerator(det).enumerate();
        assert(unoccs.size());
        assert(unoccs.back()<det.nspatorb()*2);

        Determinant excited(det.nspatorb());
        for (auto occ : occs) {
            for (auto unocc :unoccs) {
                excited = det.get_excited_det(occ, unocc);
                fill_send_buffer(send, prop, weight, det, excited);
            }
        }

        excited = det.get_excited_det(occs[0], occs[1], unoccs[0], unoccs[1]);
        fill_send_buffer(send, prop, weight, det, excited);

        VectorCombinationEnumerator occ_enumerator(occs, 2);
        defs::inds occ_inds(2);

        while (occ_enumerator.next(occ_inds)) {
            {
                VectorCombinationEnumerator unocc_enumerator(unoccs, 2);
                defs::inds unocc_inds(2);
                while (unocc_enumerator.next(unocc_inds)) {
                    excited = det.get_excited_det(occ_inds, unocc_inds);
                    fill_send_buffer(send, prop, weight, det, excited);
                }
            }
        }
    }

    void merge_recv_with_store(const Propagator &prop);

    void spawn(const Propagator &prop) {
        //#pragma omp parallel
        {
            Table thread_send_buffer(m_send_buffer.get(),
                                     m_input.buffer_temp_factor_initial * m_input.nwalker_target);
//#pragma omp for
            for (size_t irow = 0ul; irow < m_store->highwatermark()[0]; ++irow) {
                spawn_from_determinant(thread_send_buffer, prop, irow);
            }
        }
    }

    void death(const Propagator &prop);

    void evolve(Propagator &prop) {
        if (m_store->highwatermark()[0] == 1) *get_hdiag(0) = prop.m_h.get_energy(m_reference);
        spawn(prop);
        m_send_buffer->send_to(*m_recv_buffer);
        m_delta_square_norm = 0;
        m_delta_component_norm = 0;
        death(prop);
        merge_recv_with_store(prop);
        assert(m_store->highwatermark()[0] <
               integer_utils::combinatorial(prop.m_h.nspatorb() * 2, prop.m_h.nelec()));

        if (std::sqrt(m_square_norm+m_delta_square_norm)>m_input.nwalker_target){
            std::cout << "dynamically updating shift" << std::endl;
            prop.update_shift(std::sqrt(m_square_norm+m_delta_square_norm)/std::sqrt(m_square_norm));
            //prop.update_shift(1+m_delta_square_norm/m_square_norm);
        }
        m_square_norm+=m_delta_square_norm;
        m_component_norm+=m_delta_component_norm;
        std::cout << "norm: " << std::sqrt(m_square_norm) <<std::endl;
    }

};*/

#endif //M7_WAVEFUNCTION_H
