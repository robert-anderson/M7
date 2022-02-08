//
// Created by Robert John Anderson on 2020-01-18.
//

#include <src/core/enumerator/Enumerators.h>
#include "DenseHamiltonian.h"
#include "src/core/util/Foreach.h"

void DenseHamiltonian::setup_frm(const Hamiltonian &source) {

    buffered::FrmOnv bra(source.m_bd);
    buffered::FrmOnv ket(source.m_bd);

    size_t ibra = ~0ul;
    enums::FermionOnv bra_enum(source.m_bd.m_nsite, source.nelec());
    conn::FrmOnv conn(source.m_bd);

    bool spins_only = dynamic_cast<const HeisenbergFrmHam*>(source.m_frm.get());

    while (bra_enum.next(bra, ibra)) {
        if (spins_only && !bra.all_sites_single_occ()) continue;
        size_t iket = ~0ul;
        enums::FermionOnv ket_enum(source.m_bd.m_nsite, source.nelec());
        while (ket_enum.next(ket, iket)) {
            if (spins_only && !ket.all_sites_single_occ()) continue;
            conn.connect(bra, ket);
            auto h_elem = source.get_element(bra, conn);
            (*this)(ibra, iket) = h_elem;
        }
    }
}

void DenseHamiltonian::setup_frmbos(const Hamiltonian &source) {
    buffered::FrmBosOnv bra(source.m_bd);
    buffered::FrmBosOnv ket(source.m_bd);

    size_t ibra = ~0ul;
    enums::FermiBosOnv bra_enum(source.m_bd, source.nelec(), source.m_ladder->m_nboson_max);
    conn::FrmBosOnv conn(source.m_bd);

    while (bra_enum.next(bra, ibra)) {
        size_t iket = ~0ul;
        enums::FermiBosOnv ket_enum(source.m_bd, source.nelec(), source.m_ladder->m_nboson_max);
        while (ket_enum.next(ket, iket)) {
            conn.connect(bra, ket);
            auto h_elem = source.get_element(bra, conn);
            (*this)(ibra, iket) = h_elem;
        }
    }
}

void DenseHamiltonian::setup_bos(const Hamiltonian &source) {

    /*
    buffered::BosOnv bra(source.m_bd);
    buffered::BosOnv ket(source.m_bd);

    foreach::rtnd::Ordered<false, true> foreach_bra(source.m_bd.m_nmode, source.nboson());
    foreach::rtnd::Ordered<false, true> foreach_ket(source.m_bd.m_nmode, source.nboson());

    auto ket_fn = [&]() {
        foreach_ket.inds()
    };
     */
}

DenseHamiltonian::DenseHamiltonian(const Hamiltonian &source) : SquareMatrix<defs::ham_t>(source.nci()) {
    if (source.m_nboson_max) setup_frmbos(source);
    else setup_frm(source);
}
