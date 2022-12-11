//
// Created by anderson on 19/07/2022.
//

#include "HubbardBosHam.h"
#include "M7_lib/excitgen/bos/BosHubbardUniform.h"

HubbardBosHam::HubbardBosHam(ham_comp_t u, const std::shared_ptr<lattice::Lattice> &lattice) :
        BosHam(lattice), m_u(u){
    m_contribs_0011.set_nonzero(opsig::c_0011);
    m_contribs_0022.set_nonzero(opsig::c_0000);
    logging::info("Bose-Hubbard Hamiltonian initialized with U={}; {}", m_u, m_basis.m_lattice->m_info);
}

HubbardBosHam::HubbardBosHam(init_opts_t opts) :
    HubbardBosHam(opts.m_ham.m_hubbard.m_repulsion, lattice::make(opts.m_ham.m_hubbard)){}

ham_t HubbardBosHam::get_coeff_0011(uint_t a, uint_t i) const {
    // hopping coeff is always -t
    return -m_basis.m_lattice->m_sparse_inv.get(a, i);
}

ham_t HubbardBosHam::get_coeff_0022(uint_t a, uint_t b, uint_t i, uint_t j) const {
    return (a==i && b==j) ? m_u : 0.0;
}

ham_t HubbardBosHam::get_element_0000(const field::BosOnv &onv) const {
    return onv.occ_npair()*m_u;
}

ham_t HubbardBosHam::get_element_0011(const field::BosOnv &onv, const conn::BosOnv &conn) const {
    DEBUG_ASSERT_EQ(conn.size(), 2ul, "incorrect connection exsig");
    const auto imode = conn.m_ann[0].m_imode;
    const auto jmode = conn.m_cre[0].m_imode;
    int t_mat_element = -m_basis.m_lattice->m_sparse_inv.get(imode, jmode);
    if (!t_mat_element) return 0.0;
    return onv.occ_fac(conn)*t_mat_element;
}

void HubbardBosHam::log_data() const {
    BosHam::log_data();
}

uint_t HubbardBosHam::default_nboson() const {
    return m_basis.m_nmode;
}

HamOpTerm::excit_gen_list_t HubbardBosHam::make_excit_gens(PRNG& prng, const conf::Propagator& /*opts*/) const {
    excit_gen_list_t list;
    list.emplace_front(new exgen::BosHubbardUniform(*this, prng));
    return list;
}

HamOpTerm::conn_foreach_list_t HubbardBosHam::make_foreach_iters() const {
    conn_foreach_list_t list;
    list.emplace_front(new conn_foreach::bos::Hubbard);
    return list;
}
