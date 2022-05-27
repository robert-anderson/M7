//
// Created by Robert J. Anderson on 12/8/21.
//

#include "HubbardFrmHam.h"
#include "M7_lib/excitgen/frm/HubbardUniform.h"


HubbardFrmHam::HubbardFrmHam(defs::ham_t u, Lattice lattice) :
        FrmHam(lattice.nsite()),
        m_u(u), m_lattice(std::move(lattice)), m_format(m_lattice.m_spec.m_format),
        m_bcs(m_lattice.m_spec.m_bcs) {
    m_contribs_1100.set_nonzero(exsig_utils::ex_single);
    m_contribs_2200.set_nonzero(0);
    log::info("Hubbard Hamiltonian initialized with U={}; {}", m_u, m_lattice.info());
}

HubbardFrmHam::HubbardFrmHam(opt_pair_t opts) :
        HubbardFrmHam(opts.m_ham.m_hubbard.m_repulsion, lattice::make(opts.m_ham.m_hubbard)){}

defs::ham_t HubbardFrmHam::get_element_0000(const field::FrmOnv &onv) const {
    defs::ham_t h = 0.0;
    for (size_t isite = 0ul; isite < m_basis.m_nsite; ++isite)
        if (onv.get({0, isite}) && onv.get({1, isite})) h += m_u;
    return h;
}

defs::ham_t HubbardFrmHam::get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    DEBUG_ASSERT_EQ(conn.size(), 2ul, "incorrect connection exsig");
    auto isite = m_basis.isite(conn.m_ann[0]);
    auto jsite = m_basis.isite(conn.m_cre[0]);
    if (m_basis.ispin(conn.m_ann[0])!=m_basis.ispin(conn.m_cre[0])) return 0.0;
    int t_mat_element = -m_lattice.m_dense(isite, jsite);
    if (!t_mat_element) return 0.0;
    return conn.phase(onv) ? -t_mat_element : t_mat_element;
}

defs::ham_t HubbardFrmHam::get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    return 0;
}

void HubbardFrmHam::log_data() const {
    FrmHam::log_data();
}

defs::ham_t HubbardFrmHam::get_coeff_1100(size_t a, size_t i) const {
    // hopping coeff is always -t
    return -m_lattice.m_dense(a, i);
}

defs::ham_t HubbardFrmHam::get_coeff_2200(size_t a, size_t b, size_t i, size_t j) const {
    return 0.0;
}

HamOpTerm::excit_gen_list_t HubbardFrmHam::make_excit_gens(PRNG& prng, const conf::Propagator& opts) const {
    excit_gen_list_t list;
    list.emplace_front(new HubbardUniform(*this, prng));
    return list;
}

HamOpTerm::conn_foreach_list_t HubbardFrmHam::make_foreach_iters() const {
    conn_foreach_list_t list;
    list.emplace_front(new conn_foreach::frm::Hubbard(m_lattice));
    return list;
}
