//
// Created by rja on 27/02/2020.
//

#include "src/core/enumerator/ContainerCombinationEnumerator.h"
#include "src/core/basis/DecodedDeterminant.h"
#include "src/core/parallel/RankAllocator.h"
#include "FermionHamiltonian.h"


buffered::FermionOnv FermionHamiltonian::guess_reference(const int &spin_restrict) const {
    buffered::FermionOnv ref(m_nsite);
    ASSERT((size_t)abs(spin_restrict) % 2 == nelec() % 2);
    size_t n_spin_0 = (nelec() + spin_restrict) / 2;
    size_t n_spin_1 = nelec() - n_spin_0;
    for (size_t i = 0ul; i < n_spin_0; ++i) ref.set({0, i});
    for (size_t i = 0ul; i < n_spin_1; ++i) ref.set({1, i});
    ASSERT(ref.spin() == spin_restrict);
    return ref;
}
#if 0

FermionOnv FermionHamiltonian::refine_guess_reference(const DeterminantElement &ref) const {

    auto e_ref = get_energy(ref);
    /*
    * check that none of the single and double connections have a lower energy
    */

    OccupiedOrbitals occs(ref);
    VacantOrbitals vacs(ref);

    FermionOnv excited(m_nsite);

    for (auto occ : occs.m_inds) {
        for (auto vac : vacs.m_inds) {
            excited = ref;
            excited.excite(vac, occ);
            if (get_energy(excited) < e_ref) return refine_guess_reference(excited);
        }
    }

    ContainerCombinationEnumerator<defs::det_work> occ_enumerator(occs.m_inds, 2);
    defs::inds occ_inds(2);
    while (occ_enumerator.next(occ_inds)) {
        {
            ContainerCombinationEnumerator<defs::det_work> vac_enumerator(vacs.m_inds, 2);
            defs::inds vac_inds(2);
            while (vac_enumerator.next(vac_inds)) {
                excited = ref;
                excited.excite(vac_inds[0], vac_inds[1], occ_inds[0], occ_inds[1]);
                if (get_energy(excited) < e_ref) return refine_guess_reference(excited);
            }
        }
    }
    excited = ref;
    return excited;
}

FermionOnv FermionHamiltonian::choose_reference(const int &spin_level) const {
    auto ref = guess_reference(spin_level);
    ref = refine_guess_reference(ref);
    return ref;
}

void
FermionHamiltonian::generate_ci_space(WalkerTable *list, RankAllocator<DeterminantElement> &ra, const int &spin_level) const {
    size_t nalpha = nelec() / 2 + spin_level;
    size_t nbeta = nelec() / 2 - spin_level;

    defs::inds alpha_sites(nsite(), 0);
    defs::inds beta_sites(nsite(), 0);

    for (size_t i = 0; i < nsite(); ++i) {
        alpha_sites[i] = i;
        beta_sites[i] = nsite() + i;
    }
    ContainerCombinationEnumerator<defs::inds> alpha_enumerator(alpha_sites, nsite(), nalpha);
    defs::inds alpha_inds(nalpha);
    defs::inds beta_inds(nbeta);

    FermionOnv work_det(nsite());
    while (alpha_enumerator.next(alpha_inds)) {
        ContainerCombinationEnumerator<defs::inds> beta_enumerator(beta_sites, nsite(), nbeta);
        while (beta_enumerator.next(beta_inds)) {
            work_det.zero();
            work_det.set(alpha_inds);
            work_det.set(beta_inds);
            if (!mpi::i_am(ra.get_rank(work_det))) continue;
            auto irow = list->expand_push(work_det);
            auto det = list->m_determinant(irow);
            auto h_diag = list->m_hdiag(irow);
            det.set(alpha_inds);
            det.set(beta_inds);
            h_diag = get_energy(det);
        }
    }
#ifndef NDEBUG
    if (mpi::nrank()==1) {
        ASSERT(list->high_water_mark(0) ==
               integer_utils::combinatorial(nsite(), nalpha) * integer_utils::combinatorial(nsite(), nbeta))
    }
#endif
}

#endif

FermionHamiltonian::FermionHamiltonian(const size_t &nelec, const size_t &nsite, bool spin_conserving_1e,
                                       bool spin_conserving_2e, bool complex_valued, bool spin_resolved, size_t int_2e_rank) :
        m_nelec(nelec), m_nsite(nsite),
        m_spin_conserving_1e(spin_conserving_1e),
        m_spin_conserving_2e(spin_conserving_2e),
        m_complex_valued(complex_valued),
        m_int_2e_rank(int_2e_rank),
        m_int_1(nsite, spin_resolved),
        m_int_2(nsite, spin_resolved),
        m_terms(new Terms(*this))
{}

FermionHamiltonian::FermionHamiltonian(const FcidumpFileReader &file_reader) :
        FermionHamiltonian(file_reader.nelec(), file_reader.nspatorb(),
                           file_reader.spin_conserving_1e(),
                           file_reader.spin_conserving_2e(),
                           file_reader.m_complex_valued,
                           file_reader.spin_resolved(),
                           file_reader.int_2e_rank()) {

    using namespace ham_data;
    defs::inds inds(4);
    defs::ham_t value;

    // assume no contribution cases to be nonzero unless a counterexample is found
    m_nonzero_contribs.fill(false);

    log::info("Reading fermion Hamiltonian coefficients from FCIDUMP file \"" + file_reader.m_fname + "\"...");
    while (file_reader.next(inds, value)) {
        auto contrib_case = ham_data::get_coupling_contrib_case(inds);
        auto rank = c_contrib_inds[contrib_case].m_term.m_nann;
        // all integrals are particle number conserving (for now)
        ASSERT(rank == c_contrib_inds[contrib_case].m_term.m_ncre);
        m_nonzero_contribs[contrib_case] = true;
        if (rank==2) {
            if (m_on_site_only_0022 && !on_site(inds[0], inds[1], inds[2], inds[3]))
                m_on_site_only_0022 = false;
            m_int_2.set(inds, value);
        }
        else if (rank==1) {
            if (m_nn_only_1111 && !nearest_neighbors(inds[0], inds[1], false)) m_nn_only_1111 = false;
            if (m_nnp_only_1111 && !nearest_neighbors(inds[0], inds[1], true)) m_nnp_only_1111 = false;
            m_int_1.set(inds, value);
        }
        else if (rank==0) m_int_0 = value;
        else MPI_ABORT("File reader error");
    }
    mpi::barrier();
    log::info("FCIDUMP loading complete.");

    // some useful output to identify the kind of H detected
    if (!m_nonzero_contribs[contrib_0011])
        log::info("1-electron term in Hamiltonian has no diagonal contributions");
    if (!m_nonzero_contribs[contrib_1111])
        log::info("1-electron term in Hamiltonian has no single-excitation contributions");
    if (!m_nonzero_contribs[contrib_0022])
        log::info("2-electron term in Hamiltonian has no diagonal contributions");
    if (!m_nonzero_contribs[contrib_1122])
        log::info("2-electron term in Hamiltonian has no single-excitation contributions");
    if (!m_nonzero_contribs[contrib_2222])
        log::info("2-electron term in Hamiltonian has no double-excitation contributions");
    if (m_nnp_only_1111)
        log::info("single-excitation contributions to 1-electron term are periodic nearest-neighbor only");
    else if (m_nn_only_1111)
        log::info("single-excitation contributions to 1-electron term are nearest-neighbor only");
    if (m_on_site_only_0022)
        log::info("2-electron term diagonal contributions are on-site only");
}

FermionHamiltonian::FermionHamiltonian(std::string fname, bool spin_major) :
        FermionHamiltonian(FcidumpFileReader(fname, spin_major)){}

defs::ham_comp_t FermionHamiltonian::get_energy(const fields::Onv<0> &det) const {
    return consts::real(get_element_0(det));
}

defs::ham_t FermionHamiltonian::get_element_0(const defs::inds &occs, const size_t &nocc) const {
    defs::ham_t element = m_int_0;
    for (size_t i = 0ul; i < nocc; ++i) {
        auto const &occi = occs[i];
        element += m_int_1(occi, occi);
        for (size_t j = 0ul; j < i; ++j) {
            auto const &occj = occs[j];
            element += m_int_2.phys_antisym_element(occi, occj, occi, occj);
        }
    }
    return element;
}

defs::ham_t FermionHamiltonian::get_element_0(const OccupiedOrbitals &occs) const {
    return get_element_0(occs.inds(), occs.size());
}

defs::ham_t FermionHamiltonian::get_element_0(const fields::Onv<0> &fonv) const {
    // TODO: make occs member data
    OccupiedOrbitals occs(fonv);
    return get_element_0(occs.inds(), occs.size());
}

defs::ham_t FermionHamiltonian::get_element_0(const conn::Antisym<0> &connection) const {
    return get_element_0(connection.com(), connection.ncom());
}

defs::ham_t FermionHamiltonian::get_element_1(const conn::Antisym<0> &connection) const {
    const auto &cre = connection.cre(0);
    const auto &ann = connection.ann(0);
    const auto &coms = connection.com();
    const auto &ncom = connection.ncom();

    defs::ham_t element = m_int_1(cre, ann);
    for (size_t icom = 0ul; icom < ncom; ++icom)
        element += m_int_2.phys_antisym_element(cre, coms[icom], ann, coms[icom]);
    return connection.phase() ? -element : element;
}

defs::ham_t FermionHamiltonian::get_element_2(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
    return m_int_2.phys_antisym_element(i, j, k, l);
}

defs::ham_t FermionHamiltonian::get_element_2(const conn::Basic<0> &connection) const {
    return get_element_2(connection.cre(0), connection.cre(1), connection.ann(0), connection.ann(1));
}

defs::ham_t FermionHamiltonian::get_element_2(const conn::Antisym<0> &connection) const {
    const auto element = get_element_2(connection.cre(0), connection.cre(1), connection.ann(0), connection.ann(1));
    return connection.phase() ? -element : element;
}

defs::ham_t FermionHamiltonian::get_element(const conn::Antisym<0> &connection) const {
    switch (connection.nexcit()) {
        case 0:
            return get_element_0(connection);
        case 1: ASSERT(connection.ncom() + connection.nexcit() == nelec());
            return get_element_1(connection);
        case 2:
            return get_element_2(connection);
        default:
            return 0;
    }
}

defs::ham_t FermionHamiltonian::get_element(const fields::Onv<0> &bra, const fields::Onv<0> &ket) const {
    return get_element(AntisymFermionOnvConnection(ket, bra));
}

defs::ham_t FermionHamiltonian::Terms::get_element_00(const defs::inds &occs, const size_t &nocc) const {
    defs::ham_t element = m_ham.m_int_0;
    for (size_t i = 0ul; i < nocc; ++i) {
        auto const &occi = occs[i];
        element += m_ham.m_int_1(occi, occi);
        for (size_t j = 0ul; j < i; ++j) {
            auto const &occj = occs[j];
            element += m_ham.m_int_2.phys_antisym_element(occi, occj, occi, occj);
        }
    }
    return element;
}

defs::ham_t FermionHamiltonian::Terms::get_element_11(const conn::Antisym<0> &conn) const {
    ASSERT(conn.rank_label() == ham_data::encode_rank_label(1, 1));
    const auto &cre = conn.cre(0);
    const auto &ann = conn.ann(0);
    const auto &coms = conn.com();
    const auto &ncom = conn.ncom();

    defs::ham_t element = m_ham.m_int_1(cre, ann);
    for (size_t icom = 0ul; icom < ncom; ++icom)
        element += m_ham.m_int_2.phys_antisym_element(cre, coms[icom], ann, coms[icom]);
    return conn.phase() ? -element : element;
}

defs::ham_t FermionHamiltonian::Terms::get_element_22(const conn::Antisym<0> &conn) const {
    ASSERT(conn.rank_label() == ham_data::encode_rank_label(2, 2));
    const auto element = m_ham.m_int_2.phys_antisym_element(
            conn.cre(0), conn.cre(1), conn.ann(0), conn.ann(1));
    return conn.phase() ? -element : element;
}

defs::ham_t FermionHamiltonian::Terms::get_element_01(const conn::Antisym<0> &conn) const {
    return 0.0;
}

defs::ham_t FermionHamiltonian::Terms::get_element_12(const conn::Antisym<0> &conn) const {
    return 0.0;
}

defs::ham_t FermionHamiltonian::Terms::get_element_02(const conn::Antisym<0> &conn) const {
    return 0.0;
}

defs::ham_t FermionHamiltonian::Terms::get_element_10(const conn::Antisym<0> &conn) const {
    return 0.0;
}

defs::ham_t FermionHamiltonian::Terms::get_element_21(const conn::Antisym<0> &conn) const {
    return 0.0;
}

defs::ham_t FermionHamiltonian::Terms::get_element_20(const conn::Antisym<0> &conn) const {
    return 0.0;
}

defs::ham_t FermionHamiltonian::Terms::get_element(const conn::Antisym<0> &conn) const {
    auto rank_label = conn.rank_label();
    // should compile to jump table
    switch (rank_label) {
        case ham_data::encode_rank_label(0, 0):
            return get_element_00(conn.com(), conn.ncom());
        case ham_data::encode_rank_label(1, 1):
            return get_element_11(conn);
        case ham_data::encode_rank_label(2, 2):
            return get_element_22(conn);
        case ham_data::encode_rank_label(0, 1):
            return get_element_01(conn);
        case ham_data::encode_rank_label(1, 2):
            return get_element_12(conn);
        case ham_data::encode_rank_label(0, 2):
            return get_element_02(conn);
        case ham_data::encode_rank_label(1, 0):
            return get_element_10(conn);
        case ham_data::encode_rank_label(2, 1):
            return get_element_21(conn);
        case ham_data::encode_rank_label(2, 0):
            return get_element_20(conn);
        default:
            return 0.0;
    }
}

bool FermionHamiltonian::Terms::update_helement(bool get_h, bool h_nonzero_only) const {
    if (get_h || h_nonzero_only) m_helement_work = get_element(m_conn_work);
    else m_helement_work = 0.0;
    if (h_nonzero_only) return !consts::float_is_zero(m_helement_work);
    return true;
}

void FermionHamiltonian::Terms::perform_single(const fields::FermionOnv &src_onv, const size_t &occ, const size_t &vac,
                                               const FermionHamiltonian::Terms::body_fn_t &body, bool get_h,
                                               bool h_nonzero_only) const {
    m_conn_work.zero();
    m_conn_work.add(occ, vac);
    m_conn_work.apply(src_onv, m_onv_work);
    if (update_helement(get_h, h_nonzero_only)) body(m_conn_work, m_onv_work, m_helement_work);
}

void
FermionHamiltonian::Terms::perform_double(const fields::FermionOnv &src_onv, const size_t &occ1, const size_t &occ2,
                                          const size_t &vac1, const size_t &vac2,
                                          const FermionHamiltonian::Terms::body_fn_t &body, bool get_h,
                                          bool h_nonzero_only) const {
    m_conn_work.zero();
    m_conn_work.add(occ1, occ2, vac1, vac2);
    m_conn_work.apply(src_onv, m_onv_work);
    if (update_helement(get_h, h_nonzero_only)) body(m_conn_work, m_onv_work, m_helement_work);
}

void FermionHamiltonian::Terms::foreach_connection_singles(const fields::FermionOnv &src_onv, const defs::inds &occs,
                                                           const defs::inds &vacs,
                                                           const FermionHamiltonian::Terms::body_fn_t &body, bool get_h,
                                                           bool h_nonzero_only) const {
    for (size_t iocc = 0ul; iocc < occs.size(); ++iocc) {
        auto &occ = occs[iocc];
        for (size_t ivac = 0ul; ivac < vacs.size(); ++ivac) {
            auto &vac = vacs[ivac];
            perform_single(src_onv, occ, vac, body, get_h, h_nonzero_only);
        }
    }
}

void
FermionHamiltonian::Terms::foreach_connection_subset_doubles(const fields::FermionOnv &src_onv, const defs::inds &occs1,
                                                             const defs::inds &occs2, const defs::inds &vacs1,
                                                             const defs::inds &vacs2,
                                                             const FermionHamiltonian::Terms::body_fn_t &body,
                                                             bool get_h, bool h_nonzero_only) const {
    for (size_t iocc1 = 0ul; iocc1 < occs1.size(); ++iocc1) {
        auto &occ1 = occs1[iocc1];
        for (size_t ivac1 = 0ul; ivac1 < vacs1.size(); ++ivac1) {
            auto &vac1 = vacs1[ivac1];
            size_t iocc2_start = (&occs1==&occs2) ? iocc1+1 : 0ul;
            for (size_t iocc2 = iocc2_start; iocc2 < occs2.size(); ++iocc2) {
                auto &occ2 = occs2[iocc2];
                size_t ivac2_start = (&vacs1==&vacs2) ? ivac1+1 : 0ul;
                for (size_t ivac2 = ivac2_start; ivac2 < vacs2.size(); ++ivac2) {
                    auto &vac2 = vacs2[ivac2];
                    perform_double(src_onv, occ1, occ2, vac1, vac2, body, get_h, h_nonzero_only);
                }
            }
        }
    }
}

void
FermionHamiltonian::Terms::foreach_connection_subset_doubles(const fields::FermionOnv &src_onv, const defs::inds &occs,
                                                             const defs::inds &vacs,
                                                             const FermionHamiltonian::Terms::body_fn_t &body,
                                                             bool get_h, bool h_nonzero_only) const {
    foreach_connection_subset_doubles(src_onv, occs, vacs, occs, vacs, body, get_h, h_nonzero_only);
}

void FermionHamiltonian::Terms::foreach_connection_subset(const fields::FermionOnv &src_onv, const defs::inds &occs1,
                                                          const defs::inds &occs2, const defs::inds &vacs1,
                                                          const defs::inds &vacs2,
                                                          const FermionHamiltonian::Terms::body_fn_t &body, bool get_h,
                                                          bool h_nonzero_only) const {
    for (size_t iocc1 = 0ul; iocc1 < occs1.size(); ++iocc1) {
        auto &occ1 = occs1[iocc1];
        for (size_t ivac1 = 0ul; ivac1 < vacs1.size(); ++ivac1) {
            auto &vac1 = vacs1[ivac1];
            perform_single(src_onv, occ1, vac1, body, get_h, h_nonzero_only);
            size_t iocc2_start = (&occs1==&occs2) ? iocc1+1 : 0ul;
            for (size_t iocc2 = iocc2_start; iocc2 < occs2.size(); ++iocc2) {
                auto &occ2 = occs2[iocc2];
                size_t ivac2_start = (&vacs1==&vacs2) ? ivac1+1 : 0ul;
                for (size_t ivac2 = ivac2_start; ivac2 < vacs2.size(); ++ivac2) {
                    auto &vac2 = vacs2[ivac2];
                    perform_double(src_onv, occ1, occ2, vac1, vac2, body, get_h, h_nonzero_only);
                }
            }
        }
    }
}

void FermionHamiltonian::Terms::foreach_connection_subset(const fields::FermionOnv &src_onv, const defs::inds &occs,
                                                          const defs::inds &vacs,
                                                          const FermionHamiltonian::Terms::body_fn_t &body, bool get_h,
                                                          bool h_nonzero_only) const {
    foreach_connection_subset(src_onv, occs, occs, vacs, vacs, body, get_h, h_nonzero_only);
}

void FermionHamiltonian::Terms::foreach_connection(const fields::FermionOnv &src_onv,
                                                   const FermionHamiltonian::Terms::body_fn_t &body, bool get_h,
                                                   bool h_nonzero_only) const {
    m_helement_work = 0.0;
    m_occ_work.update(src_onv);
    m_vac_work.update(src_onv);
    foreach_connection_subset(src_onv, m_occ_work.inds(), m_vac_work.inds(), body, get_h, h_nonzero_only);
}

FermionHamiltonian::SpinTerms::SpinTerms(const FermionHamiltonian &ham) : Terms(ham),
                                                                          m_spin_occ_work(ham.nsite()), m_spin_vac_work(ham.nsite()){
}

void FermionHamiltonian::SpinTerms::foreach_connection(const fields::FermionOnv &src_onv,
                                                       const FermionHamiltonian::Terms::body_fn_t &body, bool get_h,
                                                       bool h_nonzero_only) const {
    m_helement_work = 0.0;
    m_spin_occ_work.update(src_onv);
    m_spin_vac_work.update(src_onv);
    // spin a->a, aa->aa
    foreach_connection_subset(src_onv, m_spin_occ_work[0], m_spin_vac_work[0], body, get_h, h_nonzero_only);
    // spin b->b, bb->bb
    foreach_connection_subset(src_onv, m_spin_occ_work[1], m_spin_vac_work[1], body, get_h, h_nonzero_only);
    // spin ab->ab
    foreach_connection_subset_doubles(src_onv, m_spin_occ_work[0], m_spin_occ_work[1],
                                      m_spin_vac_work[0], m_spin_vac_work[1], body, get_h, h_nonzero_only);
}

FermionHamiltonian::Hubbard1DTerms::Hubbard1DTerms(const FermionHamiltonian &ham) : SpinTerms(ham){}

void FermionHamiltonian::Hubbard1DTerms::foreach_connection(const fields::FermionOnv &src_onv,
                                                            const FermionHamiltonian::Terms::body_fn_t &body,
                                                            bool get_h, bool h_nonzero_only) const {
    m_helement_work = 0.0;
    m_spin_occ_work.update(src_onv);
    // spin a
    for (auto& occ: m_spin_occ_work[0]) {
        // to the left
        if (occ > 0 && !src_onv.get(occ - 1))
            perform_single(src_onv, occ, occ - 1, body, get_h, h_nonzero_only);
        // to the right
        if (occ+1<src_onv.m_nsite && !src_onv.get(occ+1))
            perform_single(src_onv, occ, occ+1, body, get_h, h_nonzero_only);
    }
    // spin b
    for (auto& occ: m_spin_occ_work[1]) {
        // to the left
        if (occ > src_onv.m_nsite && !src_onv.get(occ - 1))
            perform_single(src_onv, occ, occ - 1, body, get_h, h_nonzero_only);
        // to the right
        if (occ+1<src_onv.nbit() && !src_onv.get(occ+1))
            perform_single(src_onv, occ, occ+1, body, get_h, h_nonzero_only);
    }
}

void FermionHamiltonian::Hubbard1DPbcTerms::foreach_connection(const fields::FermionOnv &src_onv,
                                                               const FermionHamiltonian::Terms::body_fn_t &body,
                                                               bool get_h, bool h_nonzero_only) const {
    m_helement_work = 0.0;
    m_spin_occ_work.update(src_onv);
    // spin a
    for (auto& occ: m_spin_occ_work[0]) {
        // to the left
        if (occ > 0 && !src_onv.get(occ - 1))
            perform_single(src_onv, occ, occ - 1, body, get_h, h_nonzero_only);
            // to the left (wrap-around)
        else if (occ == 0 && !src_onv.get(src_onv.m_nsite))
            perform_single(src_onv, occ, src_onv.m_nsite-1, body, get_h, h_nonzero_only);
        // to the right
        if (occ+1<src_onv.m_nsite && !src_onv.get(occ+1))
            perform_single(src_onv, occ, occ+1, body, get_h, h_nonzero_only);
            // to the right (wrap-around)
        else if (occ+1 == src_onv.m_nsite && !src_onv.get(0))
            perform_single(src_onv, occ, 0, body, get_h, h_nonzero_only);
    }
    // spin b
    for (auto& occ: m_spin_occ_work[1]) {
        // to the left
        if (occ > src_onv.m_nsite && !src_onv.get(occ - 1))
            perform_single(src_onv, occ, occ - 1, body, get_h, h_nonzero_only);
            // to the left (wrap-around)
        else if (occ == src_onv.m_nsite && !src_onv.get(src_onv.nbit()-1))
            perform_single(src_onv, occ, src_onv.nbit()-1, body, get_h, h_nonzero_only);
        // to the right
        if (occ+1<src_onv.nbit() && !src_onv.get(occ+1))
            perform_single(src_onv, occ, occ+1, body, get_h, h_nonzero_only);
            // to the right (wrap-around)
        else if (occ+1 == src_onv.nbit() && !src_onv.get(src_onv.m_nsite-1))
            perform_single(src_onv, occ, src_onv.m_nsite-1, body, get_h, h_nonzero_only);
    }
}
