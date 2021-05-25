//
// Created by rja on 24/05/2021.
//

#include <src/core/io/FcidumpFileReader.h>
#include <src/core/basis/DecodedDeterminant.h>
#include <src/core/basis/Connections.h>
#include "HamiltonianParts.h"
#include "HamiltonianData.h"

buffered::FermionOnv ham_parts::Fermion::guess_reference(const int &spin_restrict) const {
    buffered::FermionOnv ref(m_nsite);
    ASSERT((size_t)abs(spin_restrict) % 2 == nelec() % 2);
    size_t n_spin_0 = (nelec() + spin_restrict) / 2;
    size_t n_spin_1 = nelec() - n_spin_0;
    for (size_t i = 0ul; i < n_spin_0; ++i) ref.set({0, i});
    for (size_t i = 0ul; i < n_spin_1; ++i) ref.set({1, i});
    ASSERT(ref.spin() == spin_restrict);
    return ref;
}

ham_parts::Fermion::Fermion(const size_t &nelec, const size_t &nsite, bool spin_conserving_1e,
                                       bool spin_conserving_2e, bool complex_valued, bool spin_resolved, size_t int_2e_rank) :
        m_nelec(nelec), m_nsite(nsite),
        m_spin_conserving_1e(spin_conserving_1e),
        m_spin_conserving_2e(spin_conserving_2e),
        m_complex_valued(complex_valued),
        m_int_2e_rank(int_2e_rank),
        m_int_1(nsite, spin_resolved),
        m_int_2(nsite, spin_resolved){}

ham_parts::Fermion::Fermion(const FcidumpFileReader &file_reader) :
        Fermion(file_reader.nelec(), file_reader.nspatorb(),
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

ham_parts::Fermion::Fermion(std::string fname, bool spin_major) :
        Fermion(FcidumpFileReader(fname, spin_major)){}

defs::ham_comp_t ham_parts::Fermion::get_energy(const fields::Onv<0> &det) const {
    return consts::real(get_element_0(det));
}

defs::ham_t ham_parts::Fermion::get_element_0(const defs::inds &occs, const size_t &nocc) const {
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

defs::ham_t ham_parts::Fermion::get_element_0(const OccupiedOrbitals &occs) const {
    return get_element_0(occs.inds(), occs.size());
}

defs::ham_t ham_parts::Fermion::get_element_0(const fields::Onv<0> &fonv) const {
    // TODO: make occs member data
    OccupiedOrbitals occs(fonv);
    return get_element_0(occs.inds(), occs.size());
}

defs::ham_t ham_parts::Fermion::get_element_0(const conn::Antisym<0> &connection) const {
    return get_element_0(connection.com(), connection.ncom());
}

defs::ham_t ham_parts::Fermion::get_element_1(const conn::Antisym<0> &connection) const {
    const auto &cre = connection.cre(0);
    const auto &ann = connection.ann(0);
    const auto &coms = connection.com();
    const auto &ncom = connection.ncom();

    defs::ham_t element = m_int_1(cre, ann);
    for (size_t icom = 0ul; icom < ncom; ++icom)
        element += m_int_2.phys_antisym_element(cre, coms[icom], ann, coms[icom]);
    return connection.phase() ? -element : element;
}

defs::ham_t ham_parts::Fermion::get_element_2(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
    return m_int_2.phys_antisym_element(i, j, k, l);
}

defs::ham_t ham_parts::Fermion::get_element_2(const conn::Basic<0> &connection) const {
    return get_element_2(connection.cre(0), connection.cre(1), connection.ann(0), connection.ann(1));
}

defs::ham_t ham_parts::Fermion::get_element_2(const conn::Antisym<0> &connection) const {
    const auto element = get_element_2(connection.cre(0), connection.cre(1), connection.ann(0), connection.ann(1));
    return connection.phase() ? -element : element;
}

defs::ham_t ham_parts::Fermion::get_element(const conn::Antisym<0> &connection) const {
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

defs::ham_t ham_parts::Fermion::get_element(const fields::Onv<0> &bra, const fields::Onv<0> &ket) const {
    return get_element(AntisymFermionOnvConnection(ket, bra));
}

defs::ham_t ham_parts::Boson::get_element_0(const fields::BosonOnv &onv) const {
    defs::ham_t res = 0;
    for (size_t imode = 0ul; imode < m_nmode; ++imode)
        res += m_omegas[imode] * onv[imode];
    return res;
}

defs::ham_comp_t ham_parts::Boson::get_energy(const fields::BosonOnv &onv) const {
    return consts::real(get_element_0(onv));
}

defs::ham_t ham_parts::Boson::get_element_0(const conn::Boson &conn) const {
    defs::ham_t res = 0;
    for (size_t imode = 0ul; imode < m_nmode; ++imode)
        res += m_omegas[imode] * conn.com(imode);
    return res;
}

defs::ham_t ham_parts::Boson::get_element_0(const conn::Antisym<1> &conn) const {
    return get_element_0(conn.m_bonvconn);
}

defs::ham_t ham_parts::Boson::get_element(const conn::Boson &bonvconn) const {
    switch (bonvconn.nchanged_mode()) {
        case 0:
            return get_element_0(bonvconn);
        default:
            return 0;
    }
}

ham_parts::Coupling::Coupling(size_t nmode, size_t nboson_cutoff, defs::ham_t v) :
        m_nmode(nmode), m_nboson_cutoff(nboson_cutoff), m_v(v){}

defs::ham_t ham_parts::Coupling::v(const size_t &p, const size_t &q, const size_t &n) const {
    return (p == q && p == n) ? m_v : 0.0;
}

defs::ham_t ham_parts::Coupling::get_element_1(const size_t &p, const size_t &imode, const size_t &com) const {
    const auto occ_fac = std::sqrt(com + 1);
    return v(p, p, imode) * occ_fac;
}

defs::ham_t ham_parts::Coupling::get_element_1(const conn::Antisym<1> &aconn) const {
    const auto imode = aconn.m_bonvconn.changed_mode(0);
    const auto change = aconn.m_bonvconn.changes(0);

    // bosons don't couple to higher fermion excitations (yet?)
    if (std::abs(change) > 1) return 0.0;

    switch (aconn.nexcit()) {
        case 0: {
            defs::ham_t res = 0;
            for (size_t iocc = 0ul; iocc < aconn.ncom(); ++iocc) {
                auto p = aconn.com()[iocc] % m_nmode;
                res += get_element_1(p, imode, aconn.m_bonvconn.com(imode));
            }
            return res;
        }
        case 1: {
            auto p = aconn.cre()[0] % m_nmode;
            auto q = aconn.ann()[0] % m_nmode;
            ASSERT(p != q) // spin conservation
            ASSERT(v(p, q, imode) == 0)
            return v(p, q, imode);
        }
        default:
            return 0;
    }
}

defs::ham_t ham_parts::Coupling::get_element(const conn::Antisym<1> &aconn) const {
    switch (aconn.m_bonvconn.nchanged_mode()) {
        case 1:
            return get_element_1(aconn);
        default:
            return 0;
    }
}
