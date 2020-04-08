//
// Created by Robert John Anderson on 2020-01-18.
//

#include <src/core/fermion/DecodedDeterminant.h>
#include <src/core/fermion/Connection.h>
#include "AbInitioHamiltonian.h"
#include "src/core/io/FcidumpFileIterator.h"

AbInitioHamiltonian::AbInitioHamiltonian(const std::string &fname) :
    Hamiltonian(FcidumpFileIterator<defs::ham_t>(fname).nsite()),
    m_file_iterator(fname),
    m_int_1(m_file_iterator.m_norb, spin_resolved()),
    m_int_2(m_file_iterator.m_norb, spin_resolved()) {
    defs::inds inds(4);
    defs::ham_t value;
    while (m_file_iterator.next(inds, value)) {
        if (m_int_2.valid_inds(inds)) m_int_2.set_from_fcidump(inds, value);
        else if (m_int_1.valid_inds(inds)) m_int_1.set_from_fcidump(inds, value);
        else if (inds[0] == ((size_t) -1)) m_int_0 = value;
    }
    m_nci = integer_utils::combinatorial(nsite() * 2, nelec());
}

size_t AbInitioHamiltonian::nelec() const { return m_file_iterator.m_nelec; }

bool AbInitioHamiltonian::spin_resolved() const { return m_file_iterator.m_spin_resolved; }

bool AbInitioHamiltonian::spin_conserving() const { return m_int_1.spin_conserving(); }

auto &AbInitioHamiltonian::int_1() const { return m_int_1; }

auto &AbInitioHamiltonian::int_2() const { return m_int_2; }

defs::ham_t AbInitioHamiltonian::get_element_0(const defs::inds &occs, const size_t &nocc) const {
    defs::ham_t element = m_int_0;
    for (size_t i=0ul; i<nocc; ++i) {
        auto const &occi = occs[i];
        element += m_int_1(occi, occi);
        for (size_t j=0ul; j<i; ++j) {
            auto const &occj = occs[j];
            element += m_int_2.phys_antisym_element(occi, occj, occi, occj);
        }
    }
    return element;
}

defs::ham_t AbInitioHamiltonian::get_element_1(const AntisymConnection &connection) const {
    const auto &cre = connection.cre(0);
    const auto &ann = connection.ann(0);
    const auto &coms = connection.com();
    const auto &ncom = connection.ncom();

    defs::ham_t element = m_int_1(cre, ann);
    for (size_t icom=0ul; icom<ncom; ++icom)
        element += m_int_2.phys_antisym_element(cre, coms[icom], ann, coms[icom]);
    return connection.phase()?-element:element;
}

defs::ham_t
AbInitioHamiltonian::get_element_2(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
    return m_int_2.phys_antisym_element(i,j,k,l);
}
