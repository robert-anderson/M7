//
// Created by Robert John Anderson on 2020-01-18.
//

#include "src/core/basis/Connection.h"
#include "src/core/io/Logging.h"
#include "AbInitioHamiltonian.h"

#if 0
const AbInitioHamiltonian::ints1_t &AbInitioHamiltonian::int_1() const { return m_int_1; }

const AbInitioHamiltonian::ints2_t &AbInitioHamiltonian::int_2() const { return m_int_2; }

defs::ham_t AbInitioHamiltonian::get_element_0(const defs::det_work &occs, const size_t &nocc) const {
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
#endif