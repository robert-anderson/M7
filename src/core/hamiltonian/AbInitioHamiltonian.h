//
// Created by Robert John Anderson on 2020-01-18.
//

#ifndef M7_ABINITIOHAMILTONIAN_H
#define M7_ABINITIOHAMILTONIAN_H


#include <cstddef>
#include <string>
#include "Hamiltonian.h"
#include "src/core/integrals/Integrals_1e.h"
#include "src/core/integrals/Integrals_2e.h"

class AbInitioHamiltonian : public Hamiltonian {
    FcidumpFileIterator<defs::ham_t> m_file_iterator;
    defs::ham_t m_int_0;
    Integrals_1e<defs::ham_t, defs::isym_1e> m_int_1;
    Integrals_2e<defs::ham_t, defs::isym_2e> m_int_2;

public:

    explicit AbInitioHamiltonian(const std::string &fname);

    using Hamiltonian::get_element_0;
    defs::ham_t get_element_0(const defs::det_work &occs, const size_t &nocc) const override;

    using Hamiltonian::get_element_1;
    defs::ham_t get_element_1(const AntisymConnection &connection) const override;

    using Hamiltonian::get_element_2;
    defs::ham_t get_element_2(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const override;

    size_t nelec() const override;

    bool spin_resolved() const override;

    bool spin_conserving() const override;

    const Integrals_1e<defs::ham_t, defs::isym_1e> &int_1() const;

    const Integrals_2e<defs::ham_t, defs::isym_2e> &int_2() const;
};


#endif //M7_ABINITIOHAMILTONIAN_H
