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
#include "src/core/io/Logging.h"

class AbInitioHamiltonian : public Hamiltonian {
    defs::ham_t m_int_0;
    Integrals_1e<defs::ham_t, defs::isym_1e> m_int_1;
    Integrals_2e<defs::ham_t, defs::isym_2e> m_int_2;

public:

    AbInitioHamiltonian(const FcidumpFileReader<defs::ham_t>& file_reader) :
        Hamiltonian(file_reader.nelec(), file_reader.nspatorb(), file_reader.spin_conserving()),
        m_int_1(file_reader.norb(), file_reader.spin_resolved()),
        m_int_2(file_reader.norb(), file_reader.spin_resolved()){
        defs::inds inds(4);
        defs::ham_t value;

        logger::write("Loading ab-initio Hamiltonian from FCIDUMP...");
        while (file_reader.next(inds, value)) {
            if (m_int_2.valid_inds(inds)) m_int_2.set_from_fcidump(inds, value);
            else if (m_int_1.valid_inds(inds)) m_int_1.set_from_fcidump(inds, value);
            else if (inds[0] == ((size_t) -1)) m_int_0 = value;
        }
        logger::write("FCIDUMP loading complete.");
        //m_nci = integer_utils::combinatorial(nsite() * 2, nelec());
    }

    explicit AbInitioHamiltonian(const std::string& fname, bool spin_major):
    AbInitioHamiltonian(FcidumpFileReader<defs::ham_t>(fname, spin_major)){}

    using Hamiltonian::get_element_0;
    defs::ham_t get_element_0(const defs::det_work &occs, const size_t &nocc) const override;

    using Hamiltonian::get_element_1;
    defs::ham_t get_element_1(const AntisymConnection &connection) const override;

    using Hamiltonian::get_element_2;
    defs::ham_t get_element_2(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const override;

    const Integrals_1e<defs::ham_t, defs::isym_1e> &int_1() const;

    const Integrals_2e<defs::ham_t, defs::isym_2e> &int_2() const;

};


#endif //M7_ABINITIOHAMILTONIAN_H
