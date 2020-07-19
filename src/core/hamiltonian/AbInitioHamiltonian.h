//
// Created by Robert John Anderson on 2020-01-18.
//

#ifndef M7_ABINITIOHAMILTONIAN_H
#define M7_ABINITIOHAMILTONIAN_H


#include <cstddef>
#include <string>
#include <src/core/io/ProgressBar.h>
#include "Hamiltonian.h"
#include "src/core/integrals/Integrals_1e.h"
#include "src/core/integrals/Integrals_2e.h"
#include "src/core/io/Logging.h"

class AbInitioHamiltonian : public Hamiltonian {
    defs::ham_t m_int_0;
    typedef Integrals_1e<defs::ham_t, defs::isym_1e> ints1_t;
    typedef Integrals_2e<defs::ham_t, defs::isym_2e> ints2_t;
    ints1_t m_int_1;
    ints2_t m_int_2;

public:

    AbInitioHamiltonian(const FcidumpFileReader<defs::ham_t>& file_reader) :
        Hamiltonian(file_reader.nelec(), file_reader.nspatorb(), file_reader.spin_conserving()),
        m_int_1(file_reader.norb(), file_reader.spin_resolved()),
        m_int_2(file_reader.norb(), file_reader.spin_resolved()){
        defs::inds inds(4);
        defs::ham_t value;

        logger::write("Loading ab-initio Hamiltonian from FCIDUMP...");
        ProgressBar progress_bar(file_reader.nline(), 100);
        while (file_reader.next(inds, value)) {
            if (ints2_t::valid_inds(inds)) m_int_2.set(inds, value);
            else if (ints1_t::valid_inds(inds)) m_int_1.set(inds, value);
            else if (inds[0] == ~0ul) m_int_0 = value;
            ++progress_bar;
            progress_bar.display();
        }
        progress_bar.done();
        mpi::barrier();
        logger::write("FCIDUMP loading complete.");
    }

    explicit AbInitioHamiltonian(const std::string& fname, bool spin_major):
    AbInitioHamiltonian(FcidumpFileReader<defs::ham_t>(fname, spin_major)){}

    using Hamiltonian::get_element_0;
    defs::ham_t get_element_0(const defs::det_work &occs, const size_t &nocc) const override;

    using Hamiltonian::get_element_1;
    defs::ham_t get_element_1(const AntisymConnection &connection) const override;

    using Hamiltonian::get_element_2;
    defs::ham_t get_element_2(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const override;

    const ints1_t &int_1() const;

    const ints2_t &int_2() const;

};


#endif //M7_ABINITIOHAMILTONIAN_H
