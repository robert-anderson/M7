//
// Created by rja on 27/02/2020.
//

#ifndef M7_HAMILTONIAN_H
#define M7_HAMILTONIAN_H

#include <cstddef>
#include <src/core/basis/DeterminantConnection.h>
#include <src/core/basis/DecodedDeterminant.h>
#include <src/core/field/Elements.h>
#include "src/core/integrals/Integrals_1e.h"
#include "src/core/integrals/Integrals_2e.h"


class Hamiltonian {
protected:
    const size_t m_nelec;
    const size_t m_nsite;
    const bool m_spin_conserving_1e, m_spin_conserving_2e;
    const bool m_complex_valued;
    defs::ham_t m_int_0;
    typedef Integrals_1e<defs::ham_t, defs::isym_1e> ints1_t;
    typedef Integrals_2e<defs::ham_t, defs::isym_2e> ints2_t;
    ints1_t m_int_1;
    ints2_t m_int_2;

public:
    Hamiltonian(const size_t &nelec, const size_t &nsite, bool spin_conserving_1e, bool spin_conserving_2e,
                bool complex_valued, bool spin_resolved) :
            m_nelec(nelec), m_nsite(nsite),
            m_spin_conserving_1e(spin_conserving_1e),
            m_spin_conserving_2e(spin_conserving_2e),
            m_complex_valued(complex_valued),
            m_int_1(nsite, spin_resolved),
            m_int_2(nsite, spin_resolved)
            {}

    Hamiltonian(const FcidumpFileReader<defs::ham_t> &file_reader) :
            Hamiltonian(file_reader.nelec(), file_reader.nspatorb(),
                        file_reader.spin_conserving_1e(),
                        file_reader.spin_conserving_2e(),
                        file_reader.m_complex_valued, file_reader.spin_resolved()) {
        defs::inds inds(4);
        defs::ham_t value;

        logger::write("Loading ab-initio Hamiltonian from FCIDUMP...");
        while (file_reader.next(inds, value)) {
            if (ints2_t::valid_inds(inds)) m_int_2.set(inds, value);
            else if (ints1_t::valid_inds(inds)) m_int_1.set(inds, value);
            else if (inds[0] == ~0ul) m_int_0 = value;
        }
        mpi::barrier();
        logger::write("FCIDUMP loading complete.");
    }


    explicit Hamiltonian(const std::string& fname, bool spin_major):
            Hamiltonian(FcidumpFileReader<defs::ham_t>(fname, spin_major)){}

    consts::component_t<defs::ham_t>::type get_energy(const views::FermionOnv &det) const {
        return consts::real(get_element_0(det));
    }

    defs::ham_t get_element_0(const defs::det_work &occs, const size_t &nocc) const {
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

    defs::ham_t get_element_0(const OccupiedOrbitals &occs) const {
        return get_element_0(occs.m_inds, occs.m_nind);
    }

    defs::ham_t get_element_0(const views::FermionOnv &det) const {
        OccupiedOrbitals occs(det);
        return get_element_0(occs.m_inds, occs.m_nind);
    }

    defs::ham_t get_element_0(const AntisymConnection &connection) const {
        return get_element_0(connection.com(), connection.ncom());
    }

    defs::ham_t get_element_1(const AntisymConnection &connection) const {
        const auto &cre = connection.cre(0);
        const auto &ann = connection.ann(0);
        const auto &coms = connection.com();
        const auto &ncom = connection.ncom();

        defs::ham_t element = m_int_1(cre, ann);
        for (size_t icom = 0ul; icom < ncom; ++icom)
            element += m_int_2.phys_antisym_element(cre, coms[icom], ann, coms[icom]);
        return connection.phase() ? -element : element;
    }

    defs::ham_t get_element_2(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
        return m_int_2.phys_antisym_element(i, j, k, l);
    }

    defs::ham_t get_element_2(const DeterminantConnection &connection) const {
        return get_element_2(connection.cre(0), connection.cre(1), connection.ann(0), connection.ann(1));
    }

    defs::ham_t get_element_2(const AntisymConnection &connection) const {
        const auto element = get_element_2(connection.cre(0), connection.cre(1), connection.ann(0), connection.ann(1));
        return connection.phase() ? -element : element;
    }

    defs::ham_t get_element(const AntisymConnection &connection) const {
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

    defs::ham_t get_element(const views::FermionOnv &bra, const views::FermionOnv &ket) const {
        return get_element(AntisymConnection(ket, bra));
    }

    size_t nci() const {
        return integer_utils::combinatorial(2 * nsite(), nelec());
    }

    const size_t &nsite() const {
        return m_nsite;
    }

    bool spin_conserving_1e() const {
        return m_spin_conserving_1e;
    }

    bool spin_conserving_2e() const {
        return m_spin_conserving_2e;
    }

    bool spin_conserving() const {
        return m_spin_conserving_1e && m_spin_conserving_2e;
    }

    const size_t &nelec() const {
        return m_nelec;
    }

    const bool &complex_valued() const {
        return m_complex_valued;
    }

//    elements::FermionOnv guess_reference(const int &spin_level) const;
//
//    elements::FermionOnv refine_guess_reference(const views::FermionOnv &ref) const;
//
//    elements::FermionOnv choose_reference(const int &spin_level) const;
//
//    class DeterminantList : public MappedList<DeterminantElement> {
//    public:
//        DeterminantField determinant;
//
//        DeterminantList(std::string name, size_t nsite, size_t nbucket) :
//                MappedList(name, determinant, nbucket),
//                determinant(this, 1, nsite){}
//    };
//
//    void generate_ci_space(WalkerList* list, RankAllocator<DeterminantElement>& ra, const int &spin_level) const;

};

#endif //M7_HAMILTONIAN_H
